 ! 
 !  @BEGIN LICENSE
 ! 
 !  Hilbert: a space for quantum chemistry plugins to Psi4 
 ! 
 !  Copyright (c) 2020 by its authors (LICENSE).
 ! 
 !  The copyrights for code used from other parties are included in
 !  the corresponding files.
 ! 
 !  This program is free software: you can redistribute it and/or modify
 !  it under the terms of the GNU Lesser General Public License as published by
 !  the Free Software Foundation, either version 3 of the License, or
 !  (at your option) any later version.
 ! 
 !  This program is distributed in the hope that it will be useful,
 !  but WITHOUT ANY WARRANTY; without even the implied warranty of
 !  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 !  GNU Lesser General Public License for more details.
 ! 
 !  You should have received a copy of the GNU Lesser General Public License
 !  along with this program.  If not, see http://www.gnu.org/licenses/.
 ! 
 !  @END LICENSE
 ! 

module focas_hessian

  use focas_data

  implicit none

  contains

    subroutine diagonal_hessian(q,z,int2,den1,den2)
      implicit none
      real(wp), intent(in) :: int2(:),den1(:),den2(:)
      real(wp), intent(in) :: q(:,:),z(:,:)
      integer :: error
      real(wp) :: t0(2),t1(2),t_diag,t_full,t_precond

      min_diag_hessian_           = 1.0e-2_wp

      num_negative_diagonal_hessian_ = 0 

      ! active-doubly occupied pairs

      if ( rot_pair_%n_ad > 0 ) error=diagonal_hessian_ad(fock_i_%occ,fock_a_%occ,q,z,int2,den1,den2)

      ! active-active pairs

      if ( rot_pair_%n_aa > 0 ) then

        ! calculate diagonal Hesian
        t0 = timer()
        error=diagonal_hessian_aa(fock_i_%occ,fock_a_%occ,q,z,int2,den1,den2)
        t1 = timer()
        t_diag = t1(1) - t0(1)

!        ! calculate entire aa block of Hessian
!        t0 = timer()
!        error=full_hessian_aa(fock_i_%occ,q,z,int2,den1,den2)
!        t1 = timer()
!        t_full = t1(1) - t0(1)
!
!        ! calculate -Hinv*g (assumes orbital gradient is availabke at this stage)
!        t0 = timer()
!        error = full_hessian_preconditioner(full_orbital_hessian_%aa,        &
!                                        & full_orbital_hessian_%HinvXg_aa)
!        t1 = timer()
!        t_precond = t1(1) - t0(1)
!
!        write(*,'(3(a,1x,f6.3,3x))')'t_diag=',t_diag,'t_full=',t_full,'t_recond=',t_precond

      end if

      ! external-doubly occupied pairs

      if ( rot_pair_%n_ed > 0 ) error=diagonal_hessian_ed(fock_i_%occ,fock_i_%ext,fock_a_%occ,fock_a_%ext,int2)

      ! external-active pairs

      if ( rot_pair_%n_ea > 0 ) error=diagonal_hessian_ea(fock_i_%ext,fock_a_%ext,q,z,int2,den1,den2)

!      ! print the orbital hessian
!      call print_orbital_hessian()

      return
    end subroutine diagonal_hessian

    integer function diagonal_hessian_ad(f_i,f_a,q,z,int2,den1,den2)
      implicit none
      ! subroutine to compute the diagonal Hessian matrix element H(ti|ti) according to Eq. 4.7c of
      ! Chaban, Schmidt, Gordon, Theor. Chem. Acc., 97, 88-95 (1997) 
      ! H(ti|ti) = 2 * [ 2 * { ( f_i(t,t) + f_a(t,t) ) - ( f_i(i,i) + f_a(i,i) ) }  + d(t,t) * ( f_i(i,i) + f_a(i,i) ) - q(t,t) -z(t) ]
      ! where z(t,t) = sum_{u} [ d(t,u) * f_i(t,u) ] ] && a \in ext & t,u \in act
      real(wp), intent(in) :: q(:,:),z(:,:),den1(:),den2(:),int2(:)
      type(matrix_block), intent(in) :: f_i(:),f_a(:)

      integer :: grad_ind,t,i,tt_den,t_sym,t_i,i_i
      real(wp) :: h_val

      ! loop over all active - doubly-occupied pairs

      do t_sym = 1 , nirrep_

        grad_ind = rot_pair_%pair_offset(t_sym,rot_pair_%act_doc_type)

        do i = first_index_(t_sym,1) , last_index_(t_sym,1)

          do t = first_index_(t_sym,2) , last_index_(t_sym,2)

            ! look up integral/density indeces
            tt_den  = dens_%gemind(t,t)

            i_i     = trans_%class_to_irrep_map(i)
            t_i     = trans_%class_to_irrep_map(t)

            ! update gradient index
            grad_ind = grad_ind + 1
 
            ! compute Hessian value
            h_val    = 2.0_wp * (                                                                   &
                       2.0_wp * ( ( f_i(t_sym)%val(t_i,t_i) + f_a(t_sym)%val(t_i,t_i) )             &
                                - ( f_i(t_sym)%val(i_i,i_i) + f_a(t_sym)%val(i_i,i_i) ) )           &
                     + den1(tt_den) * f_i(t_sym)%val(i_i,i_i) - q(t-ndoc_tot_,t) - z(t-ndoc_tot_,t) )

            if ( use_exact_hessian_diagonal_ == 0 ) then

              h_val = h_val + 2.0_wp * den1(tt_den) * f_a(t_sym)%val(i_i,i_i)

            else

              if ( df_vars_%use_df_teints == 1 ) then

                h_val = h_val + te_terms_ad_df(i,t,t_sym,int2,den1,den2)

              else

                h_val = h_val + te_terms_ad(i,t,t_sym,int2,den1,den2)

              endif

            endif

            if ( h_val < 0.0_wp ) h_val = - h_val

            if ( h_val < min_diag_hessian_ ) h_val = min_diag_hessian_

            orbital_hessian_(grad_ind) = h_val

          end do

        end do

      end do

      diagonal_hessian_ad = 0

      return

    end function diagonal_hessian_ad

    function te_terms_ad(i,t,t_sym,int2,den1,den2)

      ! function to compute exact 2-e contribution for active-doubly occupied diagonal Hessian
      ! using 4-index integrals
      ! H(it|it) = 2 * sum{u,v} [ d(tt|uv) * g(ii|uv) + 2 * d(tv|tu) * g(ui|vi)         ]  
      !          + 4 * sum{u}   [ ( delta(t,u) - d(t,u) ) * ( 3 * g(ui|ti) - g(ii|tu) ) ]

      integer, intent(in)  :: i,t,t_sym
      real(wp), intent(in) :: int2(:),den1(:),den2(:)

      integer              :: tt_den,uv_den,tv_den,tu_den,uu_den
      integer              :: ii_int,uv_int,ui_int,vi_int,ti_int,tu_int,uu_int
      integer              :: u,v,u_sym,tu_sym,int_ind,den_ind,den_offset,int_offset
      real(wp)             :: te_terms_ad,dfac,val,int_val
      
      val = 0.0_wp
 
      tt_den = dens_%gemind(t,t)
      ii_int = ints_%gemind(i,i)
      ti_int = ints_%gemind(i,t)

      do u_sym = 1 , nirrep_
      
        tu_sym     = group_mult_tab_(t_sym,u_sym)
        den_offset = dens_%offset(tu_sym)
        int_offset = ints_%offset(tu_sym)

        do u = first_index_(u_sym,2) , last_index_(u_sym,2)

          tu_den = dens_%gemind(t,u)
          uu_den = dens_%gemind(u,u)
          ui_int = ints_%gemind(u,i)
          uu_int = ints_%gemind(u,u)

          ! u > v --> factor of 2

          do v = first_index_(u_sym,2) , u - 1

            uv_den  = dens_%gemind(u,v)
            tv_den  = dens_%gemind(t,v)
            uv_int  = ints_%gemind(u,v)
            vi_int  = ints_%gemind(i,v)

            ! 2 * d(tt|uv) * g(ii|uv)
            den_ind = pq_index(tt_den,uv_den)
            int_ind = pq_index(ii_int,uv_int)
            val     = val + 2.0_wp * den2(den_ind) * int2(int_ind)

            ! 4 * d(tv|tu) * g(ui|vi)
            den_ind = pq_index(tv_den,tu_den) + den_offset
            int_ind = pq_index(ui_int,vi_int) + int_offset
            val     = val + 4.0_wp * den2(den_ind) * int2(int_ind)

          end do

          ! u = v --> factor of 1

          ! d(tt|uu) * g(ii|uu)
          den_ind = pq_index(tt_den,uu_den)
          int_ind = pq_index(ii_int,uu_int)
          val     = val + den2(den_ind) * int2(int_ind)

          ! 2 * d(tu|tu) * g(ui|ui)
          den_ind = pq_index(tu_den,tu_den) + den_offset
          int_ind = pq_index(ui_int,ui_int) + int_offset
          val     = val + 2.0_wp * den2(den_ind) * int2(int_ind)
 
        end do

      end do

      ! 4 * sum{u}   [ ( delta(t,u) - d(t,u) ) * ( 3 * g(ui|ti) - g(ii|tu) ) ]

      do u = first_index_(t_sym,2) , last_index_(t_sym,2)

        tu_den  = dens_%gemind(t,u)
        ui_int  = ints_%gemind(u,i)
        tu_int  = ints_%gemind(t,u)

        if ( t == u ) then
          dfac  = 1.0_wp - den1(tu_den)
        else
          dfac  = - den1(tu_den)
        endif

        ! 3 * g(ui|ti) 
        int_ind = pq_index(ui_int,ti_int)
        int_val = 3.0_wp * int2(int_ind)        

        ! - g(ii|tu) 
        int_ind = pq_index(ii_int,tu_int)
        int_val = int_val - int2(int_ind)

        ! only factor of 2 because of multiplication below
 
        val = val + 2.0_wp * dfac * int_val 

      end do

      ! factor of 2 because only LT elements are accessed

      te_terms_ad = 2.0_wp * val

      return
  
    end function te_terms_ad

    function te_terms_ad_df(i,t,t_sym,int2,den1,den2)

      ! function to compute exact 2-e contribution for active-doubly occupied diagonal Hessian
      ! using 3-index integrals
      ! H(it|it) = 2 * sum{u,v} [ d(tt|uv) * g(ii|uv) + 2 * d(tv|tu) * g(ui|vi)         ]  
      !          + 4 * sum{u}   [ ( delta(t,u) - d(t,u) ) * ( 3 * g(ui|ti) - g(ii|tu) ) ]

      integer, intent(in)  :: i,t,t_sym
      real(wp), intent(in) :: int2(:),den1(:),den2(:)

      integer              :: tt_den,uv_den,tv_den,tu_den,uu_den
      integer(ip)          :: ii,uv,ui,vi,ti,tu,uu,nQ
      integer              :: idf,adf,udf,vdf,tdf
      integer              :: u,v,u_sym,tu_sym,den_ind,den_offset
      real(wp)             :: te_terms_ad_df,dfac,val,int_val
      real(wp)             :: v_ui(df_vars_%nQ),v_ii(df_vars_%nQ),v_ti(df_vars_%nQ)

      val = 0.0_wp

      nQ     = int(df_vars_%nQ,kind=ip) 
      idf    = df_vars_%class_to_df_map(i)
      tdf    = df_vars_%class_to_df_map(t) 

      tt_den = dens_%gemind(t,t)
      ii     = df_pq_index(idf,idf) 
      ti     = df_pq_index(idf,tdf) 

      call my_dcopy(df_vars_%nQ,int2(ii+1:),df_vars_%Qstride,v_ii,1)
      call my_dcopy(df_vars_%nQ,int2(ti+1:),df_vars_%Qstride,v_ti,1)

      do u_sym = 1 , nirrep_

        tu_sym     = group_mult_tab_(t_sym,u_sym)
        den_offset = dens_%offset(tu_sym)

        do u = first_index_(u_sym,2) , last_index_(u_sym,2)

          tu_den = dens_%gemind(t,u)
          uu_den = dens_%gemind(u,u)

          udf    = df_vars_%class_to_df_map(u)

          ui     = df_pq_index(udf,idf) 
          uu     = df_pq_index(udf,udf) 

          call my_dcopy(df_vars_%nQ,int2(ui+1:),df_vars_%Qstride,v_ui,1)

          ! u > v --> factor of 2

          do v = first_index_(u_sym,2) , u - 1

            uv_den  = dens_%gemind(u,v)
            tv_den  = dens_%gemind(t,v)

            vdf     = df_vars_%class_to_df_map(v)

            uv      = df_pq_index(udf,vdf) 
            vi      = df_pq_index(idf,vdf) 

            ! 2 * d(tt|uv) * g(ii|uv)
            den_ind = pq_index(tt_den,uv_den)
            val     = val + 2.0_wp * den2(den_ind) * my_ddot(df_vars_%nQ,v_ii,1,int2(uv+1:),df_vars_%Qstride)

            ! 4 * d(tv|tu) * g(ui|vi)
            den_ind = pq_index(tv_den,tu_den) + den_offset
            val     = val + 4.0_wp * den2(den_ind) * my_ddot(df_vars_%nQ,v_ui,1,int2(vi+1:),df_vars_%Qstride)

          end do

          ! u = v --> factor of 1

          ! d(tt|uu) * g(ii|uu)
          den_ind = pq_index(tt_den,uu_den)
          val     = val + den2(den_ind) * my_ddot(df_vars_%nQ,v_ii,1,int2(uu+1:),df_vars_%Qstride)

          ! 2 * d(tu|tu) * g(ui|ui)
          den_ind = pq_index(tu_den,tu_den) + den_offset
          val     = val + 2.0_wp * den2(den_ind) * my_ddot(df_vars_%nQ,v_ui,1,v_ui,1)

        end do

      end do

      ! 4 * sum{u}   [ ( delta(t,u) - d(t,u) ) * ( 3 * g(ui|ti) - g(ii|tu) ) ]

      do u = first_index_(t_sym,2) , last_index_(t_sym,2)

        tu_den  = dens_%gemind(t,u)
 
        udf     = df_vars_%class_to_df_map(u)

        ui      = df_pq_index(udf,idf) 
        tu      = df_pq_index(tdf,udf) 

        if ( t == u ) then
          dfac  = 1.0_wp - den1(tu_den)
        else
          dfac  = - den1(tu_den)
        endif

        ! 3 * g(ui|ti) 
        int_val = 3.0_wp * my_ddot(df_vars_%nQ,v_ti,1,int2(ui+1:),df_vars_%Qstride)

        ! - g(ii|tu) 
        int_val = int_val - my_ddot(df_vars_%nQ,v_ii,1,int2(tu+1:),df_vars_%Qstride)

       ! only factor of 2 because of multiplication below

        val = val + 2.0_wp * dfac * int_val

      end do

      ! factor of 2 because only LT elements are accessed

      te_terms_ad_df = 2.0_wp * val

      return

    end function te_terms_ad_df

    integer function full_hessian_aa(f_i,q,z,int2,den1,den2)

      implicit none

      real(wp), intent(in) :: den1(:),int2(:),den2(:)
      real(wp), intent(in) :: z(:,:),q(:,:)
      type(matrix_block), intent(in) :: f_i(:)

      full_hessian_aa = allocate_full_hessian_data()

      full_hessian_aa = full_hessian_aa_1e(f_i,q,z,den1)

      full_hessian_aa = full_hessian_aa_2e_fast(int2,den2)

!      full_hessian_aa = check_full_hessian_aa()
!
!      stop

      return

    end function full_hessian_aa

    integer function allocate_full_hessian_data()

      implicit none

      integer :: p,q,h_p,ind

      if (allocated(full_orbital_hessian_%index_map) ) deallocate(full_orbital_hessian_%index_map)
      allocate(full_orbital_hessian_%index_map(nmo_tot_,nmo_tot_))
      full_orbital_hessian_%index_map = -1

      if (allocated(full_orbital_hessian_%npair_sym) ) deallocate(full_orbital_hessian_%npair_sym)
      allocate(full_orbital_hessian_%npair_sym(nirrep_))
      full_orbital_hessian_%npair_sym = -1

      ind = 0

      do h_p = 1 , nirrep_

        full_orbital_hessian_%npair_sym(h_p) = nactpi_(h_p) * ( nactpi_(h_p) - 1 ) / 2

        do p = first_index_(h_p,2) , last_index_(h_p,2)

          do q = p+1 , last_index_(h_p,2)

            ind = ind + 1

            full_orbital_hessian_%index_map(p,q) = ind
            full_orbital_hessian_%index_map(q,p) = ind

          end do

        end do

      end do

      h_p = sum(full_orbital_hessian_%npair_sym)
      if (allocated(full_orbital_hessian_%aa))        deallocate(full_orbital_hessian_%aa)
      allocate(full_orbital_hessian_%aa(h_p,h_p))
      full_orbital_hessian_%aa = 0.0_wp

      if (allocated(full_orbital_hessian_%HinvXg_aa)) deallocate(full_orbital_hessian_%HinvXg_aa)
      allocate(full_orbital_hessian_%HinvXg_aa(h_p))
      full_orbital_hessian_%HinvXg_aa = 0.0_wp

      return

    end function allocate_full_hessian_data

    integer function deallocate_full_hessian_data()

      implicit none

      if (allocated(full_orbital_hessian_%index_map) ) deallocate(full_orbital_hessian_%index_map)

      if (allocated(full_orbital_hessian_%npair_sym) ) deallocate(full_orbital_hessian_%npair_sym)

      if (allocated(full_orbital_hessian_%aa) )        deallocate(full_orbital_hessian_%aa)

      if (allocated(full_orbital_hessian_%HinvXg_aa) ) deallocate(full_orbital_hessian_%HinvXg_aa)

      return

    end function deallocate_full_hessian_data

    integer function full_hessian_preconditioner(full_hessian,g)

      implicit none

      real(wp), intent(inout) :: full_hessian(:,:),g(:)
      integer, allocatable   :: ipiv(:)
      integer :: mat_dim,info

      ! determine dimension of H matrix
      mat_dim = size(full_hessian,dim=1)

      ! allocate scratch array for pivot indeces
      allocate(ipiv(mat_dim))

      ! factorize H matrix
      call dgetrf(mat_dim,mat_dim,full_hessian,mat_dim,ipiv,info)

      ! collect appropriate elements of orbital gradient and store
      ! in format compatible with storage of Hessian
      full_hessian_preconditioner = collect_gH(orbital_gradient_,g,'g2H')

      ! solve matrix equation H*x=g for x
      call dgetrs('n',mat_dim,1,full_hessian,mat_dim,ipiv,g,mat_dim,info)

      ! adjust sign of step vector
      call my_dscal(mat_dim,-1.0_wp,g,1)

      ! deallocate scratch array
      deallocate(ipiv)

      full_hessian_preconditioner = info

    end function full_hessian_preconditioner

    integer function collect_gH(vec_grad,vec_hess,direction)

      implicit none
      real(wp), intent(inout) :: vec_hess(:),vec_grad(:)
      character(3) :: direction

      integer :: t_sym,grad_ind,t,u,hess_ind,dir

      dir = 1
      if ( direction == 'H2g' ) dir = -1
      if ( direction == 'h2g' ) dir = -1
      if ( direction == 'h2G' ) dir = -1
      if ( direction == 'H2G' ) dir = -1

      do t_sym = 1 , nirrep_

        grad_ind = rot_pair_%pair_offset(t_sym,rot_pair_%act_act_type)

        do u = first_index_(t_sym,2) , last_index_(t_sym,2)

          do t = u + 1 , last_index_(t_sym,2)

            grad_ind = grad_ind + 1
            hess_ind = full_orbital_hessian_%index_map(t,u)

            if ( dir == 1 ) then

              vec_hess(hess_ind) = vec_grad(grad_ind)

            else

              vec_grad(grad_ind) = vec_hess(hess_ind)

            end if

          end do

        end do

      end do

      collect_gH = 0

    end function collect_gH

    integer function full_hessian_aa_1e(f_i,q,z,den1)

      real(wp), intent(in) :: den1(:)
      real(wp), intent(in) :: z(:,:),q(:,:)
      type(matrix_block), intent(in) :: f_i(:)

      integer :: h_p
      integer :: p_e,q_e,r_e,s_e
      integer :: p_i,q_i,r_i,s_i
      integer :: pq_hess,rs_hess
      integer :: qr_den,pr_den,qs_den,ps_den

      real(wp) :: val

      full_hessian_aa_1e = 0

      do h_p = 1 , nirrep_

        do p_e = first_index_(h_p,2) , last_index_(h_p,2)

          p_i    = trans_%class_to_irrep_map(p_e)

          do q_e = p_e + 1 , last_index_(h_p,2)

            q_i     = trans_%class_to_irrep_map(q_e)
            pq_hess = full_orbital_hessian_%index_map(p_e,q_e)

            do r_e = first_index_(h_p,2) , last_index_(h_p,2)

              r_i    = trans_%class_to_irrep_map(r_e)

              do s_e = r_e + 1 , last_index_(h_p,2)

                s_i     = trans_%class_to_irrep_map(s_e)
                rs_hess = full_orbital_hessian_%index_map(r_e,s_e)

                pr_den  = dens_%gemind(p_e,r_e)
                ps_den  = dens_%gemind(p_e,s_e)
                qr_den  = dens_%gemind(q_e,r_e)
                qs_den  = dens_%gemind(q_e,s_e)

                val     = 2.0_wp * ( - f_i(h_p)%val(p_i,s_i)*den1(qr_den) &
                                   & + f_i(h_p)%val(q_i,s_i)*den1(pr_den) &
                                   & + f_i(h_p)%val(p_i,r_i)*den1(qs_den) &
                                   & - f_i(h_p)%val(q_i,r_i)*den1(ps_den) )

                if ( q_e == r_e ) val = val + (  z(p_e-ndoc_tot_,s_e) &
                                             & + q(p_e-ndoc_tot_,s_e) &
                                             & + z(s_e-ndoc_tot_,p_e) &
                                             & + q(s_e-ndoc_tot_,p_e) )

                if ( p_e == r_e ) val = val - (  z(q_e-ndoc_tot_,s_e) &
                                             & + q(q_e-ndoc_tot_,s_e) &
                                             & + z(s_e-ndoc_tot_,q_e) &
                                             & + q(s_e-ndoc_tot_,q_e) )

                if ( q_e == s_e ) val = val - (  z(p_e-ndoc_tot_,r_e) &
                                             & + q(p_e-ndoc_tot_,r_e) &
                                             & + z(r_e-ndoc_tot_,p_e) &
                                             & + q(r_e-ndoc_tot_,p_e) )

                if ( p_e == s_e ) val = val + (  z(q_e-ndoc_tot_,r_e) &
                                             & + q(q_e-ndoc_tot_,r_e) &
                                             & + z(r_e-ndoc_tot_,q_e) &
                                             & + q(r_e-ndoc_tot_,q_e) )

                if (pq_hess /= rs_hess ) cycle

                full_orbital_hessian_%aa(pq_hess,rs_hess) = val

              end do

            end do

          end do

        end do

      end do

      return

    end function full_hessian_aa_1e

    integer function full_hessian_aa_2e_fast(int2,den2)

      implicit none
 
      real(wp), intent(in) :: int2(:),den2(:)

      integer :: p,q,r,s
      integer :: pq_hess,rs_hess
      integer :: h_p,h_r
      integer :: nmo_h_p,nmo_h_r

      integer :: r_min,s_min

      real(wp) :: val

      type(matrix_block) :: int2_aa(nirrep_)

      full_hessian_aa_2e_fast = allocate_int2_aa()

      full_hessian_aa_2e_fast = precompute_int2_aa()

      do h_p = 1 , nirrep_

        nmo_h_p = nactpi_(h_p)

        if ( ( nmo_h_p < 2 ) ) cycle

        do p = first_index_(h_p,2) , last_index_(h_p,2)

          do q = p + 1 , last_index_(h_p,2)

            pq_hess = full_orbital_hessian_%index_map(p,q)

            do h_r = h_p , nirrep_

              nmo_h_r = nactpi_(h_r)

              if ( ( nmo_h_r < 2 ) ) cycle

              r_min = first_index_(h_r,2)
              if ( h_r == h_p ) r_min = p

              do r = r_min , last_index_(h_r,2)

                s_min = r + 1
                if ( r == p ) s_min = q

                do s = s_min , last_index_(h_r,2)

                  rs_hess = full_orbital_hessian_%index_map(r,s)

                  full_orbital_hessian_%aa(rs_hess,pq_hess) =       &
                      & full_orbital_hessian_%aa(rs_hess,pq_hess) + &
                      & 2.0_wp * (                                  &
                      &     y_fast(p,r,q,s,h_p,h_r)                 &
                      &   - y_fast(q,r,p,s,h_p,h_r)                 &
                      &   - y_fast(p,s,q,r,h_p,h_r)                 &
                      &   + y_fast(q,s,p,r,h_p,h_r) )

                  write(*,'(2(i3,1x),5x,4(i2,1x),5x,es20.8)')rs_hess,pq_hess,r,s,p,q, &
                    & full_orbital_hessian_%aa(rs_hess,pq_hess)

                end do

              end do

            end do

          end do

        end do

      end do

      !stop

      full_hessian_aa_2e_fast = deallocate_int2_aa()

      full_hessian_aa_2e_fast = 0

      return

      contains

        integer function precompute_int2_aa()

          implicit none

          integer :: max_ngem

          real(wp), allocatable :: scr(:,:)
          integer               :: h_p,h_q,h_pq
          integer               :: nmo_h_p,nmo_h_q,ngem_h_pq
          integer               :: p,q,q_max,p_df,q_df,pq
          integer               :: nQ
          integer(ip)           :: pq_df 

          nQ    = int(df_vars_%nQ,kind=ip)

          max_ngem = maxval(dens_%ngempi)

          allocate(scr(df_vars_%nQ,max_ngem))

          do h_pq = 1 , nirrep_

            ngem_h_pq = dens_%ngempi(h_pq)

            if ( ngem_h_pq == 0 ) cycle

            scr(:,1:ngem_h_pq) = 0.0_wp 
           
            ! gather 3-index integrals (:,pq) for this h_pq

            do h_p = 1 , nirrep_

              h_q     = group_mult_tab_(h_pq,h_p)

              if ( h_q > h_p ) cycle 

              nmo_h_p = nactpi_(h_p)
              nmo_h_q = nactpi_(h_q)
           
              do p = first_index_(h_p,2) , last_index_(h_p,2)

                p_df = df_vars_%class_to_df_map(p) 

                q_max = last_index_(h_q,2)
                if ( h_p == h_q ) q_max = p
 
                do q = first_index_(h_q,2) , q_max

                  q_df  = df_vars_%class_to_df_map(q)
                  pq_df = df_pq_index(p_df,q_df) 

                  pq   = dens_%gemind(p,q)

                  call my_dcopy(df_vars_%nQ,int2(pq_df+1:),df_vars_%Qstride,&
                        & scr(:,pq),1)
              
                end do
 
              end do
 
            end do
   
            ! calculate all 4-index integrals (rs|pq)

            call dgemm('t','n',ngem_h_pq,ngem_h_pq,nQ,1.0_wp,    &
                  & scr(:,1:ngem_h_pq),nQ,scr(:,1:ngem_h_pq),nQ, &
                  & 0.0_wp,int2_aa(h_pq)%val,ngem_h_pq)

          end do

          deallocate(scr)

          precompute_int2_aa = 0

          return

        end function precompute_int2_aa

        integer function allocate_int2_aa()

          ! function to allocate temporary matrix for active-active 2-e integrals

          implicit none

          integer :: h_pq,ngem_pq

          do h_pq = 1 , nirrep_

            ngem_pq = dens_%ngempi(h_pq)
 
            if (ngem_pq == 0) cycle
 
            allocate(int2_aa(h_pq)%val(ngem_pq,ngem_pq))

            int2_aa(h_pq)%val = 0.0_wp

          end do

          allocate_int2_aa = 0

          return

        end function allocate_int2_aa 

        integer function deallocate_int2_aa()

          ! function to allocate temporary matrix for active-active 2-e integrals 

          implicit none

          integer :: h_pq,ngem_pq

          do h_pq = 1 , nirrep_

            ngem_pq = dens_%ngempi(h_pq)

            if (ngem_pq == 0) cycle

            deallocate(int2_aa(h_pq)%val)

          end do

          deallocate_int2_aa = 0

          return

        end function deallocate_int2_aa

        function y_fast(p,q,r,s,h_1,h_2)

          ! it is implied that h_p == h_r == h_1 from above
          ! it is implies that h_q == h_s == h_2 from above

          implicit none

          integer, intent(in) :: p,q,r,s,h_1,h_2
          integer :: h_p,h_q,h_r,h_s,h_w,h_x
          integer :: h_rw,h_sx,h_rs,h_wx

          integer :: pw_int,qx_int,pq_int,wx_int
          integer :: int_row,int_col

          integer :: w,x

          integer :: sx_den,rw_den,rs_den,wx_den,pq_den,ind_den,rw_off,rs_off

          real(wp) :: y_fast

          h_p   = h_1
          h_q   = h_2
          h_r   = h_1
          h_s   = h_2

          h_rs  = group_mult_tab_(h_r,h_s)

          pq_int = dens_%gemind(p,q)
          rs_den = dens_%gemind(r,s)
          rs_off = dens_%offset(h_rs)

          y_fast = 0.0_wp

          do h_w = 1 , nirrep_

            h_rw   = group_mult_tab_(h_r,h_w)  ! same value as h_pw
            h_x    = group_mult_tab_(h_s,h_rw)
            rw_off = dens_%offset(h_rw)

            do w = first_index_(h_w,2) , last_index_(h_w,2)

              pw_int = dens_%gemind(w,p)

              do x = first_index_(h_x,2) , last_index_(h_x,2)

                rw_den  = dens_%gemind(r,w)
                sx_den  = dens_%gemind(s,x)
                ind_den = rw_off + pq_index(rw_den,sx_den)

                qx_int = dens_%gemind(q,x)
                 
                if ( qx_int > pw_int) then
                
                  int_row = qx_int
                  int_col = pw_int

                else

                  int_row = pw_int
                  int_col = qx_int 

                end if

                ! 2*g(pw|qx)*d2(rw|sx)
                y_fast = y_fast + 2.0_wp * den2(ind_den) * int2_aa(h_rw)%val(int_row,int_col)

              end do

            end do

            ! h_x = group_mult_tab_(h_rs,h_w)

            do w = first_index_(h_w,2) , last_index_(h_w,2)

              do x = first_index_(h_x,2) , last_index_(h_x,2)

                wx_den  = dens_%gemind(w,x)
                ind_den = rs_off + pq_index(rs_den,wx_den)

                wx_int = dens_%gemind(w,x)

                if ( wx_int > pq_int) then

                  int_row = wx_int
                  int_col = pq_int

                else

                  int_row = pq_int
                  int_col = wx_int

                end if                

                ! g(pq|wx)*d2(rs|wx)
                y_fast = y_fast + den2(ind_den) * int2_aa(h_rs)%val(int_row,int_col)

              end do

            end do

          end do

          return

        end function y_fast

    end function full_hessian_aa_2e_fast

    integer function full_hessian_aa_2e_slow(int2,den2)

      implicit none

      real(wp), intent(in) :: int2(:),den2(:)

      integer :: p,q,r,s
      integer :: pq_hess,rs_hess
      integer :: h_p,h_r
      integer :: nmo_h_p,nmo_h_r

      integer :: r_min,s_min

      real(wp) :: val

      do h_p = 1 , nirrep_

        nmo_h_p = nactpi_(h_p)

        if ( ( nmo_h_p < 2 ) ) cycle

        do p = first_index_(h_p,2) , last_index_(h_p,2)

          do q = p + 1 , last_index_(h_p,2)

            pq_hess = full_orbital_hessian_%index_map(p,q)

            do h_r = h_p , nirrep_

              nmo_h_r = nactpi_(h_r)

              if ( ( nmo_h_r < 2 ) ) cycle

              r_min = first_index_(h_r,2)
              if ( h_r == h_p ) r_min = p

              do r = r_min , last_index_(h_r,2)

                s_min = r + 1
                if ( r == p ) s_min = q

                do s = s_min , last_index_(h_r,2)

                  rs_hess = full_orbital_hessian_%index_map(r,s)

!                  if (rs_hess /= pq_hess) cycle

                  full_orbital_hessian_%aa(rs_hess,pq_hess) =       &
                      & full_orbital_hessian_%aa(rs_hess,pq_hess) + &
                      & 2.0_wp * (                                  &
                      &     y_slow(p,r,q,s,h_p,h_r)                 &
                      &   - y_slow(q,r,p,s,h_p,h_r)                 &
                      &   - y_slow(p,s,q,r,h_p,h_r)                 &
                      &   + y_slow(q,s,p,r,h_p,h_r) )

!                  write(*,'(2(i2,1x),4x,2(i2,1x),4x,es20.12)')r,s,p,q,full_orbital_hessian_%aa(rs_hess,pq_hess)

                end do

              end do

            end do

          end do

        end do

      end do

!        ! off-diagonal symmetry block
!
!        do h_r = h_p + 1 , nirrep_
!
!          nmo_h_r = nactpi_(h_r)
!
!          if ( ( nmo_h_p < 2 ) .or. ( nmo_h_r < 2 ) ) cycle          
!
!          do p = first_index_(h_p,2) , last_index_(h_p,2)
!
!            do q = p + 1 , last_index_(h_p,2)
!
!              pq_hess = full_orbital_hessian_%index_map(p,q)
!
!
!              do r = first_index_(h_r,2) , last_index_(h_r,2)
!
!                do s = r + 1 , last_index_(h_r,2)
!
!                  rs_hess = full_orbital_hessian_%index_map(r,s)
!
!                  if ( rs_hess < pq_hess ) cycle
!
!                  full_orbital_hessian_%aa(pq_hess,rs_hess) =       &
!                      & full_orbital_hessian_%aa(pq_hess,rs_hess) + &
!                      & 2.0_wp * (                                  &
!                      &     y_slow(p,r,q,s,h_p,h_r)                 &
!                      &   - y_slow(q,r,p,s,h_p,h_r)                 &
!                      &   - y_slow(p,s,q,r,h_p,h_r)                 &
!                      &   + y_slow(q,s,p,r,h_p,h_r) )
!
!                end do
!
!              end do
!
!            end do
!
!          end do
!
!        end do
!
!      end do

      full_hessian_aa_2e_slow = 1

      return

      contains

        function y_slow(p,q,r,s,h_1,h_2)

          ! it is implied that h_p == h_r == h_1 from above
          ! it is implies that h_q == h_s == h_2 from above

          implicit none

          integer, intent(in) :: p,q,r,s,h_1,h_2
          integer :: h_p,h_q,h_r,h_s,h_w,h_x
          integer :: h_rw,h_sx,h_rs,h_wx

          integer :: pdf,qdf,rdf,sdf,wdf,xdf
          integer(ip) :: pw_df,qx_df,pq_df,wx_df,nQ

          integer :: w,x

          integer :: sx_den,rw_den,rs_den,wx_den,pq_den,ind_den,rw_off,rs_off

          real(wp) :: y_slow

          real(wp) :: v_pw(df_vars_%nQ),v_pq(df_vars_%nQ)

          nQ    = int(df_vars_%nQ,kind=ip)

          h_p   = h_1
          h_q   = h_2
          h_r   = h_1
          h_s   = h_2

          h_rs  = group_mult_tab_(h_r,h_s)

          pdf   = df_vars_%class_to_df_map(p)
          qdf   = df_vars_%class_to_df_map(q)
          rdf   = df_vars_%class_to_df_map(r)
          sdf   = df_vars_%class_to_df_map(s)

          pq_df = df_pq_index(pdf,qdf)
          call my_dcopy(df_vars_%nQ,int2(pq_df+1:),df_vars_%Qstride,&
                        & v_pq,1)

          rs_den = dens_%gemind(r,s)
          rs_off = dens_%offset(group_mult_tab_(h_r,h_s))

          y_slow = 0.0_wp

          do h_w = 1 , nirrep_

            h_rw   = group_mult_tab_(h_r,h_w)  ! same value as h_pw
            h_x    = group_mult_tab_(h_s,h_rw)
            rw_off = dens_%offset(h_rw)

            do w = first_index_(h_w,2) , last_index_(h_w,2)


              wdf    = df_vars_%class_to_df_map(w)
              pw_df  = df_pq_index(pdf,wdf)

              call my_dcopy(df_vars_%nQ,int2(pw_df+1:),df_vars_%Qstride,&
                        & v_pw,1)

              do x = first_index_(h_x,2) , last_index_(h_x,2)

                rw_den  = dens_%gemind(r,w)
                sx_den  = dens_%gemind(s,x)
                ind_den = rw_off + pq_index(rw_den,sx_den)

                xdf     = df_vars_%class_to_df_map(x)
                qx_df   = df_pq_index(qdf,xdf)

                ! 2*g(pw|qx)*d2(rw|sx)
                y_slow = y_slow + 2.0_wp * den2(ind_den) * &
                       & my_ddot(df_vars_%nQ,int2(qx_df+1:),df_vars_%Qstride,v_pw,1)

              end do

            end do

            ! h_x = group_mult_tab_(h_rs,h_w)

            do w = first_index_(h_w,2) , last_index_(h_w,2)

              wdf = df_vars_%class_to_df_map(w)

              do x = first_index_(h_x,2) , last_index_(h_x,2)

                wx_den  = dens_%gemind(w,x)
                ind_den = rs_off + pq_index(rs_den,wx_den)

                xdf = df_vars_%class_to_df_map(x)
                wx_df  = df_pq_index(wdf,xdf)

                ! g(pq|wx)*d2(rs|wx)
                y_slow = y_slow + den2(ind_den) * &
                       & my_ddot(df_vars_%nQ,int2(wx_df+1:),df_vars_%Qstride,v_pq,1)

              end do

            end do

          end do

          return

        end function y_slow

    end function full_hessian_aa_2e_slow

    integer function check_full_hessian_aa()

      implicit none

      integer :: t_sym,grad_ind,hess_ind,t,u

      real(wp) :: diag,full

      do t_sym = 1 , nirrep_

        grad_ind = rot_pair_%pair_offset(t_sym,rot_pair_%act_act_type)

        do u = first_index_(t_sym,2) , last_index_(t_sym,2)

          do t = u + 1 , last_index_(t_sym,2)

            grad_ind = grad_ind + 1
            hess_ind = full_orbital_hessian_%index_map(t,u)

            full = full_orbital_hessian_%aa(hess_ind,hess_ind)
            diag = orbital_hessian_(grad_ind)

            if ( abs(full - diag) < 1.0e-12_wp ) cycle

            write(*,'(2(i4,1x),4x,es10.3,5x,2(es14.6,1x))')t,u,diag-full,diag,full

          end do

        end do

      end do

      write(*,*)'done checking aa block of diagonal hessian'

      check_full_hessian_aa = 0

      return

    end function check_full_hessian_aa

    integer function diagonal_hessian_aa(f_i,f_a,q,z,int2,den1,den2)
      implicit none
      ! subroutine to compute the exact diagonal Hessian matrix element H(tu|tu)
      ! according to Eq. A.3 in Jensen, Chem. Phys. 104, 229, (1982)
      ! H(tu|tu) = 2 * ( d(t,t) * f_i(u,u) + d(u,u) * f_i(t,t) - 2 * d(u,t) * f_i(u,t) 
      !                  - ( q(u,u) + q(t,t) + z(t,t) +z(u,u) ) )
      ! some 2-e contributions are neglected
      real(wp), intent(in) :: den1(:),int2(:),den2(:)
      real(wp), intent(in) :: z(:,:),q(:,:)
      type(matrix_block), intent(in) :: f_i(:),f_a(:)
      integer :: t_sym,t,u,grad_ind
      integer :: t_i,u_i
      integer :: tt_den,uu_den,ut_den
      integer :: tt_int,uu_int,ut_int
      real(wp) :: h_val

      ! loop over all active - active pairs

      do t_sym = 1 , nirrep_

        grad_ind = rot_pair_%pair_offset(t_sym,rot_pair_%act_act_type)

        do u = first_index_(t_sym,2) , last_index_(t_sym,2)

          do t = u + 1 , last_index_(t_sym,2)

            grad_ind = grad_ind + 1

            ! compute Hessian value

            ! Fock matrix addressing
            u_i    = trans_%class_to_irrep_map(u)
            t_i    = trans_%class_to_irrep_map(t)

            ! Density addressing
            tt_den = dens_%gemind(t,t)
            uu_den = dens_%gemind(u,u)
            ut_den = dens_%gemind(u,t)

            if ( use_exact_hessian_diagonal_ == 1 ) then

              h_val    = 2.0_wp * ( den1(tt_den) * f_i(t_sym)%val(u_i,u_i) +            &
                                 &  den1(uu_den) * f_i(t_sym)%val(t_i,t_i)              &
                                 & - (  2.0_wp * den1(ut_den) * f_i(t_sym)%val(t_i,u_i) &
                                 &    + q(t-ndoc_tot_,t) + q(u-ndoc_tot_,u)             &
                                 &    + z(t-ndoc_tot_,t) + z(u-ndoc_tot_,u) ) )

              if ( df_vars_%use_df_teints == 1 ) then

                h_val = h_val + te_terms_aa_df(t,u,t_sym,int2,den2)

              else

                h_val = h_val + te_terms_aa(t,u,t_sym,int2,den2)

              endif

            else

              h_val = 2.0_wp * (                                                             &
                    &   den1(tt_den) * ( f_i(t_sym)%val(u_i,u_i) + f_i(t_sym)%val(u_i,u_i) ) &
                    & + den1(uu_den) * ( f_i(t_sym)%val(t_i,t_i) + f_i(t_sym)%val(t_i,t_i) ) &                      
                    & - ( q(t-ndoc_tot_,t) + z(t-ndoc_tot_,t) )                              &
                    & - ( q(u-ndoc_tot_,u) + z(u-ndoc_tot_,u) )                              &
                    & - 2.0_wp * den1(ut_den) * f_i(t_sym)%val(u_i,t_i)                      &
                    & )

            end if

            if ( h_val < 0.0_wp ) h_val = - h_val

            if ( h_val < min_diag_hessian_ ) h_val = min_diag_hessian_

            orbital_hessian_(grad_ind) = h_val

          end do

        end do

      end do

      diagonal_hessian_aa = 0

      return

    end function diagonal_hessian_aa

    function te_terms_aa_df(t,u,t_sym,int2,den2)

      ! function to accumulate 2-e contributions to exact Hessian active-active diagonal element
      ! using 4-index integrals according to Eq. A3 in Jensen Chem. Phys. 104, 229 (1982)
      ! val = sum(x,y) [   4 * d(tx|ty) * g(ux|uy) + 4 * d(ux|uy) * g(tx|ty)
      !                  + 2 * d(tt|xy) * g(uu|xy) + 2 * d(uu|xy) * g(tt|xy)
      !                  - 8 * d(ux|ty) * g(ux|ty) - 4 * d(ut|xy) * g(ut|xy) ]  
      ! the operations are performed in two steps
      ! i)  x>y --> factor of 2 (use coefficients above)
      ! ii) x=y --> factor of 1 (use 1/2 the coefficients above)
      ! and then multiplying the result by 2 at the end

      real(wp)             :: te_terms_aa
      real(wp), intent(in) :: int2(:),den2(:)
      integer, intent(in)  :: t,u,t_sym
      integer              :: x,y
      integer              :: ux_den,uy_den,tx_den,ty_den,tt_den,uu_den,ut_den,xy_den,xx_den
      integer              :: udf,tdf,xdf,ydf
      integer(ip)          :: ux,uy,tx,ty,uu,tt,ut,xy,xx,nQ
      integer              :: den_offset,den_ind
      integer              :: x_sym,xt_sym
      
      real(wp)             :: te_terms_aa_df,val
      real(wp)             :: v_uu(df_vars_%nQ),v_tt(df_vars_%nQ),v_ut(df_vars_%nQ)
      real(wp)             :: v_ux(df_vars_%nQ),v_tx(df_vars_%nQ),v_xy(df_varS_%nQ)     
      real(wp)             :: v_ty(df_vars_%nQ)

      nQ     = int(df_vars_%nQ,kind=ip)

      val    = 0.0_wp

      ! integral/density indeces

      tt_den = dens_%gemind(t,t)
      uu_den = dens_%gemind(u,u)
      ut_den = dens_%gemind(u,t)

      tdf    = df_vars_%class_to_df_map(t)
      udf    = df_vars_%class_to_df_map(u)

      tt     = df_pq_index(tdf,tdf) 
      uu     = df_pq_index(udf,udf) 
      ut     = df_pq_index(udf,tdf) 

      call my_dcopy(df_vars_%nQ,int2(tt+1:),df_vars_%Qstride,v_tt,1)
      call my_dcopy(df_vars_%nQ,int2(uu+1:),df_vars_%Qstride,v_uu,1)
      call my_dcopy(df_vars_%nQ,int2(ut+1:),df_vars_%Qstride,v_ut,1)

      ! loop over symmetries for x

      do x_sym = 1 , nirrep_

        ! loop over active x indeces

        do x = first_index_(x_sym,2) , last_index_(x_sym,2)

          xdf    = df_vars_%clasS_to_df_map(x)   

          xt_sym = group_mult_tab_(x_sym,t_sym)

          ! integral/density indeces
          ux_den = dens_%gemind(u,x)
          tx_den = dens_%gemind(t,x)
          xx_den = dens_%gemind(x,x)

          ux     = df_pq_index(udf,xdf) 
          tx     = df_pq_index(tdf,xdf) 
          xx     = df_pq_index(xdf,xdf) 

          call my_dcopy(df_vars_%nQ,int2(ux+1:),df_vars_%Qstride,v_ux,1)
          call my_dcopy(df_vars_%nQ,int2(tx+1:),df_vars_%Qstride,v_tx,1)

          den_offset = dens_%offset(xt_sym)

          ! x > y --> factor of 2

          do y = first_index_(x_sym,2) , x - 1

            ydf     = df_vars_%class_to_df_map(y)

            ! integral/density indeces
            uy_den  = dens_%gemind(u,y)
            ty_den  = dens_%gemind(t,y)
            xy_den  = dens_%gemind(x,y)

            uy      = df_pq_index(udf,ydf) 
            ty      = df_pq_index(tdf,ydf) 
            xy      = df_pq_index(xdf,ydf) 

            call my_dcopy(df_vars_%nQ,int2(xy+1:),df_vars_%Qstride,v_xy,1)
            call my_dcopy(df_vars_%nQ,int2(ty+1:),df_vars_%Qstride,v_ty,1)

            ! 4 * d(tx|ty) * g(ux|uy)
            den_ind = pq_index(tx_den,ty_den) + den_offset
            val     = val + 4.0_wp * den2(den_ind) * my_ddot(df_vars_%nQ,v_ux,1,int2(uy+1:),df_vars_%Qstride)

            ! 4 * d(ux|uy) * g(tx|ty)
            den_ind = pq_index(ux_den,uy_den) + den_offset
            val     = val + 4.0_wp * den2(den_ind) * my_ddot(df_vars_%nQ,v_tx,1,v_ty,1)

            ! 2 * d(tt,xy) * g(uu|xy)
            den_ind = pq_index(tt_den,xy_den)
            val     = val + 2.0_wp * den2(den_ind) * my_ddot(df_vars_%nQ,v_uu,1,v_xy,1)

            ! 2 * d(uu,xy) * g(tt|xy)
            den_ind = pq_index(uu_den,xy_den)
            val     = val + 2.0_wp * den2(den_ind) * my_ddot(df_vars_%nQ,v_tt,1,v_xy,1)

            ! - 8 * d(ux|ty) * g(ux|ty)
            den_ind = pq_index(ux_den,ty_den) + den_offset
            val     = val - 8.0_wp * den2(den_ind) * my_ddot(df_vars_%nQ,v_ux,1,v_ty,1)

            ! - 4 * d(tu|xy) * g(ut|xy)
            den_ind = pq_index(ut_den,xy_den)
            val     = val - 4.0_wp * den2(den_ind) * my_ddot(df_vars_%nQ,v_ut,1,v_xy,1)

          end do

          ! x == y

          call my_dcopy(df_vars_%nQ,int2(xx+1:),df_vars_%Qstride,v_xy,1)

          ! 2 * d(tx|tx) * g(ux|ux)
          den_ind = pq_index(tx_den,tx_den) + den_offset
          val     = val + 2.0_wp * den2(den_ind) * my_ddot(df_vars_%nQ,v_ux,1,v_ux,1)

          ! 2 * d(ux|ux) * g(tx|tx)
          den_ind = pq_index(ux_den,ux_den) + den_offset
          val     = val + 2.0_wp * den2(den_ind) * my_ddot(df_vars_%nQ,v_tx,1,v_tx,1)

          ! d(tt,xx) * g(uu|xx)
          den_ind = pq_index(tt_den,xx_den)
          val     = val + den2(den_ind) * my_ddot(df_vars_%nQ,v_uu,1,v_xy,1) 

          ! d(uu,xx) * g(tt|xx)
          den_ind = pq_index(uu_den,xx_den)
          val     = val + den2(den_ind) * my_ddot(df_vars_%nQ,v_tt,1,v_xy,1)

          ! - 4 * d(ux|tx) * g(ux|tx)

          den_ind = pq_index(ux_den,tx_den) + den_offset
          val     = val - 4.0_wp * den2(den_ind) * my_ddot(df_vars_%nQ,v_ux,1,v_tx,1)

          ! - 2 * d(tu|xx) * g(ut|xx)
          den_ind = pq_index(ut_den,xx_den)
          val     = val - 2.0_wp * den2(den_ind) * my_ddot(df_vars_%nQ,v_ut,1,v_xy,1)

        end do

      end do

      ! factor of 2 because only LT elements are accessed 

      te_terms_aa_df = 2.0_wp * val

      return

    end function te_terms_aa_df

    function te_terms_aa(t,u,t_sym,int2,den2)

      ! function to accumulate 2-e contributions to exact Hessian active-active diagonal element
      ! using 4-index integrals according to Eq. A3 in Jensen Chem. Phys. 104, 229 (1982)
      ! val = sum(x,y) [   4 * d(tx|ty) * g(ux|uy) + 4 * d(ux|uy) * g(tx|ty)
      !                  + 2 * d(tt|xy) * g(uu|xy) + 2 * d(uu|xy) * g(tt|xy)
      !                  - 8 * d(ux|ty) * g(ux|ty) - 4 * d(ut|xy) * g(ut|xy) ]  
      ! the operations are performed in two steps
      ! i)  x>y --> factor of 2 (use coefficients above)
      ! ii) x=y --> factor of 1 (use 1/2 the coefficients above)
      ! and then multiplying the result by 2 at the end

      real(wp)             :: te_terms_aa
      real(wp), intent(in) :: int2(:),den2(:)
      integer, intent(in)  :: t,u,t_sym
      integer              :: x,y
      integer              :: ux_den,uy_den,tx_den,ty_den,tt_den,uu_den,ut_den,xy_den,xx_den
      integer              :: ux_int,uy_int,tx_int,ty_int,uu_int,tt_int,ut_int,xy_int,xx_int
      integer              :: den_offset,int_offset,den_ind,int_ind
      integer              :: x_sym,xt_sym
      real(wp)             :: val
 
      ! initialize 2-e contribution

      val = 0.0_wp

      ! integral/density indeces

      tt_den = dens_%gemind(t,t)
      uu_den = dens_%gemind(u,u)
      ut_den = dens_%gemind(u,t)
      tt_int = ints_%gemind(t,t)
      uu_int = ints_%gemind(u,u)
      ut_int = ints_%gemind(u,t)

      ! loop over symmetries for x

      do x_sym = 1 , nirrep_

        ! loop over active x indeces

        do x = first_index_(x_sym,2) , last_index_(x_sym,2)

          xt_sym     = group_mult_tab_(x_sym,t_sym)

          ! integral/density indeces
          ux_den = dens_%gemind(u,x)
          tx_den = dens_%gemind(t,x)
          xx_den = dens_%gemind(x,x)
          ux_int = ints_%gemind(u,x)
          tx_int = ints_%gemind(t,x)
          xx_int = ints_%gemind(x,x)

          den_offset = dens_%offset(xt_sym) 
          int_offset = ints_%offset(xt_sym)

          ! x > y --> factor of 2

          do y = first_index_(x_sym,2) , x - 1

            ! integral/density indeces
            uy_den = dens_%gemind(u,y)
            ty_den = dens_%gemind(t,y)
            xy_den = dens_%gemind(x,y)
            uy_int = ints_%gemind(u,y)
            ty_int = ints_%gemind(t,y)
            xy_int = ints_%gemind(x,y)

            ! 4 * d(tx|ty) * g(ux|uy)

            den_ind = pq_index(tx_den,ty_den) + den_offset
            int_ind = pq_index(ux_int,uy_int) + int_offset
            val = val + 4.0_wp * den2(den_ind) * int2(int_ind)

            ! 4 * d(ux|uy) * g(tx|ty)

            den_ind = pq_index(ux_den,uy_den) + den_offset
            int_ind = pq_index(tx_int,ty_int) + int_offset
            val = val + 4.0_wp * den2(den_ind) * int2(int_ind)

            ! 2 * d(tt,xy) * g(uu|xy)
            den_ind = pq_index(tt_den,xy_den) 
            int_ind = pq_index(uu_int,xy_int)
            val = val + 2.0_wp * den2(den_ind) * int2(int_ind)

            ! 2 * d(uu,xy) * g(tt|xy)
            den_ind = pq_index(uu_den,xy_den) 
            int_ind = pq_index(tt_int,xy_int)
            val = val + 2.0_wp * den2(den_ind) * int2(int_ind)

            ! - 8 * d(ux|ty) * g(ux|ty)

            den_ind = pq_index(ux_den,ty_den) + den_offset
            int_ind = pq_index(ux_int,ty_int) + int_offset
            val = val - 8.0_wp * den2(den_ind) * int2(int_ind)

            ! - 4 * d(tu|xy) * g(ut|xy)
            den_ind = pq_index(ut_den,xy_den)
            int_ind = pq_index(ut_int,xy_int) 
            val = val - 4.0_wp * den2(den_ind) * int2(int_ind)

          end do

          ! x == y
 
          ! 2 * d(tx|tx) * g(ux|ux)

          den_ind = pq_index(tx_den,tx_den) + den_offset
          int_ind = pq_index(ux_int,ux_int) + int_offset
          val = val + 2.0_wp * den2(den_ind) * int2(int_ind)

          ! 2 * d(ux|ux) * g(tx|tx)

          den_ind = pq_index(ux_den,ux_den) + den_offset
          int_ind = pq_index(tx_int,tx_int) + int_offset
          val = val + 2.0_wp * den2(den_ind) * int2(int_ind)

          ! d(tt,xx) * g(uu|xx)
          den_ind = pq_index(tt_den,xx_den)
          int_ind = pq_index(uu_int,xx_int)
          val = val + den2(den_ind) * int2(int_ind)

          ! d(uu,xx) * g(tt|xx)
          den_ind = pq_index(uu_den,xx_den)
          int_ind = pq_index(tt_int,xx_int)
          val = val + den2(den_ind) * int2(int_ind)

          ! - 4 * d(ux|tx) * g(ux|tx)

          den_ind = pq_index(ux_den,tx_den) + den_offset
          int_ind = pq_index(ux_int,tx_int) + int_offset
          val = val - 4.0_wp * den2(den_ind) * int2(int_ind)

          ! - 2 * d(tu|xx) * g(ut|xx)
          den_ind = pq_index(ut_den,xx_den)
          int_ind = pq_index(ut_int,xx_int)
          val = val - 2.0_wp * den2(den_ind) * int2(int_ind)

        end do

      end do

      ! factor of 2 because only LT elements are accessed      

      te_terms_aa = 2.0_wp * val

      return

    end function te_terms_aa

    integer function diagonal_hessian_ed(f_i,f_i_ext,f_a,f_a_ext,int2)
      implicit none
      ! subroutine to compute the diagonal Hessian matrix element H(ia|ia)according to Eq. 4.7a of
      ! Chaban, Schmidt, Gordon, Theor. Chem. Acc., 97, 88-95 (1997) 
      ! H(ai|ai) = 4 * (f_i(a,a) + f_a(a,a) ) - 4 * ( f_i(i,i) + f_a(i,i)) 
      ! a \in ext & i \in act
      real(wp), intent(in)           :: int2(:)
      type(vector_block), intent(in) :: f_i_ext(:),f_a_ext(:)
      type(matrix_block), intent(in) :: f_i(:),f_a(:)

      integer :: a,i,grad_ind,a_sym,i_i,a_i
      real(wp) :: h_val

      ! loop over all external - doubly-occupied pairs

      do a_sym = 1 , nirrep_

        grad_ind = rot_pair_%pair_offset(a_sym,rot_pair_%ext_doc_type)

        do i = first_index_(a_sym,1) , last_index_(a_sym,1)

          do a = first_index_(a_sym,3) , last_index_(a_sym,3)

            ! update Hessian index

            grad_ind = grad_ind + 1

            ! inactive/active fock matrix indeces

            i_i     = trans_%class_to_irrep_map(i)
            a_i     = trans_%class_to_irrep_map(a)-ndocpi_(a_sym)-nactpi_(a_sym)

            ! calculate diagonal Hessian element

            h_val    = 4.0 * ( f_i_ext(a_sym)%val(a_i) + f_a_ext(a_sym)%val(a_i) - &
                             & f_i(a_sym)%val(i_i,i_i) - f_a(a_sym)%val(i_i,i_i) )

            if ( use_exact_hessian_diagonal_ == 1 ) then

              if ( df_vars_%use_df_teints == 1 ) then

                h_val = h_val + te_terms_ed_df(i,a,int2)

              else

                h_val = h_val + te_terms_ed(i,a,int2)

              endif

            endif

            if ( h_val < 0.0_wp ) h_val = - h_val

            if ( h_val < min_diag_hessian_ ) h_val = min_diag_hessian_

            orbital_hessian_(grad_ind) = h_val

          end do

        end do
 
      end do

      diagonal_hessian_ed = 0

      return

    end function diagonal_hessian_ed

    function te_terms_ed(i,a,int2)

      ! function to compute the 2-e contributions to the external-doubly occupied diagonal Hessian
      ! using 4-index integrals element according to
      ! te_(i,a) = 4 (3 g(ai|ai) - (aa|ii) )

      integer, intent(in)  :: i,a
      real(wp), intent(in) :: int2(:)

      integer              :: aa_int,ii_int,ai_int,int_ind
      real(wp)             :: te_terms_ed,val

      aa_int      = ints_%gemind(a,a)
      ii_int      = ints_%gemind(i,i)
      ai_int      = ints_%gemind(a,i)

      ! 3 * g(ai|ai)
      int_ind     = pq_index(ai_int,ai_int) 
      val         = 2.0_wp * int2(int_ind)

      ! - g(aa|ii)
      int_ind     = pq_index(aa_int,ii_int)
      val         = val - int2(int_ind)
 
      ! factor of 4 from overall formula

      te_terms_ed = 4.0_wp * val

      return

    end function te_terms_ed

    function te_terms_ed_df(i,a,int2)

      ! function to compute the 2-e contributions to the external-doubly occupied diagonal Hessian
      ! using 3-index integrals element according to
      ! te_(i,a) = 4 (3 g(ai|ai) - (aa|ii) )

      integer, intent(in)  :: i,a
      real(wp), intent(in) :: int2(:)

      integer(ip)          :: ii,ai,aa,nQ
      integer              :: adf,idf
      real(wp)             :: te_terms_ed_df,val

      nQ      = int(df_vars_%nQ,kind=ip)

      adf     = df_vars_%class_to_df_map(a)
      idf     = df_vars_%class_to_df_map(i)

      aa      = df_pq_index(adf,adf) 
      ii      = df_pq_index(idf,idf) 
      ai      = df_pq_index(adf,idf) 

      ! 3 * g(ai|ai)
      val     = 3.0_wp * my_ddot(df_vars_%nQ,int2(ai+1:),df_vars_%Qstride,int2(ai+1:),df_vars_%Qstride)

      ! - g(aa|ii)
      val     = val - my_ddot(df_vars_%nQ,int2(aa+1:),df_vars_%Qstride,int2(ii+1:),df_vars_%Qstride)

      ! factor of 4 from overall formula
 
      te_terms_ed_df = 4.0_wp * val

      return

    end function te_terms_ed_df

    integer function diagonal_hessian_ea(f_i_ext,f_a_ext,q,z,int2,den1,den2)
      implicit none
      ! subroutine to compute the diagonal Hessian matrix element H(ai|ai) according to Eq. 4.7b of
      ! Chaban, Schmidt, Gordon, Theor. Chem. Acc., 97, 88-95 (1997) 
      ! H(at|at) = 2 * [  d(t,t) * ( f_i(a,a) + f_a(a,a) ) - q(t,t) - z(t,t) ] 
      ! where z(t,t) = sum_{u} [ d(t,u) * f_i(t,u) ]  &&  a \in ext & t,u \in act 
      real(wp), intent(in) :: z(:,:),q(:,:),int2(:),den1(:),den2(:)
      type(vector_block), intent(in) :: f_i_ext(:),f_a_ext(:)

      integer :: a,t,grad_ind,tt_den,a_sym,a_i
      real(wp) :: h_val

      ! loop over all external - doubly-occupied pairs

      do a_sym = 1 , nirrep_

        grad_ind = rot_pair_%pair_offset(a_sym,rot_pair_%ext_act_type)

        do t = first_index_(a_sym,2) , last_index_(a_sym,2)

          do a = first_index_(a_sym,3) , last_index_(a_sym,3)

            ! integral/density indeces
  
            tt_den  = dens_%gemind(t,t)

            a_i     = trans_%class_to_irrep_map(a)-ndocpi_(a_sym)-nactpi_(a_sym)

            ! update diagonal Hessian index
  
            grad_ind = grad_ind + 1

            ! calculate diagonal Hessian element

            h_val    = 2.0_wp * ( den1(tt_den) * f_i_ext(a_sym)%val(a_i) - q(t-ndoc_tot_,t) - z(t-ndoc_tot_,t) )

            if ( use_exact_hessian_diagonal_ == 0 ) then

              h_val = h_val + 2.0_wp * den1(tt_den) * f_a_ext(a_sym)%val(a_i)

            else

              if ( df_vars_%use_df_teints == 1 ) then

                h_val = h_val + 2.0_wp * den1(tt_den) * f_a_ext(a_sym)%val(a_i) 

              else

                h_val = h_val + te_terms_ea(t,a,a_sym,int2,den2)

              endif

            endif

            if ( h_val < 0.0_wp ) h_val = - h_val

            if ( h_val < min_diag_hessian_ ) h_val = min_diag_hessian_

            orbital_hessian_(grad_ind) = h_val

          end do

        end do

      end do

      diagonal_hessian_ea = 0

      return

    end function diagonal_hessian_ea

    function te_terms_ea(t,a,a_sym,int2,den2)

      ! function to compute 2-e approximation to the active-external digaonal Hessian
      ! using 4-index integrals according to
      ! te = 2* sum{u,v} [ d(tt|uv) * g(aa|uv) + 2 * d(tv|tu) * g(av|au) ] 

      integer, intent(in)  :: t,a,a_sym
      real(wp), intent(in) :: int2(:),den2(:)

      integer              :: u,v,u_sym,au_sym
      integer              :: tt_den,uv_den,tv_den,tu_den,uu_den
      integer              :: aa_int,uv_int,av_int,au_int,uu_int
      integer              :: den_offset,int_offset,int_ind,den_ind

      real(wp)             :: te_terms_ea,val

      val    = 0.0_wp

      tt_den = dens_%gemind(t,t)
      aa_int = ints_%gemind(a,a)

      do u_sym = 1 , nirrep_

        au_sym     = group_mult_tab_(a_sym,u_sym)

        den_offset = dens_%offset(au_sym)
        int_offset = ints_%offset(au_sym)
 
        do u = first_index_(u_sym,2) , last_index_(u_sym,2)

          tu_den = dens_%gemind(t,u)
          uu_den = dens_%gemind(u,u)
          au_int = ints_%gemind(a,u)
          uu_int = ints_%gemind(u,u)

          ! u > v --> factor of 2

          do v = first_index_(u_sym,2) , u - 1

            uv_den  = dens_%gemind(u,v)
            tv_den  = dens_%gemind(t,v)
            uv_int  = ints_%gemind(u,v)
            av_int  = ints_%gemind(a,v)

            ! d(tt|uv) * g(aa|uv)
            int_ind = pq_index(aa_int,uv_int) 
            den_ind = pq_index(tt_den,uv_den)
            val     = val + den2(den_ind) * int2(int_ind)

            ! 2 * d(tu|tv) * g(au|av)
            int_ind = pq_index(au_int,av_int) + int_offset
            den_ind = pq_index(tu_den,tv_den) + den_offset
            val     = val + 2.0_wp * den2(den_ind) * int2(int_ind)

          end do

          ! u = v --> factor of 1
 
          ! d(tt|uu) * g(aa|uu)
          int_ind = pq_index(aa_int,uu_int)
          den_ind = pq_index(tt_den,uu_den)
          val     = val + 0.5_wp * den2(den_ind) * int2(int_ind)

          ! 2 * d(tu|tu) * g(au|au)
          int_ind = pq_index(au_int,au_int) + int_offset
          den_ind = pq_index(tu_den,tu_den) + den_offset
          val     = val + den2(den_ind) * int2(int_ind)

        end do

      end do

      ! factor of 2 because only LT elements are accessed 
      ! factor of 2 from the overall formula

      te_terms_ea = 4.0_wp * val

      return

    end function te_terms_ea

    function te_terms_ea_df(t,a,a_sym,int2,den2)

      integer, intent(in)  :: t,a,a_sym
      real(wp), intent(in) :: int2(:),den2(:)

      integer              :: u,v,u_sym,au_sym
      integer              :: tt_den,uv_den,tv_den,tu_den,uu_den
      integer              :: aa,uv,av,au,uu,nQ
      integer              :: adf,udf,vdf,tdf
      integer              :: den_offset,den_ind

      real(wp)             :: te_terms_ea_df,val
      real(wp)             :: v_aa(df_vars_%nQ),v_au(df_vars_%nQ)

      val    = 0.0_wp

      nQ     = int(df_vars_%nQ,kind=ip)

      tt_den = dens_%gemind(t,t)

      adf    = df_vars_%class_to_df_map(a)
      tdf    = df_vars_%class_to_df_map(t)

      aa     = df_pq_index(adf,adf) 

      call my_dcopy(df_vars_%nQ,int2(aa+1:),df_vars_%Qstride,v_aa,1)

      do u_sym = 1 , nirrep_

        au_sym     = group_mult_tab_(a_sym,u_sym)

        den_offset = dens_%offset(au_sym)

        do u = first_index_(u_sym,2) , last_index_(u_sym,2)

          tu_den = dens_%gemind(t,u)
          uu_den = dens_%gemind(u,u)

          udf    = df_vars_%class_to_df_map(u)

          au     = df_pq_index(adf,udf) 
          uu     = df_pq_index(udf,udf) 

          call my_dcopy(df_vars_%nQ,int2(au+1:),df_vars_%Qstride,v_au,1)

          ! u > v --> factor of 2

          do v = first_index_(u_sym,2) , u - 1

            uv_den  = dens_%gemind(u,v)
            tv_den  = dens_%gemind(t,v)

            vdf     = df_vars_%clasS_to_df_map(v)

            uv      = df_pq_index(udf,vdf) 
            av      = df_pq_index(adf,vdf) 

            ! d(tt|uv) * g(aa|uv)
            den_ind = pq_index(tt_den,uv_den)
            val     = val + den2(den_ind) * my_ddot(df_vars_%nQ,v_aa,1,int2(uv+1:),df_vars_%Qstride)

            ! 2 * d(tu|tv) * g(au|av)
            den_ind = pq_index(tu_den,tv_den) + den_offset
            val     = val + 2.0_wp * den2(den_ind) * my_ddot(df_vars_%nQ,v_au,1,int2(av+1:),df_vars_%Qstride)

          end do

          ! u = v --> factor of 1

          ! d(tt|uu) * g(aa|uu)
          den_ind = pq_index(tt_den,uu_den)
          val     = val + 0.5_wp * den2(den_ind) * my_ddot(df_vars_%nQ,v_aa,1,int2(uu+1:),df_vars_%Qstride)

          ! 2 * d(tu|tu) * g(au|au)
          den_ind = pq_index(tu_den,tu_den) + den_offset
          val     = val + den2(den_ind) * my_ddot(df_vars_%nQ,v_aa,1,v_au,1)

        end do

      end do

      ! factor of 2 because only LT elements are accessed 
      ! factor of 2 from the overall formula

      te_terms_ea_df = 4.0_wp * val

      return

    end function te_terms_ea_df

    subroutine allocate_hessian_data()
      implicit none
      if (allocated(orbital_hessian_)) call deallocate_hessian_data()

      allocate(orbital_hessian_(rot_pair_%n_tot))

      orbital_hessian_ = 0.0_wp

      return
    end subroutine allocate_hessian_data

    subroutine deallocate_hessian_data
      implicit none

      if (allocated(orbital_hessian_)) deallocate(orbital_hessian_)
      return
    end subroutine deallocate_hessian_data

  subroutine print_orbital_hessian()
    implicit none
    integer :: ij_pair,i,j,newline,i_class,j_class,j_class_start,j_start,i_sym
    character(1) :: ityp,jtyp
    ! loop over rotation pairs

    write(fid_,'(a)')'orbital hessian:'

    newline=0

    ! the loop structure below cycles through the possible rotation pairs
    ! rotation pairs are sorted according to symmetry and for each irrep,
    ! the rotation pairs are sorted according to orbital classes: ad,ed,aa,ea
    ! for each pair j>i; for more details, see subroutine setup_rotation_indeces in focas_main.F90

    ij_pair = 0

    do i_sym = 1 , nirrep_

      do i_class = 1 , 3

        j_class_start = i_class + 1

        if ( ( include_aa_rot_ == 1 ) .and. ( i_class == 2 ) ) j_class_start = i_class

        do j_class = j_class_start , 3

          do i = first_index_(i_sym,i_class) , last_index_(i_sym,i_class)

            j_start = first_index_(i_sym,j_class)

            if ( i_class == j_class ) j_start = i + 1

            do j = j_start , last_index_(i_sym,j_class)

              if ( i_class == 1 ) then

                ityp = 'd'

                jtyp = 'a'

                if ( j_class == 3 ) jtyp = 'e'

              else

                ityp = 'a'

                jtyp = 'a'

                if ( j_class == 3 ) jtyp = 'e'

              end if

              ij_pair = ij_pair + 1

              write(fid_,'(a,a,a,a,i3,a,i3,a,1x,es10.3,4x)',advance='no')jtyp,'-',ityp,' (',j,',',i,')',orbital_hessian_(ij_pair)

              newline = newline + 1

              if ( mod(newline,4) == 0 ) write(fid_,*)

            end do

          end do

        end do

      end do

    end do
    if ( mod(newline,4) /= 0 ) write(fid_,*)

    return

  end subroutine print_orbital_hessian

end module focas_hessian
