!!
 !@BEGIN LICENSE
 !
 ! v2RDM-CASSCF, a plugin to:
 !
 ! Psi4: an open-source quantum chemistry software package
 !
 ! This program is free software; you can redistribute it and/or modify
 ! it under the terms of the GNU General Public License as published by
 ! the Free Software Foundation; either version 2 of the License, or
 ! (at your option) any later version.
 !
 ! This program is distributed in the hope that it will be useful,
 ! but WITHOUT ANY WARRANTY; without even the implied warranty of
 ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 ! GNU General Public License for more details.
 !
 ! You should have received a copy of the GNU General Public License along
 ! with this program; if not, write to the Free Software Foundation, Inc.,
 ! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 !
 !@END LICENSE
 !
 !!

module focas_gradient
  use focas_data
 
  implicit none

  contains

  subroutine orbital_gradient(int1,int2,den1,den2)
    implicit none
    real(wp), intent(in) :: int1(:),int2(:),den1(:),den2(:)
    real(wp) :: t0,t1
    type(fock_info) :: fock
    integer :: i
    real(wp), allocatable :: tq(:,:)
   
    ! calculate inactive Fock matrix
    if ( df_vars_%use_df_teints == 0 ) then
       call compute_f_i(int1,int2)
    else
      call compute_f_i_df_coulomb(int1,int2)
      call compute_f_i_df_exchange_fast(int2)
!      call compute_f_i_df_exchange(int2)
    endif
    call transpose_matrix(fock_i_)

    ! calculate active Fock matrix
    if ( df_vars_%use_df_teints == 0 ) then
      call compute_f_a(den1,int2)
    else
      call compute_f_a_df_coulomb(den1,int2)
!      call compute_f_a_df_exchange(den1,int2)
      call compute_f_a_df_exchange_fast(den1,int2)
    endif
    call transpose_matrix(fock_a_)

    ! calculate auxiliary q matrix
    if ( df_vars_%use_df_teints == 0 ) then
      call compute_q(den2,int2)
    else
      call compute_q_df(den2,int2)
    endif

    ! calculate auxiliary z matrix
    call compute_z(den1)

    ! compute gradient
    call compute_orbital_gradient()

!    write(fid_,*)
!    write(fid_,*)'F_i matrix'
!    call print_f_matrix(fock_i_)
!
!    write(fid_,*)
!    write(fid_,*)'F_a matrix'
!    call print_f_matrix(fock_a_)
!
!    write(fid_,*)
!    write(fid_,*)'Q matrix'
!    call print_q_matrix(q_)
!
!    write(fid_,*)'orbital gradient vector'
!    call print_vector(orbital_gradient_)

    ! determine value,type, and orbital indices for largest gradient element
    call check_max_gradient()

    return
  end subroutine orbital_gradient

  subroutine transpose_matrix(fock)

    implicit none
    
    type(fock_info) :: fock

    integer :: num_p,p,q,p_sym

    ! subroutine to symmetrize the inactive-inactive and active-active 
    ! blocks of the Fock matrix
        
    do p_sym = 1 , nirrep_

      num_p = ndocpi_(p_sym) + nactpi_(p_sym)

      if ( num_p == 0 ) cycle

      do p = 1 , num_p

        do q = p+1 , num_p

          fock%occ(p_sym)%val(p,q) = fock%occ(p_sym)%val(q,p)

        end do

      end do

    end do

  end subroutine transpose_matrix

  subroutine check_max_gradient()

    implicit none

    integer  :: grad_ind,a_sym,t_sym,i,a,t,u
    real(wp) :: abs_grad_val,grad_tol

    ! ********************************
    ! external - doubly-occupied pairs
    ! ********************************

    max_grad_val_ = 0.0_wp
    max_grad_ind_ = 0
    max_grad_typ_ = 0
    max_grad_sym_ = 0

    do a_sym = 1 , nirrep_

      grad_ind = rot_pair_%pair_offset(a_sym,rot_pair_%ext_doc_type)

      do i = first_index_(a_sym,1) , last_index_(a_sym,1)

        do a = first_index_(a_sym,3) , last_index_(a_sym,3)

          grad_ind = grad_ind + 1

          abs_grad_val = abs(orbital_gradient_(grad_ind))

          if ( abs_grad_val < max_grad_val_ ) cycle

          max_grad_val_    = abs_grad_val 

          max_grad_ind_(1) = a
          max_grad_ind_(2) = i

          max_grad_typ_    = 2

          max_grad_sym_    = a_sym

        end do

      end do

    end do

    ! ******************************
    ! acitve - doubly-occupied pairs
    ! ******************************

    do t_sym = 1 , nirrep_

      grad_ind = rot_pair_%pair_offset(t_sym,rot_pair_%act_doc_type)

      do i = first_index_(t_sym,1) , last_index_(t_sym,1)

        do t = first_index_(t_sym,2) , last_index_(t_sym,2)

          grad_ind = grad_ind + 1

          abs_grad_val = abs(orbital_gradient_(grad_ind))

          if ( abs_grad_val < max_grad_val_ ) cycle

          max_grad_val_    = abs_grad_val

          max_grad_ind_(1) = t
          max_grad_ind_(2) = i

          max_grad_typ_    = 1

          max_grad_sym_    = t_sym

        end do

      end do

    end do

    ! ***********************
    ! external - active pairs
    ! ***********************

    do a_sym = 1 , nirrep_

      grad_ind = rot_pair_%pair_offset(a_sym,rot_pair_%ext_act_type)

      do t = first_index_(a_sym,2) , last_index_(a_sym,2)

        do a = first_index_(a_sym,3) , last_index_(a_sym,3)

          grad_ind = grad_ind + 1

          abs_grad_val = abs(orbital_gradient_(grad_ind))

          if ( abs_grad_val < max_grad_val_ ) cycle

          max_grad_val_    = abs_grad_val

          max_grad_ind_(1) = a
          max_grad_ind_(2) = t

          max_grad_typ_    = 4

          max_grad_sym_    = a_sym

        end do

      end do

    end do

    ! *********************
    ! active - active pairs
    ! *********************

    if ( include_aa_rot_ == 1 ) then

      do t_sym = 1 , nirrep_

        grad_ind = rot_pair_%pair_offset(t_sym,rot_pair_%act_act_type)

        do u = first_index_(t_sym,2) , last_index_(t_sym,2)

          do t = u + 1 , last_index_(t_sym,2)

            grad_ind = grad_ind + 1

            abs_grad_val = abs(orbital_gradient_(grad_ind))

            if ( abs_grad_val < max_grad_val_ ) cycle

            max_grad_val_    = abs_grad_val

            max_grad_ind_(1) = t
            max_grad_ind_(2) = u

            max_grad_typ_    = 3

            max_grad_sym_    = t_sym

          end do

        end do

      end do

    end if

    return

  end subroutine check_max_gradient

  subroutine zero_frozen_docc_vector_elements()

    implicit none
 
    integer  :: a,i,t
    integer  :: a_i,i_i,t_i
    integer  :: a_sym,t_sym
    integer  :: grad_ind

    ! ********************************
    ! external - doubly-occupied pairs
    ! ********************************

    do a_sym = 1 , nirrep_

      grad_ind = rot_pair_%pair_offset(a_sym,rot_pair_%ext_doc_type)

      do i = first_index_(a_sym,1) , last_index_(a_sym,1)

        do a = first_index_(a_sym,3) , last_index_(a_sym,3)

          grad_ind = grad_ind + 1

          i_i      = trans_%class_to_irrep_map(i)
          a_i      = trans_%class_to_irrep_map(a)

          if ( i_i > nfzcpi_(a_sym) ) cycle

          orbital_gradient_(grad_ind) = 0.0_wp

        end do

      end do

    end do

    ! ******************************
    ! acitve - doubly-occupied pairs
    ! ******************************

    do t_sym = 1 , nirrep_

      grad_ind = rot_pair_%pair_offset(t_sym,rot_pair_%act_doc_type)

      do i = first_index_(t_sym,1) , last_index_(t_sym,1)

        do t = first_index_(t_sym,2) , last_index_(t_sym,2)

          grad_ind = grad_ind + 1

          i_i      = trans_%class_to_irrep_map(i)
          t_i      = trans_%class_to_irrep_map(t)

          if ( i_i > nfzcpi_(t_sym) ) cycle

          orbital_gradient_(grad_ind) = 0.0_wp

        end do

      end do

    end do


  end subroutine zero_frozen_docc_vector_elements

  subroutine compute_orbital_gradient()

    ! subroutine to compute the orbital gradient without explicit storage
    ! of the generalized Fock matrix

    implicit none

    integer :: grad_ind
    integer :: a_sym,t_sym
    integer :: i,t,u,a
    integer :: ia,it,i_i,a_i,t_i

    orbital_gradient_ = 0.0_wp

    ! ********************************
    ! external - doubly-occupied pairs
    ! ********************************

    do a_sym = 1 , nirrep_

      grad_ind = rot_pair_%pair_offset(a_sym,rot_pair_%ext_doc_type)

      do i = first_index_(a_sym,1) , last_index_(a_sym,1)

        do a = first_index_(a_sym,3) , last_index_(a_sym,3)
  
          grad_ind = grad_ind + 1

          i_i      = trans_%class_to_irrep_map(i)
          a_i      = trans_%class_to_irrep_map(a)

          ia       = ints_%gemind(i,a)

          orbital_gradient_(grad_ind) = 4.0_wp*( fock_i_%occ(a_sym)%val(a_i,i_i) &
                                       &       + fock_a_%occ(a_sym)%val(a_i,i_i) )

        end do

      end do
 
    end do  

    ! ******************************
    ! acitve - doubly-occupied pairs
    ! ******************************

    do t_sym = 1 , nirrep_

      grad_ind = rot_pair_%pair_offset(t_sym,rot_pair_%act_doc_type)

      do i = first_index_(t_sym,1) , last_index_(t_sym,1)

        do t = first_index_(t_sym,2) , last_index_(t_sym,2)

          grad_ind = grad_ind + 1

          i_i      = trans_%class_to_irrep_map(i)
          t_i      = trans_%class_to_irrep_map(t)

          it       = ints_%gemind(i,t)

          orbital_gradient_(grad_ind) = 4.0_wp*( fock_i_%occ(t_sym)%val(t_i,i_i)    & 
                                      &  +       fock_a_%occ(t_sym)%val(t_i,i_i) )  &
                                      & - 2.0_wp * ( q_(t - ndoc_tot_,i) + z_(t - ndoc_tot_,i))

        end do

      end do

    end do

    ! ***********************
    ! external - active pairs
    ! ***********************

    do a_sym = 1 , nirrep_

      grad_ind = rot_pair_%pair_offset(a_sym,rot_pair_%ext_act_type)

      do t = first_index_(a_sym,2) , last_index_(a_sym,2)

        do a = first_index_(a_sym,3) , last_index_(a_sym,3)

          grad_ind = grad_ind + 1

          orbital_gradient_(grad_ind) = 2.0_wp * ( q_(t - ndoc_tot_,a) + z_(t - ndoc_tot_,a))

        end do
    
      end do

    end do

    ! *********************
    ! active - active pairs
    ! *********************

    if ( include_aa_rot_ == 1 ) then

      do t_sym = 1 , nirrep_

        grad_ind = rot_pair_%pair_offset(t_sym,rot_pair_%act_act_type)

        do u = first_index_(t_sym,2) , last_index_(t_sym,2)

          do t = u + 1 , last_index_(t_sym,2)

            grad_ind = grad_ind + 1

            orbital_gradient_(grad_ind) = 2.0_wp * ( q_(u - ndoc_tot_,t) +         &
                 & z_(u - ndoc_tot_,t) - q_(t - ndoc_tot_,u) - z_(t - ndoc_tot_,u))

          end do
 
        end do

      end do

    end if

    if ( nfzc_tot_ > 0 ) call zero_frozen_docc_vector_elements()

    grad_norm_ = sqrt(my_ddot(rot_pair_%n_tot,orbital_gradient_,1,orbital_gradient_,1))

    return 

  end subroutine compute_orbital_gradient

  subroutine compute_z(den1)
    implicit none
    ! function to compute contraction of the density with the inactive Fock matrix according to
    ! z(m,t) = sum_u { den1(tu) * f_i(mu) } t,u \in A m \in D,A,E
    real(wp), intent(in) :: den1(:)
    integer :: t,u,tu_den,mu_int,t_sym,m,m_class,m_i,u_i
    real(wp) :: val

    ! initialize output

    z_ = 0.0_wp

    ! loop over symmetries for t

    do t_sym = 1 , nirrep_

      ! loop over u indeces

      do t = first_index_(t_sym,2) , last_index_(t_sym,2)
 
        ! m_class < u_class ==> m < u

        ! loop over m orbital classes

        do m_class = 1 , 3
 
          do m = first_index_(t_sym,m_class) , last_index_(t_sym,m_class)

            m_i = trans_%class_to_irrep_map(m)

            ! initialize contraction value

            val = 0.0_wp

            do u = first_index_(t_sym,2) , last_index_(t_sym,2) 

              ! integral/density addressing

              tu_den = dens_%gemind(t,u)
              u_i    = trans_%class_to_irrep_map(u)

              val = val + den1(tu_den) * fock_i_%occ(t_sym)%val(m_i,u_i)

            end do ! end u loop

            z_( t - ndoc_tot_ , m ) = val

          end do ! end m loop

        end do ! end m_class loop 

      end do ! end t loop

    end do ! end t_sym loop

    return

  end subroutine compute_z

  subroutine compute_q_df(den2,int2)

    implicit none

    real(wp), intent(in) :: den2(:),int2(:)

    integer :: tu_sym,t_sym,u_sym,v_sym,w_sym,p_sym
    integer :: p_class
    integer :: t,u,v,w,p,tu_den,vw_den,tuvw
    integer :: vdf,wdf,pdf,udf,den_off,den_ind,tu_int
    integer(ip) :: vw_df,pu_df,nQ
  
    real(wp) :: val
 
    nQ = int ( df_vars_%nQ , kind = ip )

    ! initialize
 
    q_ = 0.0_wp
 
    do tu_sym = 1 , nirrep_

      if ( dens_%ngempi(tu_sym) == 0 ) cycle

      qint_%tuQ(tu_sym)%val = 0.0_wp

    end do

    ! assemble intermediates

    tu_sym = 1

    do t_sym = 1 , nirrep_

      do t = first_index_(t_sym,2) , last_index_(t_sym,2)

!$omp parallel shared(t_sym,t,tu_sym,first_index_,last_index_, &
!$omp int2,den2,df_vars_,dens_,qint_) num_threads(nthread_use_)
!$omp do private(v,vdf,w,wdf,vw_df,vw_den,den_ind,u,tu_den,v_sym)

        do u = first_index_(t_sym,2) , t

          tu_den = dens_%gemind(t,u)

          do v_sym = 1 , nirrep_

            do v = first_index_(v_sym,2) , last_index_(v_sym,2)

              vdf = df_vars_%class_to_df_map(v)

              do w = first_index_(v_sym,2) , v - 1

                wdf = df_vars_%class_to_df_map(w)

                vw_df =  df_pq_index(vdf,wdf) 

                vw_den = dens_%gemind(v,w)

                den_ind = pq_index(tu_den,vw_den) 
                ! update array

                call my_daxpy(df_vars_%nQ,2.0_wp*den2(den_ind),int2(vw_df+1:),df_vars_%Qstride,qint_%tuQ(tu_sym)%val(:,tu_den),1)

              end do ! end w loop

              vw_df =  df_pq_index(vdf,vdf) 

              vw_den = dens_%gemind(v,v)

              den_ind = pq_index(tu_den,vw_den)

              ! update array

              call my_daxpy(df_vars_%nQ,den2(den_ind),int2(vw_df+1:),df_vars_%Qstride,qint_%tuQ(tu_sym)%val(:,tu_den),1) 

            end do ! end v loop

          end do ! end v_sym loop

        end do ! end u loop

!$omp end do nowait
!$omp end parallel

      end do ! end t loop

    end do ! end t_sym loop

    do tu_sym = 2 , nirrep_
 
      den_off = dens_%offset(tu_sym)
 
      do t_sym = 1 , nirrep_

        u_sym=group_mult_tab_(tu_sym,t_sym)
         
        if ( u_sym > t_sym ) cycle

        do t = first_index_(t_sym,2) , last_index_(t_sym,2)

!$omp parallel shared(t_sym,u_sym,tu_sym,t,den_off,      &
!$omp first_index_,last_index_,int2,den2,df_vars_,dens_, &
!$omp qint_) num_threads(nthread_use_)
!$omp do private(v,vdf,w,wdf,vw_df,vw_den,den_ind,tu_den,u,v_sym,w_sym)

          do u = first_index_(u_sym,2) , last_index_(u_sym,2)

            tu_den  = dens_%gemind(t,u)

            do v_sym = 1 , nirrep_

              w_sym = group_mult_tab_(tu_sym,v_sym)

              if ( w_sym > v_sym ) cycle

              do v = first_index_(v_sym,2) , last_index_(v_sym,2)

                vdf = df_vars_%class_to_df_map(v)

                do w = first_index_(w_sym,2) , last_index_(w_sym,2)

                  wdf = df_vars_%class_to_df_map(w)

                  vw_df =  df_pq_index(vdf,wdf) 

                  vw_den = dens_%gemind(v,w)

                  den_ind = pq_index(tu_den,vw_den) + den_off

                  ! update array

                  call my_daxpy(df_vars_%nQ,2.0_wp*den2(den_ind),int2(vw_df+1:),df_vars_%Qstride,qint_%tuQ(tu_sym)%val(:,tu_den),1)

                end do ! end w loop

              end do ! end v loop

            end do ! end v_sym loop

          end do ! end u_loop

!$omp end do nowait
!$omp end parallel

        end do ! end t_loop

      end do ! end t_sym loop
  
    end do ! end tu_sym loop

    ! compute Q-elements

    do p_sym = 1 , nirrep_

      do p_class = 1 , 3

        do p = first_index_(p_sym,p_class) , last_index_(p_sym,p_class)

          pdf = df_vars_%class_to_df_map(p)

          do t = first_index_(p_sym,2) , last_index_(p_sym,2)

            val = 0.0_wp

            do tu_sym = 1 , nirrep_

              u_sym = group_mult_tab_(tu_sym,p_sym)

              do u = first_index_(u_sym,2) , last_index_(u_sym,2)

                tu_den = dens_%gemind(t,u)

                udf = df_vars_%class_to_df_map(u)

                pu_df = df_pq_index(pdf,udf) 

                val   = val + my_ddot(df_vars_%nQ,int2(pu_df+1:),df_vars_%Qstride,qint_%tuQ(tu_sym)%val(:,tu_den),1)

              end do ! end u loop

            end do ! end tu_sym loop

            q_(t - ndoc_tot_ , p ) = val

          end do ! end t loop

        end do ! end p loop

      end do ! and p_class loop

    end do ! end p_sym_loop

    return
  
  end subroutine compute_q_df

  subroutine compute_q(den2,int2)
    implicit none
    ! this subroutine computes the "auxiliary Q" matrix according to Eq. 12.5.14 in Helgaker on page 622
    ! Q(v,m) = \SUM[w,x,y \in A] { d2(vw|xy) * g(mw|xy) }
    ! Because only those D2 elements with all four indeces \in A are nonzero, we have v \in A
    ! Also, m is a general orbital index && sym_m == sym_v
    real(wp), intent(in) :: den2(:),int2(:)   
    integer :: mw_sym,x_sym,y_sym,w_sym,m_sym,m_class,m,v,w,x,y
    integer :: xy_den,xy_int,mw,vw
    integer :: den_ind,den_sym_offset,int_ind,int_sym_offset
    real(wp) :: val

    ! initialize
    q_ = 0.0_wp

    ! loop irreps for m and v

    do m_sym = 1 , nirrep_
 
      ! loop over irreps for w

      do w_sym = 1, nirrep_

        ! mw//vw symmetry
        mw_sym = group_mult_tab_(m_sym,w_sym)
        ! offsets for integral/density addressing
        den_sym_offset = dens_%offset(mw_sym)
        int_sym_offset = ints_%offset(mw_sym)

        ! loop over irreps for x

        do x_sym = 1 , nirrep_

          ! correspoing irrep for y
          y_sym = group_mult_tab_(x_sym,mw_sym)

          ! at this point, we have mw_sym == xy_sym && m_sym == v_sym
          ! loop over m_class

          do m_class = 1 , 3

            ! loop over v indeces
   
            do m = first_index_(m_sym,m_class) , last_index_(m_sym,m_class)

              ! loop over v \in A

              do v = first_index_(m_sym,2) , last_index_(m_sym,2)

                ! initialize q matrix element
                val     = 0.0_wp

                ! loop over w indeces

                do w = first_index_(w_sym,2) , last_index_(w_sym,2)
                  
                  ! save geminal indeces for integral/density addressing
                  mw = ints_%gemind(m,w)
                  vw = dens_%gemind(v,w)

                  ! loop over x indeces
 
                  do x = first_index_(x_sym,2) , last_index_(x_sym,2)

                    ! loop over y indeces

                    do y = first_index_(y_sym,2) , last_index_(y_sym,2) 

                      ! save geminal indeces for integral/density addressing
                      xy_int  = ints_%gemind(x,y)
                      xy_den  = dens_%gemind(x,y)
                      ! integral/density addresses
                      int_ind = pq_index(xy_int,mw) + int_sym_offset                      
                      den_ind = pq_index(xy_den,vw) + den_sym_offset
                      ! update temporary value
                      val     = val + den2(den_ind) * int2(int_ind)
 
                    end do ! end y loop

                  end do ! end x loop
                   
                end do ! end w loop

                ! update q matrix element
                q_( v - ndoc_tot_ , m ) = q_( v - ndoc_tot_ , m ) + val

              end do ! end v loop

            end do ! end m loop 

          end do ! end m_class loop

        end do ! end x_sym loop

      end do ! end w_sym loop

    end do ! end m_sym loop

    return
  end subroutine compute_q

  subroutine compute_f_a_df_exchange_fast(den1,int2)

    implicit none
    ! subroutine to compute exchange contributions to the active Fock matrix

    ! Fa(p,q) = Fa(p,q) - 0.5_wp * SUM_[t,u \in A] (pt|qu) * D1(t|u)

    ! The subroutine performs the comptations in 3 steps
    !    1) collect ALL three-index integrals and sort according to p_class,p_sym,t_sym
    !    2) evaluate required 4-index integrals using DGEMMs 
    !    3) update Fock matrix elements

    real(wp), intent(in) :: int2(:)
    real(wp), intent(in) :: den1(:)

    integer     :: p_class,q_class
    integer     :: p_sym,t_sym
    integer     :: p_max
    integer     :: p,q,t,u
    integer     :: p_i,q_i,t_i,u_i
    integer     :: pt_s,qu_s
    integer     :: p_off,t_off
    integer     :: pdf,tdf
    integer     :: nmo_p,nmo_t,nmo_q
    integer     :: pt_ind,qu_ind,tu_den
    integer(ip) :: pt
    real(wp)    :: d_val

    ! allocate temporary scratch arrays to hold three-index integrals
    call allocate_fock_scr(fa_scr_,fa_scr_4index_,nactpi_,2)

    ! **************************************
    ! gather three-index integrals integrals
    ! **************************************

    do p_class = 1 , 3

      do p_sym = 1 , nirrep_

        nmo_p = last_index_(p_sym,p_class)-first_index_(p_sym,p_class)+1

        if ( nmo_p == 0 ) cycle

        do t_sym = 1 , nirrep_

          nmo_t = nactpi_(t_sym)

          if ( nmo_t == 0 ) cycle

          pt_ind = 0

          do t = first_index_(t_sym,2) , last_index_(t_sym,2)

            tdf = df_vars_%class_to_df_map(t)

            do p = first_index_(p_sym,p_class) , last_index_(p_sym,p_class)

              pdf = df_vars_%class_to_df_map(p)

              pt = df_pq_index(pdf,tdf)

              pt_ind = pt_ind + 1

              call my_dcopy(df_vars_%nQ,int2(pt+1:),df_vars_%Qstride, &
              & fa_scr_(p_class)%p_sym(p_sym)%t_sym(t_sym)%val(:,pt_ind),1)

            end do ! t loop

          end do ! p loop

        end do ! t_sym loop

      end do ! p_sym loop

    end do ! p_class loop

    ! ****************************************************
    ! off diagonal elements (p_class > q_class --> p > q )
    ! ****************************************************

    do p_class = 1 , 3

      do q_class = 1 , p_class - 1

        do p_sym = 1 , nirrep_

          nmo_p = last_index_(p_sym,p_class) - first_index_(p_sym,p_class) + 1

          if ( nmo_p == 0 ) cycle

          nmo_q = last_index_(p_sym,q_class) - first_index_(p_sym,q_class) + 1

          if ( nmo_q == 0 ) cycle

          fa_scr_4index_ = 0.0_wp

          do t_sym = 1 , nirrep_

            nmo_t    = nactpi_(t_sym)

            if ( nmo_t == 0 ) cycle

            t_off = ndocpi_(t_sym)

            do u = first_index_(t_sym,2) , last_index_(t_sym,2)

              qu_s = ( trans_%class_to_irrep_map(u) - 1 -t_off ) * nmo_q

              do t = first_index_(t_sym,2) , last_index_(t_sym,2)

                pt_s = ( trans_%class_to_irrep_map(t) - 1 -t_off ) * nmo_p

                tu_den = dens_%gemind(t,u)

                d_val = 0.5_wp * den1(tu_den)

                ! *********************************************************************
                ! calculate ALL 4index integrals for this p_class,p_sym, t_sym, t and u
                ! (pt|qu) = fa_scr_4index_(pt_ind,qu_ind)
                ! *********************************************************************

                call dgemm('t','n',nmo_p,nmo_q,df_vars_%nQ,d_val,                  &
                     & fa_scr_(p_class)%p_sym(p_sym)%t_sym(t_sym)%val(:,pt_s+1:pt_s+nmo_p),df_vars_%nQ, &
                     & fa_scr_(q_class)%p_sym(p_sym)%t_sym(t_sym)%val(:,qu_s+1:qu_s+nmo_q),df_vars_%nQ, &
                     & 1.0_wp,fa_scr_4index_(1:nmo_p,1:nmo_q),nmo_p)

              end do ! t loop

            end do ! u loop
 
          end do ! t_sym loop

          qu_ind = 0

          do q = first_index_(p_sym,q_class) , last_index_(p_sym,q_class)

            qu_ind = qu_ind + 1

            q_i = trans_%class_to_irrep_map(q)

            pt_ind = 0

            do p = first_index_(p_sym,p_class) , last_index_(p_sym,p_class)

              p_i = trans_%class_to_irrep_map(p)

              pt_ind = pt_ind + 1

              fock_a_%occ(p_sym)%val(p_i,q_i) = fock_a_%occ(p_sym)%val(p_i,q_i) - &
                   &  fa_scr_4index_(pt_ind,qu_ind)

            end do ! p loop

          end do ! q loop

        end do ! p_sym

      end do ! q_class

    end do ! p_class

    ! *******************************************************
    ! diagonal elements (p_class = q_class // enforce p > q )
    ! doubly-occupied and active orbitals (need all elements) 
    ! *******************************************************

    do p_class = 1 , 2

      do p_sym = 1 , nirrep_

        p_off = 0

        if ( p_class == 2 ) p_off = ndocpi_(p_sym)

        nmo_p = last_index_(p_sym,p_class) - first_index_(p_sym,p_class) + 1

        if ( nmo_p == 0 ) cycle

        fa_scr_4index_ = 0.0_wp

        do t_sym = 1 , nirrep_

          nmo_t = nactpi_(t_sym)

          if ( nmo_t == 0 ) cycle

          t_off = ndocpi_(t_sym)

          do u = first_index_(t_sym,2) , last_index_(t_sym,2)

            qu_s = ( trans_%class_to_irrep_map(u) - 1 - t_off ) * nmo_p

            do t = first_index_(t_sym,2) , last_index_(t_sym,2)

              pt_s = ( trans_%class_to_irrep_map(t) - 1 - t_off ) * nmo_p

              tu_den = dens_%gemind(t,u)

              d_val = 0.5_wp * den1(tu_den)

              ! *********************************************************************
              ! calculate ALL 4index integrals for this p_class,p_sym, t_sym, t and u
              ! (pt|qu) = fa_scr_4index_(pt_ind,qu_ind)
              ! *********************************************************************

              call dgemm('t','n',nmo_p,nmo_p,df_vars_%nQ,d_val,                                       &
                   & fa_scr_(p_class)%p_sym(p_sym)%t_sym(t_sym)%val(:,pt_s+1:pt_s+nmo_p),df_vars_%nQ, &
                   & fa_scr_(p_class)%p_sym(p_sym)%t_sym(t_sym)%val(:,qu_s+1:qu_s+nmo_p),df_vars_%nQ, &
                   & 1.0_wp,fa_scr_4index_(1:nmo_p,1:nmo_p),nmo_p)

            end do ! t loop

          end do ! u loop

        end do ! t_sym loop

        p_max = last_index_(p_sym,p_class)

        qu_ind = 0

        do q = first_index_(p_sym,p_class) , last_index_(p_sym,p_class)

          qu_ind = qu_ind + 1

          q_i = trans_%class_to_irrep_map(q)

          pt_ind = qu_ind - 1

          do p = q , p_max

            pt_ind =pt_ind + 1

            p_i = trans_%class_to_irrep_map(p)

            fock_a_%occ(p_sym)%val(p_i,q_i) = fock_a_%occ(p_sym)%val(p_i,q_i) - &
                 & fa_scr_4index_(pt_ind,qu_ind)

          end do ! p loop

        end do ! q loop

      end do ! p_sym_loop

    end do ! p_class loop

    ! *******************************************************
    ! diagonal elements (p_class = q_class // enforce p > q )
    ! external orbitals (only need diagonal p = q elements) 
    ! *******************************************************

    do p_sym = 1 , nirrep_

      p_off = ndocpi_(p_sym) + nactpi_(p_sym)

      nmo_p = last_index_(p_sym,3) - first_index_(p_sym,3) + 1

      if ( nmo_p == 0 ) cycle

      do t_sym = 1 , nirrep_

        nmo_t = nactpi_(t_sym)

        if ( nmo_t == 0 ) cycle

        t_off = ndocpi_(t_sym)

        do u = first_index_(t_sym,2) , last_index_(t_sym,2)

          qu_s = ( trans_%class_to_irrep_map(u) - 1 - t_off ) * nmo_p

          do t = first_index_(t_sym,2) , last_index_(t_sym,2)

            pt_ind = ( trans_%class_to_irrep_map(t) - 1 - t_off ) * nmo_p

            tu_den = dens_%gemind(t,u)

            d_val = 0.5_wp * den1(tu_den)

            qu_ind = qu_s

            do p = first_index_(p_sym,3) , last_index_(p_sym,3)
  
              p_i = trans_%class_to_irrep_map(p)

              pt_ind = pt_ind + 1
  
              qu_ind = qu_ind + 1

              fock_a_%ext(p_sym)%val(p_i-p_off) = fock_a_%ext(p_sym)%val(p_i-p_off) -                &
                 & d_val * my_ddot(df_vars_%nQ,fa_scr_(3)%p_sym(p_sym)%t_sym(t_sym)%val(:,pt_ind),1, & 
                 &                 fa_scr_(3)%p_sym(p_sym)%t_sym(t_sym)%val(:,qu_ind),1)

            end do ! p loop

          end do ! t loop

        end do ! u loop

      end do ! t_sym loop

    end do ! p_sym_loop

    ! deallocate temporary scratch arrays
    call deallocate_fock_scr(fa_scr_,fa_scr_4index_)

    return

  end subroutine compute_f_a_df_exchange_fast

  subroutine compute_f_a_df_coulomb(den1,int2)

    implicit none

    real(wp), intent(in) :: den1(:), int2(:)

    integer :: p_class,q_class,q_class_max
    integer :: tu_den,qt_int
    integer :: t_sym,u_sym,p_sym,tu_sym,q_sym,qt_sym
    integer :: t,u,p,q,q_max,q_min
    integer :: p_i,q_i
    integer(ip) :: tu_df,pq_df,qu_df,pt_df,nQ
    integer :: tdf,udf,pdf,qdf

    real(wp) :: val

    nQ = df_vars_%nQ

    ! initialize
    do p_sym = 1 , nirrep_

      if ( allocated( fock_a_%occ(p_sym)%val ) ) fock_a_%occ(p_sym)%val = 0.0_wp
      if ( allocated( fock_a_%ext(p_sym)%val ) ) fock_a_%ext(p_sym)%val = 0.0_wp

    end do

    ! *** Coulomb terms ***
    
    ! initialize
    qint_%tuQ(1)%val(:,1) = 0.0_wp 

    ! gather intermediate
 
    do t_sym = 1 , nirrep_

      do t = first_index_(t_sym,2) , last_index_(t_sym,2)

        tdf = df_vars_%class_to_df_map(t)

        do u = first_index_(t_sym,2) , t - 1

          udf = df_vars_%class_to_df_map(u)

          tu_den = dens_%gemind(t,u)

          tu_df  = df_pq_index(tdf,udf) 

          call my_daxpy(df_vars_%nQ,2.0_wp*den1(tu_den),int2(tu_df+1:),df_vars_%Qstride,qint_%tuQ(1)%val(:,1),1)        

        end do ! end u loop

        tu_den = dens_%gemind(t,t)

        tu_df = df_pq_index(tdf,tdf) 

        call my_daxpy(df_vars_%nQ,den1(tu_den),int2(tu_df+1:),df_vars_%Qstride,qint_%tuQ(1)%val(:,1),1)

      end do ! end t loop

    end do ! end t_sym loop

    ! calculate fock matrix element

    do p_class = 1 , 3

      do q_class = 1 , 3

        if ( q_class > p_class ) cycle

        do p_sym = 1 , nirrep_

          do p = first_index_(p_sym,p_class) , last_index_(p_sym,p_class)
 
            pdf = df_vars_%class_to_df_map(p)

            q_max = last_index_(p_sym,q_class)

            q_min = first_index_(p_sym,q_class)

            if ( p_class == q_class ) q_max = p
   
            if ( q_class == 3 ) q_min = p

            do q = q_min , q_max

              qdf = df_vars_%class_to_df_map(q)

              pq_df = df_pq_index(pdf,qdf) 

              val = my_ddot(df_vars_%nQ,int2(pq_df+1:),df_vars_%Qstride,qint_%tuQ(1)%val(:,1),1)

              if ( ( p_class == 3 ) .and. ( q_class == 3 ) ) then

                p_i = trans_%class_to_irrep_map(p)-ndocpi_(p_sym)-nactpi_(p_sym)
                fock_a_%ext(p_sym)%val(p_i) = val

              else

                p_i = trans_%class_to_irrep_map(p)
                q_i = trans_%class_to_irrep_map(q)
                fock_a_%occ(p_sym)%val(p_i,q_i) = val

              endif
              
            end do ! end q loop

          end do ! end p loop
          
        end do ! end p_sym loop     

      end do ! end q_class loop

    end do ! end p_class loop

  end subroutine compute_f_a_df_coulomb

  subroutine compute_f_a_df_exchange(den1,int2)
    implicit none
    ! this subroutine computes the "active fock matrix" according to Eq. 12.5.13 in Helgaker on page 622
    ! F_a(m|n) = h(m|n) + SUM[v,w \in A] { g(mn|vw) - 0.5 g(mw|vn) }
    ! Apart from the symmetry constraint m_sym == n_sym and m >= n, m and n can belong to either of the three
    ! orbital classes.
    ! v and w belong to the active orbital class
    real(wp), intent(in) :: den1(:),int2(:)
    integer :: m_class,n_class,m_sym,m,n,v,w,w_sym,mdf,ndf,wdf,vdf
    integer :: den_ind,n_first,m_i,n_i
    integer(ip) :: mw,vn,wn,ww,nQ 
    real(wp) :: val,dval
    real(wp) :: v_mw(df_vars_%nQ)

    nQ = int(df_vars_%nQ,kind=ip)

    ! loop over orbital classes for m

    do m_class = 1 , 3

      ! loop over irreps (m_sym == n_sym for F(m,n) != 0)

      do m_sym = 1 , nirrep_

        ! loop over m indeces

        do m = first_index_(m_sym,m_class) , last_index_(m_sym,m_class)

          ! orbital index in df order
          mdf = df_vars_%class_to_df_map(m)

          ! loop over n indeces ( m_class == n_class && m >= n && m_sym == n_sym )

          n_first = first_index_(m_sym,m_class)

          ! only compute diagonal elements for external-external block

          if ( m_class == 3 ) n_first = m

!$omp parallel shared(ints_,nQ,dens_,m_class,m_sym,m,mdf,first_index_,  &
!$omp last_index_,df_vars_,den1,int2,n_first,fock_a_,trans_) num_threads(nthread_use_)
!$omp do private(n,ndf,val,v_mw,w_sym,w,wdf,v,den_ind,dval,vdf, &
!$omp ww,mw,wn,vn,n_i,m_i)

          do n = n_first , m

            ! orbital index in df order
            ndf = df_vars_%class_to_df_map(n)

            ! initialize Fock matrix element
            val = 0.0_wp

            ! loop over irreps for v

            do w_sym = 1 , nirrep_

              ! loop over w indeces

              do w = first_index_(w_sym,2) , last_index_(w_sym,2)

                ! orbital index in df order
                wdf = df_vars_%class_to_df_map(w)

                mw      = df_pq_index(mdf,wdf) 

                call my_dcopy(df_vars_%nQ,int2(mw+1:),df_vars_%Qstride,v_mw,1)

                ! loop over v indeces

                do v = first_index_(w_sym,2) , last_index_(w_sym,2) !w - 1 

                  ! 1-e density element
                  den_ind = dens_%gemind(v,w)
                  dval    = den1(den_ind)

                  ! orbital index in df order
                  vdf     = df_vars_%class_to_df_map(v)

                  ! 2-e exchange contribution - g(mw|vn)
                  vn      = df_pq_index(vdf,ndf) 

                  ! contract with integral/density matrix elements
                  val     = val - dval * my_ddot(df_vars_%nQ,v_mw,1,int2(vn+1:),df_vars_%Qstride)

                end do ! end v loop

              end do ! end w loop

            end do ! end w_sym loop

            ! save Fock matrix element

            if ( m_class == 3 ) then

              m_i = trans_%class_to_irrep_map(m)-ndocpi_(m_sym)-nactpi_(m_sym)

              fock_a_%ext(m_sym)%val(m_i) = fock_a_%ext(m_sym)%val(m_i) + 0.5_wp * val

            else

              m_i = trans_%class_to_irrep_map(m)
              n_i = trans_%class_to_irrep_map(n)

              fock_a_%occ(m_sym)%val(m_i,n_i) = fock_a_%occ(m_sym)%val(m_i,n_i) + 0.5_wp * val

            endif

          end do ! end n loop

!$omp end do nowait
!$omp end parallel

          ! loop over orbital classes for n (n_class < m_class --> n < m) 

          do n_class = 1 , m_class - 1

            ! loop over n indeces

!$omp parallel shared(n_class,nQ,ints_,dens_,m_class,m_sym,m,mdf,first_index_, &
!$omp last_index_,df_vars_,den1,int2,fock_a_,trans_) num_threads(nthread_use_)
!$omp do private(n,ndf,val,v_mw,w_sym,w,wdf,v,den_ind,dval,vdf,&
!$omp mw,wn,n_i,m_i)

            do n = first_index_(m_sym,n_class) , last_index_(m_sym,n_class)

              ! orbital index in df order
              ndf = df_vars_%class_to_df_map(n)

              ! initialize Fock matrix element
              val = 0.0_wp

              ! loop over irreps for v

              do w_sym = 1 , nirrep_

                ! loop over w indeces

                do w = first_index_(w_sym,2) , last_index_(w_sym,2)

                  ! orbital index in df order
                  wdf = df_vars_%class_to_df_map(w)

                  mw      = df_pq_index(mdf,wdf)              

                  call my_dcopy(df_vars_%nQ,int2(mw+1:),df_vars_%Qstride,v_mw,1)

                  ! loop over v indeces

                  do v = first_index_(w_sym,2) , last_index_(w_sym,2) !w -1 

                    ! 1-e density element
                    den_ind = dens_%gemind(v,w)
                    dval    = den1(den_ind)

                    ! orbital index in df order
                    vdf     = df_vars_%class_to_df_map(v)

                    ! 2-e exchange contribution - g(mw|vn)
                    vn      = df_pq_index(vdf,ndf) 

                    ! contract with integral/density matrix elements
                    val     = val - dval * my_ddot(df_vars_%nQ,v_mw,1,int2(vn+1:),df_vars_%Qstride)

                  end do ! end v loop

                end do ! end w loop

              end do ! end w_sym loop

              ! save Fock matrix element

              m_i = trans_%class_to_irrep_map(m)
              n_i = trans_%class_to_irrep_map(n)

              fock_a_%occ(m_sym)%val(m_i,n_i) = fock_a_%occ(m_sym)%val(m_i,n_i) + 0.5_wp * val

            end do ! end n loop

!$omp end do nowait
!$omp end parallel

          end do ! and n_class loop

        end do ! end m loop

      end do ! end m_sym loop

    end do ! end m_class loop

    return
  end subroutine compute_f_a_df_exchange

  subroutine compute_f_a(den1,int2)
    implicit none
    ! this subroutine computes the "active fock matrix" according to Eq. 12.5.13 in Helgaker on page 622
    ! F_a(m|n) = h(m|n) + SUM[v,w \in A] { g(mn|vw) - 0.5 g(mw|vn) }
    ! Apart from the symmetry constraint m_sym == n_sym and m >= n, m and n can belong to either of the three
    ! orbital classes.
    ! v and w belong to the active orbital class
    real(wp), intent(in) :: den1(:),int2(:)
    integer :: m_class,n_class,m_sym,m,n,v,w,w_sym,mw_sym
    integer :: mn,vw,mw,vn,den_ind,n_first,m_i,n_i
    integer :: int_ind,sym_offset
    real(wp) :: val,ival,dval

    ! initialize
    do m_sym = 1 , nirrep_

      if ( allocated( fock_a_%occ(m_sym)%val ) ) fock_a_%occ(m_sym)%val = 0.0_wp
      if ( allocated( fock_a_%ext(m_sym)%val ) ) fock_a_%ext(m_sym)%val = 0.0_wp

    end do

    ! loop over orbital classes for m

    do m_class = 1 , 3

      ! loop over irreps (m_sym == n_sym for F(m,n) != 0)

      do m_sym = 1 , nirrep_

        ! loop over m indeces

        do m = first_index_(m_sym,m_class) , last_index_(m_sym,m_class)

          ! loop over n indeces ( m_class == n_class && m >= n && m_sym == n_sym )

          n_first = first_index_(m_sym,m_class)
        
          if ( m_class == 3 ) n_first = m

          do n = n_first , m

            ! mn-geminal index
            mn  = ints_%gemind(m,n)
 
            ! initialize Fock matrix element
            val = 0.0_wp
 
            ! loop over irreps for v

            do w_sym = 1 , nirrep_
  
              mw_sym     = group_mult_tab_(m_sym,w_sym)
              sym_offset = ints_%offset(mw_sym)


              ! loop over w indeces

              do w = first_index_(w_sym,2) , last_index_(w_sym,2)

                ! loop over v indeces

                do v = first_index_(w_sym,2) , last_index_(w_sym,2)

                  ! 1-e density element
                  den_ind    = dens_%gemind(v,w)
                  dval       = den1(den_ind)

                  ! 2-e coulomb contribution 2 g(mn|vw)
                  vw         = ints_%gemind(v,w)
                  int_ind    = pq_index(mn,vw)
                  ival       = int2(int_ind)

                  ! 2-e exchange contribution - g(mw|vn)
                  mw         = ints_%gemind(m,w)
                  vn         = ints_%gemind(v,n)
                  int_ind    = sym_offset + pq_index(mw,vn)
                  ival       = ival - 0.5_wp * int2(int_ind)

                  ! contract with integral/density matrix elements
                  val        = val + dval * ival

                end do ! end v loop

              end do ! end w loop

            end do ! end w_sym loop

            ! save Fock matrix element

            if ( m_class == 3 ) then

              m_i                         = trans_%class_to_irrep_map(m)-&
                                          & ndocpi_(m_sym)-nactpi_(m_sym)
              fock_a_%ext(m_sym)%val(m_i) = val

            else

              m_i                             = trans_%class_to_irrep_map(m)
              n_i                             = trans_%class_to_irrep_map(n)
              fock_a_%occ(m_sym)%val(m_i,n_i) = val

            endif

          end do ! end n loop

          ! loop over orbital classes for n (n_class < m_class --> n < m) 

          do n_class = 1 , m_class - 1

            ! loop over n indeces

            do n = first_index_(m_sym,n_class) , last_index_(m_sym,n_class)

              ! mn-geminal index
              mn  = ints_%gemind(m,n)

              ! initialize Fock matrix element
              val = 0.0_wp

              ! loop over irreps for v

              do w_sym = 1 , nirrep_
  
                mw_sym     = group_mult_tab_(m_sym,w_sym)
                sym_offset = ints_%offset(mw_sym)

                ! loop over w indeces

                do w = first_index_(w_sym,2) , last_index_(w_sym,2)

                  ! loop over v indeces

                  do v = first_index_(w_sym,2) , last_index_(w_sym,2)

                    ! 1-e density element
                    den_ind    = dens_%gemind(v,w)
                    dval       = den1(den_ind)

                    ! 2-e coulomb contribution 2 g(mn|vw)
                    vw         = ints_%gemind(v,w)
                    int_ind    = pq_index(mn,vw)
                    ival       = int2(int_ind)

                    ! 2-e exchange contribution - g(mw|vn)
                    mw         = ints_%gemind(m,w)
                    vn         = ints_%gemind(v,n)
                    int_ind    = sym_offset + pq_index(mw,vn)
                    ival       = ival - 0.5_wp * int2(int_ind)

                    ! contract with integral/density matrix elements
                    val        = val + dval * ival

                  end do ! end v loop

                end do ! end w loop
  
              end do ! end w_sym loop

              ! save Fock matrix element

              m_i = trans_%class_to_irrep_map(m)
              n_i = trans_%class_to_irrep_map(n)

              fock_a_%occ(m_sym)%val(m_i,n_i) = val

            end do ! end n loop

          end do ! and n_class loop

        end do ! end m loop

      end do ! end m_sym loop

    end do ! end m_class loop
   
    return
  end subroutine compute_f_a

  subroutine compute_f_i_df_coulomb(int1,int2)

    implicit none

    real(wp), intent(in) :: int1(:),int2(:)

    integer :: i_sym,p_sym
    integer :: i,p,q,p_i,q_i,q_min,q_max
    integer :: pq_int
    integer :: p_class,q_class
    integer :: idf,pdf,qdf
    integer(ip) :: ii_df,pq_df,nQ
    real(wp) :: val

    nQ = int(df_vars_%nQ,kind=ip)

    ! initialize
    do p_sym = 1 , nirrep_

      if ( allocated( fock_i_%occ(p_sym)%val ) ) fock_i_%occ(p_sym)%val = 0.0_wp
      if ( allocated( fock_i_%ext(p_sym)%val ) ) fock_i_%ext(p_sym)%val = 0.0_wp

    end do

    qint_%tuQ(1)%val(:,1) = 0.0_wp

    ! calculate intermediate

    do i_sym = 1 , nirrep_

      do i = first_index_(i_sym,1) , last_index_(i_sym,1)

        idf   = df_vars_%class_to_df_map(i)

        ii_df = df_pq_index(idf,idf) 

        call my_daxpy(df_vars_%nQ,2.0_wp,int2(ii_df+1:),df_vars_%Qstride,qint_%tuQ(1)%val(:,1),1)

      end do ! end i loop

    end do ! end i_sym loop

    ! calculate Fock matrix element

    do p_class = 1 , 3

      do q_class = 1 , 3

        if ( q_class > p_class ) cycle

        do p_sym = 1 , nirrep_

          do p = first_index_(p_sym,p_class) , last_index_(p_sym,p_class)

            pdf = df_vars_%class_to_df_map(p)

            q_max = last_index_(p_sym,q_class)

            q_min = first_index_(p_sym,q_class)

            if ( p_class == q_class ) q_max = p

            if ( q_class == 3 ) q_min = p

            do q = q_min , q_max

              qdf    = df_vars_%class_to_df_map(q)

              pq_df  = df_pq_index(pdf,qdf) 

              pq_int = ints_%gemind(p,q) 

              val = int1(pq_int) + & 
                  & my_ddot(df_vars_%nQ,int2(pq_df+1:),df_vars_%Qstride,qint_%tuQ(1)%val(:,1),1)

              if ( ( p_class == 3 ) .and. ( q_class == 3 ) ) then

                p_i = trans_%class_to_irrep_map(p)-ndocpi_(p_sym)-nactpi_(p_sym)
                fock_i_%ext(p_sym)%val(p_i) = val

              else

                p_i = trans_%class_to_irrep_map(p)
                q_i = trans_%class_to_irrep_map(q)
                fock_i_%occ(p_sym)%val(p_i,q_i) = val

              endif

            end do ! end q loop

          end do ! end p loop

        end do ! end p_sym loop     

      end do ! end q_class loop

    end do ! end p_class loop

    return

  end subroutine compute_f_i_df_coulomb

  subroutine compute_f_i_df_exchange(int2)
    ! this subroutine computes the "inactive fock matrix" according to Eq. 12.5.12 in Helgaker on page 622
    ! F_i(m|n) = h(m|n) + SUM[i \in D] { 2 g(mn|ii) - g(mi|in) }
    ! using 3-index 2-e integrals
    ! Apart from the symmetry constraint m_sym == n_sym and m >= n, m and n can belong to either of the three
    ! orbital classes. 
    implicit none
    real(wp), intent(in) :: int2(:)
    integer :: m_class,n_class,m_sym,m,n,i,i_sym,mdf,ndf,idf
    integer :: n_first,m_i,n_i
    integer(ip) :: mi,in,nQ
    real(wp) :: val,int_val

    nQ = int(df_vars_%nQ,kind=ip) 

    ! loop over orbital classes for m

    do m_class = 1 , 3

      ! loop over irreps (m_sym == n_sym for F(m,n) != 0)

      do m_sym = 1 , nirrep_

        ! loop over m indeces

        do m = first_index_(m_sym,m_class) , last_index_(m_sym,m_class)

          ! orbital index in df order
          mdf = df_vars_%class_to_df_map(m)

          ! loop over n indeces ( m_class == n_class && m >= n && m_sym == n_sym )

          n_first = first_index_(m_sym,m_class)

          ! only compute diagonal elements for external-external block

          if ( m_class == 3 ) n_first = m

!$omp parallel shared(m,mdf,first_index_,last_index_,fock_i_,int2,  &
!$omp ints_,df_vars_,m_class,m_sym,n_first,trans_) num_threads(nthread_use_)
!$omp do private(n,val,ndf,idf,mi,in,m_i,n_i,i,i_sym)

          do n = n_first , m

            ! orbital index in df order
            ndf = df_vars_%class_to_df_map(n)

            val = 0.0_wp

            ! loop over irreps for i

            do i_sym = 1 , nirrep_

              ! loop over i indeces

              do i = first_index_(i_sym,1) , last_index_(i_sym,1)

                ! orbital index in df order
                idf        = df_vars_%class_to_df_map(i) 

                ! geminal indeces in df order
                mi         = df_pq_index(mdf,idf) 
                in         = df_pq_index(idf,ndf) 

                ! 2-e exchange contribution - g(mi|in)

                val        = val -  my_ddot(df_vars_%nQ,int2(mi+1:),df_vars_%Qstride,int2(in+1:),df_vars_%Qstride)

              end do ! end i loop

            end do ! end i_sym loop

            ! save Fock matrix element

            if ( m_class == 3 ) then

              m_i = trans_%class_to_irrep_map(m)-ndocpi_(m_sym)-nactpi_(m_sym)
              fock_i_%ext(m_sym)%val(m_i) = fock_i_%ext(m_sym)%val(m_i) + val

            else

              m_i = trans_%class_to_irrep_map(m)
              n_i = trans_%class_to_irrep_map(n)
              fock_i_%occ(m_sym)%val(m_i,n_i) = fock_i_%occ(m_sym)%val(m_i,n_i) + val

            endif

          end do ! end n loop

!$omp end do nowait
!$omp end parallel

          ! loop over orbital classes for n (n_class < m_class --> n < m) 

          do n_class = 1 , m_class - 1

            ! loop over n indeces

!$omp parallel shared(n_class,m,mdf,first_index_,last_index_,fock_i_,int2, &
!$omp ints_,df_vars_,m_class,m_sym,trans_) num_threads(nthread_use_)
!$omp do private(n,val,ndf,idf,mi,in,m_i,n_i,i_sym,i)

            do n = first_index_(m_sym,n_class) , last_index_(m_sym,n_class)

              ! orbital index in df order
              ndf     = df_vars_%class_to_df_map(n)

              val = 0.0_wp

              ! loop over irreps for i

              do i_sym = 1 , nirrep_

                ! loop over i indeces

                do i = first_index_(i_sym,1) , last_index_(i_sym,1)

                  ! orbital index in df order
                  idf        = df_vars_%class_to_df_map(i)

                  ! geminal indeces in df order
                  mi         = df_pq_index(mdf,idf) 
                  in         = df_pq_index(idf,ndf) 
 
                  ! 2-e exchange contribution - g(mi|in)

                  val        = val - my_ddot(df_vars_%nQ,int2(mi+1:),df_vars_%Qstride,int2(in+1:),df_vars_%Qstride)

                end do ! end i loop

              end do ! end i_sym loop

              ! save Fock matrix element

              m_i              = trans_%class_to_irrep_map(m)
              n_i              = trans_%class_to_irrep_map(n)

              fock_i_%occ(m_sym)%val(m_i,n_i) = fock_i_%occ(m_sym)%val(m_i,n_i) + val

            end do ! end n_loop

!$omp end do nowait
!$omp end parallel

          end do ! end n_class loop

        end do ! end m_loop

      end do ! end m_sym loop

    end do ! end m_class loop

    return

  end subroutine compute_f_i_df_exchange

  subroutine compute_f_i_df_exchange_fast(int2)

    implicit none
    ! subroutine to compute exchange contributions to the active Fock matrix

    ! Fi(p,q) = Fi(p,q) - SUM_[i \in D] (pi|qi)

    ! The subroutine performs the comptations in 3 steps
    !    1) collect ALL three-index integrals and sort according to p_class,p_sym,t_sym
    !    2) evaluate required 4-index integrals using DGEMMs 
    !    3) update Fock matrix elements

    real(wp), intent(in) :: int2(:)

    integer :: p_class,q_class
    integer :: p_sym,i_sym
    integer :: p_max
    integer :: p,q,i
    integer :: p_i,q_i,i_i
    integer :: p_off
    integer :: pdf,qdf,idf
    integer :: nmo_p,nmo_i,nmo_q,ngem_pi,ngem_qi
    integer :: pi_s,qi_s,qi_ind,pi_ind
    integer(ip) :: pi,qi
    integer     :: max_gen
    real(wp), allocatable :: pi_int(:,:),qi_int(:,:),fi_scr_4index(:,:)

    ! *************************
    ! allocate temporary arrays
    ! *************************

    max_gen = max(maxval(ndocpi_),maxval(nactpi_),maxval(nextpi_))

    allocate(pi_int(df_vars_%nQ,max_gen))

    allocate(qi_int(df_vars_%nQ,max_gen))

    allocate(fi_scr_4index(max_gen,max_gen))

    ! ****************************************************
    ! off diagonal elements (p_class > q_class --> p > q )
    ! ****************************************************

    do p_class = 1 , 3

      do q_class = 1 , p_class - 1

        do p_sym = 1 , nirrep_

          nmo_p = last_index_(p_sym,p_class) - first_index_(p_sym,p_class) + 1

          if ( nmo_p == 0 ) cycle

          nmo_q = last_index_(p_sym,q_class) - first_index_(p_sym,q_class) + 1

          if ( nmo_q == 0 ) cycle

          fi_scr_4index = 0.0_wp

          do i_sym = 1 , nirrep_

            nmo_i    = ndocpi_(i_sym)

            if ( nmo_i == 0 ) cycle

            do i = first_index_(i_sym,1) , last_index_(i_sym,1)

              idf    = df_vars_%class_to_df_map(i)

              pi_ind = 0

              do p = first_index_(p_sym,p_class) , last_index_(p_sym,p_class)

                pdf = df_vars_%class_to_df_map(p)

                pi  = df_pq_index(pdf,idf)

                pi_ind = pi_ind + 1

                call my_dcopy(df_vars_%nQ,int2(pi+1:),df_vars_%Qstride, &
                   & pi_int(:,pi_ind),1)

              end do

              qi_ind = 0

              do q = first_index_(p_sym,q_class) , last_index_(p_sym,q_class)

                qdf = df_vars_%class_to_df_map(q)

                qi  = df_pq_index(qdf,idf)

                qi_ind = qi_ind + 1

                call my_dcopy(df_vars_%nQ,int2(qi+1:),df_vars_%Qstride, &
                   & qi_int(:,qi_ind),1)

              end do

              call dgemm('t','n',nmo_p,nmo_q,df_vars_%nQ,1.0_wp,  &
                     & pi_int(:,1:nmo_p),df_vars_%nQ,             &
                     & qi_int(:,1:nmo_q),df_vars_%nQ,             &
                     & 1.0_wp,fi_scr_4index(1:nmo_p,1:nmo_q),nmo_p)

            end do ! i loop

          end do ! i_sym loop
 
          qi_ind = 0

          do q = first_index_(p_sym,q_class) , last_index_(p_sym,q_class)

            q_i = trans_%class_to_irrep_map(q)

            qi_ind = qi_ind + 1

            pi_ind = 0

            do p = first_index_(p_sym,p_class) , last_index_(p_sym,p_class)

              p_i = trans_%class_to_irrep_map(p)

              pi_ind = pi_ind + 1

              fock_i_%occ(p_sym)%val(p_i,q_i) = fock_i_%occ(p_sym)%val(p_i,q_i) - &
                    & fi_scr_4index(pi_ind,qi_ind)

            end do ! p loop

          end do ! q loop

        end do ! p_sym

      end do ! q_class

    end do ! p_class

    ! *******************************************************
    ! diagonal elements (p_class = q_class // enforce p > q )
    ! doubly-occupied and active orbitals (need all elements)
    ! *******************************************************

    do p_class = 1 , 2

      do p_sym = 1 , nirrep_

        p_off = 0

        if ( p_class == 2 ) p_off = ndocpi_(p_sym)

        nmo_p = last_index_(p_sym,p_class) - first_index_(p_sym,p_class) + 1

        if ( nmo_p == 0 ) cycle

        fi_scr_4index = 0.0_wp

        do i_sym = 1 , nirrep_

          nmo_i = ndocpi_(i_sym)

          if ( nmo_i == 0 ) cycle

          do i = first_index_(i_sym,1) , last_index_(i_sym,1)

            idf    = df_vars_%class_to_df_map(i)

            pi_ind = 0

            do p = first_index_(p_sym,p_class) , last_index_(p_sym,p_class)

              pdf = df_vars_%class_to_df_map(p)

              pi  = df_pq_index(pdf,idf)

              pi_ind = pi_ind + 1

              call my_dcopy(df_vars_%nQ,int2(pi+1:),df_vars_%Qstride, &
                 & pi_int(:,pi_ind),1)

            end do

            call dgemm('t','n',nmo_p,nmo_p,df_vars_%nQ,1.0_wp, &
                   & pi_int(:,1:nmo_p),df_vars_%nQ,            &
                   & pi_int(:,1:nmo_p),df_vars_%nQ,            &
                   & 1.0_wp,fi_scr_4index(1:nmo_p,1:nmo_p),nmo_p)

          end do ! i loop

        end do ! i_sym loop

        p_max = last_index_(p_sym,p_class)

        qi_ind = 0

        do q = first_index_(p_sym,p_class) , last_index_(p_sym,p_class)

          q_i = trans_%class_to_irrep_map(q)

!          if ( p_class == 3 ) p_max = q

          qi_ind = qi_ind + 1

          pi_ind = qi_ind - 1

          do p = q , p_max

            p_i = trans_%class_to_irrep_map(p)

            pi_ind = pi_ind + 1

!            if ( p_class == 3 ) then
!
!              fock_i_%ext(p_sym)%val(p_i - p_off) = fock_i_%ext(p_sym)%val(p_i - p_off) - &
!                   & fi_scr_4index(pi_ind,qi_ind)
!
!            else

              fock_i_%occ(p_sym)%val(p_i,q_i) = fock_i_%occ(p_sym)%val(p_i,q_i) - &
                   & fi_scr_4index(pi_ind,qi_ind)

!            end if

          end do ! p loop

        end do ! q loop

      end do ! _sym_loop

    end do ! p_class loop

    ! *******************************************************
    ! diagonal elements (p_class = q_class // enforce p = q )
    ! external orbitals (only need diagonal p = q elements)
    ! *******************************************************

    do p_sym = 1 , nirrep_

      p_off = ndocpi_(p_sym) + nactpi_(p_sym)

      nmo_p = last_index_(p_sym,3) - first_index_(p_sym,3) + 1

      if ( nmo_p == 0 ) cycle

      do i_sym = 1 , nirrep_

        nmo_i = ndocpi_(i_sym)

        if ( nmo_i == 0 ) cycle

        do i = first_index_(i_sym,1) , last_index_(i_sym,1)

          idf    = df_vars_%class_to_df_map(i)

          do p = first_index_(p_sym,3) , last_index_(p_sym,3)

            pdf = df_vars_%class_to_df_map(p)

            pi  = df_pq_index(pdf,idf)

            p_i = trans_%class_to_irrep_map(p)

            fock_i_%ext(p_sym)%val(p_i - p_off) = fock_i_%ext(p_sym)%val(p_i - p_off) -            &
                   & my_ddot(df_vars_%nQ,int2(pi+1:),df_vars_%Qstride,int2(pi+1:),df_vars_%Qstride)

          end do ! p loop

        end do ! i loop

      end do ! i_sym loop

    end do ! p_sym loop

    ! ***************************
    ! deallocate temporary arrays
    ! ***************************

    deallocate(pi_int,qi_int,fi_scr_4index)

    return

  end subroutine compute_f_i_df_exchange_fast

  subroutine compute_f_i(int1,int2)
    ! this subroutine computes the "inactive fock matrix" according to Eq. 12.5.12 in Helgaker on page 622
    ! F_i(m|n) = h(m|n) + SUM[i \in D] { 2 g(mn|ii) - g(mi|in) }
    ! Apart from the symmetry constraint m_sym == n_sym and m >= n, m and n can belong to either of the three
    ! orbital classes. 
    implicit none
    real(wp), intent(in) :: int1(:),int2(:)
    integer :: m_class,n_class,m_sym,m,n,i,i_sym,mi_sym
    integer :: mn,ii,mi,in,m_i,n_i,n_first
    integer :: int_ind,sym_offset
    real(wp) :: val

    ! initialize
    do m_sym = 1 , nirrep_

      if ( allocated( fock_i_%occ(m_sym)%val ) )     fock_i_%occ(m_sym)%val = 0.0_wp
      if ( allocated( fock_i_%ext(m_sym)%val ) ) fock_i_%ext(m_sym)%val = 0.0_wp

    end do

    ! loop over orbital classes for m

    do m_class = 1 , 3

      ! loop over irreps (m_sym == n_sym for F(m,n) != 0)

      do m_sym = 1 , nirrep_

        ! loop over m indeces

        do m = first_index_(m_sym,m_class) , last_index_(m_sym,m_class)


          ! loop over n indeces ( m_class == n_class && m >= n && m_sym == n_sym )

          n_first = first_index_(m_sym,m_class)

          ! only compute diagonal elements for external-external block

          if ( m_class == 3 ) n_first = m

          do n = n_first , m

            ! mn-geminal index
            mn = ints_%gemind(m,n)
            ! 1-e contribution h(m|n)
            val = int1(mn)

            ! loop over irreps for i

            do i_sym = 1 , nirrep_

              mi_sym     = group_mult_tab_(m_sym,i_sym)
              sym_offset = ints_%offset(mi_sym)

              ! loop over i indeces

              do i = first_index_(i_sym,1) , last_index_(i_sym,1)

                ! 2-e coulomb contribution 2 g(mn|ii)
                ii         = ints_%gemind(i,i)
                int_ind    = pq_index(ii,mn)
                val        = val + 2.0_wp * int2(int_ind)

                ! 2-e exchange contribution - g(mi|in)
                mi         = ints_%gemind(m,i)
                in         = ints_%gemind(i,n)
                int_ind    = sym_offset + pq_index(mi,in)
                val        = val - int2(int_ind) 

              end do ! end i loop

            end do ! end i_sym loop

            ! save Fock matrix element

            if ( m_class == 3 ) then

              m_i                             = trans_%class_to_irrep_map(m)-ndocpi_(m_sym)-nactpi_(m_sym)
              fock_i_%ext(m_sym)%val(m_i)     = val

            else

              m_i                             = trans_%class_to_irrep_map(m)
              n_i                             = trans_%class_to_irrep_map(n)
              fock_i_%occ(m_sym)%val(m_i,n_i) = val

            endif

          end do ! end n loop

          ! loop over orbital classes for n (n_class < m_class --> n < m) 

          do n_class = 1 , m_class - 1

            ! loop over n indeces

            do n = first_index_(m_sym,n_class) , last_index_(m_sym,n_class)

              mn = ints_%gemind(m,n)
              ! 1-e contribution h(m|n)
              val = int1(mn) 

              ! loop over irreps for i

              do i_sym = 1 , nirrep_

                mi_sym     = group_mult_tab_(m_sym,i_sym)
                sym_offset = ints_%offset(mi_sym)

                ! loop over i indeces

                do i = first_index_(i_sym,1) , last_index_(i_sym,1)

                  ! 2-e coulomb contribution 2 g(mn|ii)
                  ii         = ints_%gemind(i,i)
                  int_ind    = pq_index(ii,mn)
                  val        = val + 2.0_wp * int2(int_ind)

                  ! 2-e exchange contribution - g(mi|in)
                  mi         = ints_%gemind(m,i)
                  in         = ints_%gemind(i,n)
                  int_ind    = sym_offset + pq_index(mi,in)
                  val        = val - int2(int_ind)

                end do ! end i loop
 
              end do ! end i_sym loop

              ! save Fock matrix element

              m_i              = trans_%class_to_irrep_map(m)
              n_i              = trans_%class_to_irrep_map(n)

              fock_i_%occ(m_sym)%val(m_i,n_i) = val

            end do ! end n_loop

          end do ! end n_class loop

        end do ! end m_loop

      end do ! end m_sym loop
  
    end do ! end m_class loop

    return
  end subroutine compute_f_i

  subroutine print_f_matrix(mat)

    implicit none

    type(fock_info) :: mat

    integer :: m_sym,n_class,m_class,m,n,m_i,n_i
    integer :: iprint

    do m_sym = 1 , nirrep_

      if ( ndocpi_(m_sym) + nactpi_(m_sym) > 0 ) then !allocated ( mat%occ(m_sym)%val ) ) then

        write(fid_,'(a,1x,i1,2(i3,1x))')'lower-triangular elements for irrep',m_sym,&
              & ndocpi_(m_sym) + nactpi_(m_sym) + nextpi_(m_sym) , ndocpi_(m_sym) + nactpi_(m_sym)

        iprint=0
 
        do m_class = 1 , 3 
 
          do n_class = 1 , 2

            do m = first_index_(m_sym,m_class) , last_index_(m_sym,m_class) 

              do n = first_index_(m_sym,n_class) , min( m,last_index_(m_sym,n_class) )

                iprint = iprint + 1

                m_i              = trans_%class_to_irrep_map(m)
                n_i              = trans_%class_to_irrep_map(n)
 
                write(fid_,'(2(i3,1x),es20.13,1x)',advance='no')m_i,n_i,mat%occ(m_sym)%val(m_i,n_i)

                if ( iprint < 4 ) cycle

                iprint = 0
                write(fid_,*)

              end do

            end do
  
          end do

        end do

        if ( iprint > 0 ) write(fid_,*)

      end if

      if ( nextpi_(m_sym) > 0 ) then

         write(fid_,'(a,1x,i1,2(i3,1x))')'diagonal elements for irrep',m_sym,nextpi_(m_sym)

         iprint = 0

         do m = first_index_(m_sym,3) , last_index_(m_sym,3)
           
            iprint = iprint + 1
 
            m_i = trans_%class_to_irrep_map(m)-ndocpi_(m_sym)-nactpi_(m_sym)

            write(fid_,'(2(i3,1x),es20.13,1x)',advance='no')m_i,n_i,mat%ext(m_sym)%val(m_i)

            if ( iprint < 4 ) cycle

            iprint = 0
            write(fid_,*)

         end do

         if ( iprint > 0 ) write(fid_,*)

      end if

    end do

    return

  end subroutine print_f_matrix

  subroutine print_q_matrix(mat)

    implicit none

    real(wp), intent(in) :: mat(:,:)
 
    integer :: m,n,m_i,n_i,m_class,m_sym,iprint

    do m_sym = 1 , nirrep_

      if ( nactpi_(m_sym) == 0 ) cycle

      write(fid_,'(a,1x,i1,2(i3,1x))')'matrix elements for irrep',m_sym,&
         & nactpi_(m_sym) , ndocpi_(m_sym) + nactpi_(m_sym) + nextpi_(m_sym)

      iprint = 0

      do m_class = 1 , 3

        do m = first_index_(m_sym,m_class) , last_index_(m_sym,m_class)

          do n = first_index_(m_sym,2) , last_index_(m_sym,2)

            iprint = iprint + 1

            m_i              = trans_%class_to_irrep_map(m)
            n_i              = trans_%class_to_irrep_map(n)

            write(fid_,'(2(i3,1x),es20.13,1x)',advance='no')m_i,n_i,mat(m-ndocpi_(m_sym),n)

            if ( iprint < 4 ) cycle

            iprint = 0
            write(fid_,*)

          end do

        end do

      end do

      if ( iprint > 0 ) write(fid_,*)

    end do

    return

  end subroutine print_q_matrix

  subroutine print_vector(vector)

    implicit none

    real(wp), intent(in) :: vector(:)

    integer :: a_sym,t_sym,i,a,t,u,grad_ind,newline
    character(1) :: a_typ,i_typ,t_typ,u_typ

    newline = 0

    ! ******************************
    ! doubly occupied - active pairs
    ! ******************************

    i_typ='d'
    t_typ='a'

    do t_sym = 1 , nirrep_

      grad_ind = rot_pair_%pair_offset(t_sym,rot_pair_%act_doc_type)

      do i = first_index_(t_sym,1) , last_index_(t_sym,1)

        do t = first_index_(t_sym,2) , last_index_(t_sym,2)

          grad_ind = grad_ind + 1

          write(fid_,'(a,a,a,a,i3,a,i3,a,1x,es10.3,4x)',advance='no')t_typ,'-',i_typ,&
               & ' (',t,',',i,')',vector(grad_ind)

          newline = newline + 1

          if ( mod(newline,4) == 0 ) write(fid_,*)

        end do

      end do

    end do

    ! ********************************
    ! doubly occupied - external pairs
    ! ********************************

    i_typ='d'
    a_typ='e'

    do a_sym = 1 , nirrep_

      grad_ind = rot_pair_%pair_offset(a_sym,rot_pair_%ext_doc_type)

      do i = first_index_(a_sym,1) , last_index_(a_sym,1)

        do a = first_index_(a_sym,3) , last_index_(a_sym,3)

          grad_ind = grad_ind + 1

          write(fid_,'(a,a,a,a,i3,a,i3,a,1x,es10.3,4x)',advance='no')a_typ,'-',i_typ,&
               & ' (',a,',',i,')',vector(grad_ind)

          newline = newline + 1

          if ( mod(newline,4) == 0 ) write(fid_,*)

        end do

      end do

    end do

    ! *********************
    ! active - active pairs
    ! *********************

    t_typ='a'
    u_typ='a'

    if ( include_aa_rot_ == 1 ) then

      do t_sym = 1 , nirrep_

        grad_ind = rot_pair_%pair_offset(t_sym,rot_pair_%act_act_type)

        do u = first_index_(t_sym,2) , last_index_(t_sym,2)

          do t = u + 1 , last_index_(t_sym,2)

            grad_ind = grad_ind + 1

            write(fid_,'(a,a,a,a,i3,a,i3,a,1x,es10.3,4x)',advance='no')t_typ,'-',u_typ,&
                 & ' (',t,',',u,')',vector(grad_ind)

            newline = newline + 1

            if ( mod(newline,4) == 0 ) write(fid_,*)

          end do

        end do

      end do

    end if

    ! ***********************
    ! active - external pairs
    ! ***********************

    t_typ='a'
    a_typ='e'

    do a_sym = 1 , nirrep_

      grad_ind = rot_pair_%pair_offset(a_sym,rot_pair_%ext_act_type)

      do t = first_index_(a_sym,2) , last_index_(a_sym,2)

        do a = first_index_(a_sym,3) , last_index_(a_sym,3)

          grad_ind = grad_ind + 1

          write(fid_,'(a,a,a,a,i3,a,i3,a,1x,es10.3,4x)',advance='no')a_typ,'-',t_typ,&
               & ' (',a,',',t,')',vector(grad_ind)

          newline = newline + 1

          if ( mod(newline,4) == 0 ) write(fid_,*)

        end do

      end do

    end do

    if ( mod(newline,4) /= 0 ) write(fid_,*)

    return

  end subroutine print_vector

  subroutine print_orbital_gradient()
    implicit none
    integer :: ij_pair,i,j,newline,i_class,j_class,j_class_start,j_start,i_sym
    character(1) :: ityp,jtyp
    ! loop over rotation pairs

    write(fid_,'(a)')'orbital gradient:'

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

              write(fid_,'(a,a,a,a,i3,a,i3,a,1x,es10.3,4x)',advance='no')jtyp,'-',ityp,' (',j,',',i,')',orbital_gradient_(ij_pair)

              newline = newline + 1

              if ( mod(newline,4) == 0 ) write(fid_,*)

            end do

          end do

        end do

      end do

    end do
    if ( mod(newline,4) /= 0 ) write(fid_,*)

    write(fid_,'(a,1x,es10.3)')'gradient norm:',grad_norm_

    return

  end subroutine print_orbital_gradient

  subroutine allocate_temporary_fock_matrices()
    implicit none
    integer :: i_sym,ntot,nocc
    ! the total number of nonzero LT elements in the active/inactive Fock matrix 
    ! is equal to the number of geminals in the totally symmetric irrep

    ! allocate Fock matrices for the external block, only the diagonal elements are needed)

    allocate(fock_i_%occ(nirrep_),fock_i_%ext(nirrep_))

    do i_sym = 1 , nirrep_

      nocc = nactpi_(i_sym)+ndocpi_(i_sym)
      ntot = nocc + nextpi_(i_sym)

      allocate(fock_i_%occ(i_sym)%val(ntot,nocc)) 

      allocate(fock_i_%ext(i_sym)%val(nextpi_(i_sym)))

    end do

    allocate(fock_a_%occ(nirrep_),fock_a_%ext(nirrep_))

    do i_sym = 1 , nirrep_

      nocc = nactpi_(i_sym)+ndocpi_(i_sym)
      ntot = nocc + nextpi_(i_sym)

      allocate(fock_a_%occ(i_sym)%val(ntot,nocc))

      allocate(fock_a_%ext(i_sym)%val(nextpi_(i_sym)))

    end do

    ! the row index here will be between 1-nact, while the column index is arbitrary

    allocate(q_(nact_tot_,nmo_tot_))

    allocate(z_(nact_tot_,nmo_tot_))

    return
  end subroutine allocate_temporary_fock_matrices

  subroutine deallocate_temporary_fock_matrices()
    implicit none

    integer :: i_sym

    if ( allocated(fock_i_%occ) ) then

      do i_sym = 1 , nirrep_

        if ( .not. allocated(fock_i_%occ(i_sym)%val) ) cycle

        deallocate(fock_i_%occ(i_sym)%val)

      end do

      deallocate(fock_i_%occ)

    endif

    if ( allocated(fock_i_%ext) ) then

      do i_sym = 1 , nirrep_

        if ( .not. allocated(fock_i_%ext(i_sym)%val) ) cycle

        deallocate(fock_i_%ext(i_sym)%val)

      end do

      deallocate(fock_i_%ext)

    end if

    if ( allocated(fock_a_%occ) ) then

      do i_sym = 1 , nirrep_

        if ( .not. allocated(fock_a_%occ(i_sym)%val) ) cycle

        deallocate(fock_a_%occ(i_sym)%val)

      end do

      deallocate(fock_a_%occ)

    endif

    if ( allocated(fock_a_%ext) ) then

      do i_sym = 1 , nirrep_

        if ( .not. allocated(fock_a_%ext(i_sym)%val) ) cycle

        deallocate(fock_a_%ext(i_sym)%val)

      end do

      deallocate(fock_a_%ext)

    end if

    if (allocated(q_))        deallocate(q_)
    if (allocated(z_))        deallocate(z_)
    return
  end subroutine deallocate_temporary_fock_matrices

  subroutine allocate_qint()
    implicit none

    integer :: i_sym

    allocate(qint_%tuQ(nirrep_))

    do i_sym = 1 , nirrep_

      allocate(qint_%tuQ(i_sym)%val(df_vars_%nQ,dens_%ngempi(i_sym)))
!      allocate(qint_%tuQ(i_sym)%val(df_vars_%nQ,df_vars_%noccgempi(i_sym)))

    end do

    return
  end subroutine allocate_qint

  subroutine deallocate_qint()
    implicit none

    integer :: i_sym

    if ( allocated(qint_%tuQ) ) then

      do i_sym = 1 , nirrep_

        if ( .not. allocated(qint_%tuQ(i_sym)%val) ) cycle

        deallocate(qint_%tuQ(i_sym)%val)

      end do

      deallocate(qint_%tuQ)

    end if

  end subroutine deallocate_qint

  subroutine allocate_fa_scr()

    implicit none

    integer :: p_class,p_sym,t_sym,nmo_p,nmo_t
    integer :: max_nact,max_gen

    do p_class = 1 , 3

      allocate(fa_scr_(p_class)%p_sym(nirrep_))

      do p_sym = 1 , nirrep_

        allocate(fa_scr_(p_class)%p_sym(p_sym)%t_sym(nirrep_))

        if ( p_class == 1 ) then

          nmo_p = ndocpi_(p_sym)

        elseif ( p_class == 2 ) then

          nmo_p = nactpi_(p_sym)

        else

          nmo_p = nextpi_(p_sym)

        end if

        if ( nmo_p == 0 ) cycle

        do t_sym = 1 , nirrep_

          nmo_t = nactpi_(t_sym)

          if ( nmo_t == 0 ) cycle

          allocate(fa_scr_(p_class)%p_sym(p_sym)%t_sym(t_sym)%val(df_vars_%nQ,nmo_p*nmo_t))

        end do

      end do

    end do

    max_nact = maxval(nactpi_)

    max_gen  = max(maxval(ndocpi_),max_nact,maxval(nextpi_))

    allocate(fa_scr_4index_(max_gen*max_nact,max_gen*max_nact))

    return

  end subroutine allocate_fa_scr

  subroutine deallocate_fa_scr()

    implicit none

    integer :: p_class,p_sym,t_sym

    do p_class = 1 , 3

      do p_sym = 1 , nirrep_

        do t_sym = 1 , nirrep_

          if ( .not. allocated(fa_scr_(p_class)%p_sym(p_sym)%t_sym(t_sym)%val ) ) cycle

          deallocate(fa_scr_(p_class)%p_sym(p_sym)%t_sym(t_sym)%val)

        end do

        deallocate(fa_scr_(p_class)%p_sym(p_sym)%t_sym)

      end do

      deallocate(fa_scr_(p_class)%p_sym)

    end do

    if ( allocated(fa_scr_4index_) ) deallocate(fa_scr_4index_)

    return

  end subroutine deallocate_fa_scr

  subroutine allocate_fock_scr(f_scr,f_scr_4index,nmopi_contract,contract_class)

    implicit none
 
    type(fock_scr_info)  :: f_scr(:)
    real(wp),allocatable :: f_scr_4index(:,:) 
    integer, intent(in)  :: nmopi_contract(:)
    integer, intent(in)  :: contract_class

    integer :: p_class,p_sym,t_sym,nmo_p,nmo_t
    integer :: max_nmo_contract,max_gen

    do p_class = 1 , 3

      allocate(f_scr(p_class)%p_sym(nirrep_))

      do p_sym = 1 , nirrep_

        allocate(f_scr(p_class)%p_sym(p_sym)%t_sym(nirrep_))

        if ( p_class == 1 ) then

          nmo_p = ndocpi_(p_sym)

        elseif ( p_class == 2 ) then

          nmo_p = nactpi_(p_sym)

        else

          nmo_p = nextpi_(p_sym)

        end if

        if ( nmo_p == 0 ) cycle

        do t_sym = 1 , nirrep_

          nmo_t = nmopi_contract(t_sym)

          if ( nmo_t == 0 ) cycle

          allocate(f_scr(p_class)%p_sym(p_sym)%t_sym(t_sym)%val(df_vars_%nQ,nmo_p*nmo_t))

        end do

      end do

    end do

    max_gen         = max(maxval(ndocpi_),maxval(nactpi_),maxval(nextpi_))

    allocate(f_scr_4index(max_gen,max_gen))

    return

  end subroutine allocate_fock_scr

  subroutine deallocate_fock_scr(f_scr,f_scr_4index)

    implicit none

    type(fock_scr_info)   :: f_scr(:) 
    real(wp), allocatable :: f_scr_4index(:,:) 

    integer :: p_class,p_sym,t_sym

    do p_class = 1 , 3

      do p_sym = 1 , nirrep_

        do t_sym = 1 , nirrep_

          if ( .not. allocated(f_scr(p_class)%p_sym(p_sym)%t_sym(t_sym)%val ) ) cycle

          deallocate(f_scr(p_class)%p_sym(p_sym)%t_sym(t_sym)%val)

        end do

        deallocate(f_scr(p_class)%p_sym(p_sym)%t_sym)

      end do

      deallocate(f_scr(p_class)%p_sym)

    end do

    if ( allocated(f_scr_4index) ) deallocate(f_scr_4index)

    return

  end subroutine deallocate_fock_scr

end module focas_gradient
