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

module focas_semicanonical

  use focas_data
  use focas_driver
  use focas_energy, only     : compute_energy
  use focas_transform_driver 
  use focas_redundant, only  : diagonalize_opdm_block

  implicit none

  contains

    subroutine compute_semicanonical_mos(den1,den2,int1,int2,den1_nnz,den2_nnz,int1_nnz,int2_nnz, &
                                       & mo_coeff,ndocpi,nactpi,nextpi,nirrep,orbopt_data,fname)

      implicit none

      integer, intent(in)     :: nirrep
      integer, intent(in)     :: den1_nnz
      integer, intent(in)     :: den2_nnz
      integer, intent(in)     :: int1_nnz
      integer(ip), intent(in) :: int2_nnz

      real(wp), intent(in)    :: den1(den1_nnz)
      real(wp), intent(in)    :: den2(den2_nnz)
      real(wp), intent(in)    :: int1(int1_nnz)
      real(wp), intent(in)    :: int2(int2_nnz)
 
      integer, intent(in)     :: ndocpi(nirrep)
      integer, intent(in)     :: nactpi(nirrep)
      integer, intent(in)     :: nextpi(nirrep)

      real(wp), intent(inout) :: orbopt_data(14)
      real(wp), intent(inout) :: mo_coeff(:,:)

      character(120)          :: fname
            
      logical :: fexist 
      integer :: error

      integer :: h,m,n,offset

      nthread_use_ = int(orbopt_data(1))
      log_print_   = int(orbopt_data(6)) 
      df_vars_%use_df_teints      = int(orbopt_data(10))

      if ( log_print_ == 1 ) then

        inquire(file=fname,exist=fexist)

        if (fexist) then

          open(fid_,file=fname,status='old',position='append')

        else

          open(fid_,file=fname,status='new')

        endif

        write(fid_,'(a)')'computing semicanonical orbitals'

      endif

      ! allocate variables and set up mapping arrays
      call allocate_semicanonical_initial(ndocpi,nactpi,nextpi,nirrep,int2_nnz)

      ! *** DEBUG STUFF

      call compute_energy(int1,int2,den1,den2)

      ! *** END DEBUG STUFF

      ! compute the fock matrices
      call compute_gen_fock(int1,int2,den1)

      ! compute semicanonical orbitals and orbital energies
      error=diagonalize_gen_fock()

      ! adjust phase of each transformation vecto so that the element
      ! with the largest magnitude if positive
      ! error=adjust_phase()

      ! copy MO coefficients into U
      error=copy_semicanonical_mos()

      ! mo_coeff matrix
      error = transform_mocoeff(mo_coeff)
      if ( error /= 0 ) call abort_print(50)

      ! *** DEBUG STUFF 

      ! 1-e integrals
      error = transform_oeints(int1)
      if ( error /= 0 ) call abort_print(30)

      ! 2-e integrals
      if ( df_vars_%use_df_teints == 0 ) then
        error = transform_teints(int2)
      else
        error = transform_teints_df(int2)
      end if
      if ( error /= 0 ) call abort_print(31)

      call compute_energy(int1,int2,den1,den2)
      
      ! *** END DEBUG STUFF

      ! deallocate variables
      call deallocate_semicanonical_final()

      return

    end subroutine compute_semicanonical_mos

    integer function adjust_phase()

      implicit none
 
      adjust_phase = adjust_phase_block(gen_f_%doc,ndocpi_)

      adjust_phase = adjust_phase_block(gen_f_%act,nactpi_)

      adjust_phase = adjust_phase_block(gen_f_%ext,nextpi_)

      adjust_phase = 0

      return

    end function adjust_phase

    integer function adjust_phase_block(f_block,dims)

      implicit none

      type(matrix_block)  :: f_block(nirrep_)
      integer, intent(in) :: dims(nirrep_)

      integer :: nmo,p_sym,p,q
      real(wp) :: max_val

      do p_sym = 1 , nirrep_

        nmo = dims(p_sym)

        if ( nmo == 0 ) cycle

        do p = 1 , nmo

          max_val = 0.0_wp

          do q = 1 , nmo

            if ( abs(max_val) > abs(f_block(p_sym)%val(q,p)) ) cycle

            max_val = f_block(p_sym)%val(q,p)

          end do

          if ( max_val > 0.0_wp ) cycle

          f_block(p_sym)%val(:,p) = - f_block(p_sym)%val(:,p)

        end do

      end do
        
      adjust_phase_block = 0 
 
      return

    end function adjust_phase_block

    subroutine compute_gen_fock(int1,int2,den1)

      implicit none

      real(wp), intent(in)  :: den1(:)
      real(wp), intent(in)  :: int1(:)
      real(wp), intent(in)  :: int2(:)

      real(wp), allocatable :: coulomb(:) 

      integer :: error
      integer :: offset(nirrep_)

      if ( df_vars_%use_df_teints == 1 ) then

        allocate(coulomb(df_vars_%nQ))

        ! precompute Coulob terms
        error = precompute_coulomb()
        if ( error /= 0 ) call abort_print(510)

        ! inactive-inactive block
        offset = 0
        if ( sum(ndocpi_) > 0 ) error=compute_gen_fock_block_df(1,gen_f_%doc,ndocpi_)
        if ( error /= 0 ) call abort_print(511)

        ! active-active block
        offset = ndocpi_
        if ( sum(nactpi_) > 0 ) error=compute_gen_fock_block_df(2,gen_f_%act,nactpi_)
        if ( error /= 0 ) call abort_print(512)
      
        ! external-external block
        offset = ndocpi_ + nactpi_
        if ( sum(nextpi_) > 0 ) error=compute_gen_fock_block_df(3,gen_f_%ext,nextpi_)
        if ( error /= 0 ) call abort_print(513)

        deallocate(coulomb)

      else

        ! inactive-inactive block
        offset = 0
        if ( sum(ndocpi_) > 0 ) error=compute_gen_fock_block(1,gen_f_%doc,ndocpi_)
        if ( error /= 0 ) call abort_print(521)

        ! active-active block
        offset = ndocpi_
        if ( sum(nactpi_) > 0 ) error=compute_gen_fock_block(2,gen_f_%act,nactpi_)
        if ( error /= 0 ) call abort_print(522)

        ! external-external block
        offset = ndocpi_ + nactpi_
        if ( sum(nextpi_) > 0 ) error=compute_gen_fock_block(3,gen_f_%ext,nextpi_)     
        if ( error /= 0 ) call abort_print(523)

      end if

      return

      contains

        integer function compute_gen_fock_block(p_class,f_block,dims)

          implicit none
   
          integer, intent(in) :: p_class
          integer, intent(in) :: dims(nirrep_)
 
          type(matrix_block)  :: f_block(nirrep_)

          integer  :: p_sym,i_sym,t_sym
          integer  :: p_i,q_i,off
          integer  :: p,q,i,t,u
          integer  :: ii,pq,pi,qi,tu,pt,qu
          integer  :: tu_den
          integer  :: pqii,piqi,pqtu,ptqu
          integer  :: int_sym_offset
          real(wp) :: f_val,f_tmp

          compute_gen_fock_block = 1 

          do p_sym = 1 , nirrep_

            if ( dims(p_sym) == 0 ) cycle

            if ( allocated(f_block(p_sym)%val) ) f_block(p_sym)%val = 0.0_wp

            off                = offset(p_sym)

            do p = first_index_(p_sym,p_class) , last_index_(p_sym,p_class)

              p_i = trans_%class_to_irrep_map(p) - off

              do q = p , last_index_(p_sym,p_class)

                q_i   = trans_%class_to_irrep_map(q) - off

                pq    = ints_%gemind(p,q)

                ! 1-e contribution

                f_val = int1(pq)

                ! inactive 2-e contribution

                do i_sym = 1 , nirrep_

                  int_sym_offset=ints_%offset(group_mult_tab_(p_sym,i_sym))

                  do i = first_index_(i_sym,1) , last_index_(i_sym,1)

                    ! coulomb term
                    ii    = ints_%gemind(i,i)
                    pqii  = pq_index(pq,ii)

                    f_val = f_val + 2.0_wp * int2(pqii)

                    ! exchange term
                    pi    = ints_%gemind(p,i)
                    qi    = ints_%gemind(q,i)
                    piqi  = pq_index(pi,qi) + int_sym_offset

                    f_val = f_val - int2(piqi)

                  end do
  
                end do

                ! active 2-e contribution
 
                do t_sym = 1 , nirrep_

                  int_sym_offset = ints_%offset(group_mult_tab_(p_sym,t_sym))

                  do t = first_index_(t_sym,2) , last_index_(t_sym,2)

                    pt = ints_%gemind(p,t)

                    do u = first_index_(t_sym,2) , last_index_(t_sym,2)

                      ! coulomb term
                      tu     = ints_%gemind(t,u)
                      pqtu   = pq_index(pq,tu)
 
                      f_tmp  = int2(pqtu) 

                      ! exchange term
                      qu     = ints_%gemind(q,u)
                      ptqu   = pq_index(pt,qu) + int_sym_offset

                      f_tmp  = f_tmp - 0.5_wp * int2(ptqu)

                      ! density value
                      tu_den = dens_%gemind(t,u)

                      ! update fock matrix element
                      f_val  = f_val + den1(tu_den) * f_tmp 

                    end do

                  end do

                end do

                f_block(p_sym)%val(q_i,p_i) = f_val
                f_block(p_sym)%val(p_i,q_i) = f_val

              end do

            end do

          end do

          compute_gen_fock_block = 0 

          return
          
        end function compute_gen_fock_block

        integer function precompute_coulomb()

          implicit none

          integer     :: i,t,u
          integer     :: i_sym,t_sym
          integer     :: idf,tdf,udf
          integer     :: tu_den
          integer(ip) :: ii_df,tu_df
          integer(ip) :: nQ

          precompute_coulomb = 1

          nQ = int ( df_vars_%nQ , kind = ip )

          ! precompute the Coulb terms

          coulomb = 0.0_wp

          do i_sym = 1 , nirrep_

            do i = first_index_(i_sym,1) , last_index_(i_sym,1)

              idf   = df_vars_%class_to_df_map(i)

              ii_df = df_pq_index(idf,idf) 

              call my_daxpy(df_vars_%nQ,2.0_wp,int2(ii_df+1:),df_vars_%Qstride,coulomb,1)

            end do

          end do

          do t_sym = 1 , nirrep_

            do t = first_index_(t_sym,2) , last_index_(t_sym,2)

              tdf = df_vars_%class_to_df_map(t)

              do u = first_index_(t_sym,2) , t - 1

                udf    = df_vars_%class_to_df_map(u)

                tu_df  = df_pq_index(tdf,udf) 

                tu_den = dens_%gemind(t,u)

                call my_daxpy(df_vars_%nQ,2.0_wp*den1(tu_den),int2(tu_df+1:),df_vars_%Qstride,coulomb,1)

              end do

              tu_df  = df_pq_index(tdf,tdf) 

              tu_den = dens_%gemind(t,t)

              call my_daxpy(df_vars_%nQ,den1(tu_den),int2(tu_df+1:),df_vars_%Qstride,coulomb,1)

            end do

          end do

          precompute_coulomb = 0

          return 

        end function precompute_coulomb

        integer function compute_gen_fock_block_df(p_class,f_block,dims)

          implicit none

          integer, intent(in)  :: p_class
          integer, intent(in)  :: dims(nirrep_)

          type(matrix_block)  :: f_block(nirrep_)

          integer      :: p_sym,i_sym,t_sym
          integer      :: p,q,i,t,u
          integer      :: pdf,qdf,idf,tdf,udf
          integer(ip)  :: pq_df,pi_df,qi_df,pt_df,qu_df
          integer(ip)  :: nQ
          integer      :: tu_den,pq_int
          integer      :: p_i,q_i
          integer      :: off
          real(wp)     :: f_val

          compute_gen_fock_block_df = 1

          nQ = int( df_vars_%nQ , kind = ip )

          do p_sym = 1 , nirrep_

            if ( dims(p_sym) == 0 ) cycle

            if ( allocated(f_block(p_sym)%val) ) f_block(p_sym)%val = 0.0_wp

            off                = offset(p_sym)

            do p = first_index_(p_sym,p_class) , last_index_(p_sym,p_class)

              p_i  = trans_%class_to_irrep_map(p) - off

              pdf  = df_vars_%class_to_df_map(p) 

              do q = p , last_index_(p_sym,p_class)

                q_i    = trans_%class_to_irrep_map(q) - off

                qdf    = df_vars_%class_to_df_map(q) 

                ! 1-e contribution
                pq_int = ints_%gemind(p,q)

                f_val  = int1(pq_int)

                ! 2-e coulomb contributions
                pq_df  = df_pq_index(pdf,qdf)

                f_val = f_val + my_ddot(df_vars_%nQ,int2(pq_df+1:pq_df+nQ),1,coulomb,1)  

                ! 2-e inactive exchange contribution

                do i_sym = 1 , nirrep_

                  do i = first_index_(i_sym,1) , last_index_(i_sym,1)
 
                    idf   = df_vars_%class_to_df_map(i)

                    pi_df = df_pq_index(pdf,idf) 
                    qi_df = df_pq_index(qdf,idf) 

                    f_val = f_val - my_ddot(df_vars_%nQ,int2(pi_df+1:),df_vars_%Qstride,int2(qi_df+1:),df_vars_%Qstride) 

                  end do

                end do

                ! 2-e active exchange contribution

                do t_sym = 1 , nirrep_
 
                  do t = first_index_(t_sym,2) , last_index_(t_sym,2)

                    tdf   = df_vars_%class_to_df_map(t)

                    pt_df = df_pq_index(pdf,tdf) 
 
                    do u = first_index_(t_sym,2) , last_index_(t_sym,2)

                      udf    = df_vars_%class_to_df_map(u)
                      
                      qu_df  = df_pq_index(qdf,udf) 

                      tu_den = dens_%gemind(t,u)

                      f_val = f_val - 0.5_wp*den1(tu_den)* &
                              & my_ddot(df_vars_%nQ,int2(pt_df+1:),df_vars_%Qstride,int2(qu_df+1:),df_vars_%Qstride)

                    end do

                  end do

                end do

                f_block(p_sym)%val(q_i,p_i) = f_val
                f_block(p_sym)%val(p_i,q_i) = f_val

              end do

            end do

          end do

          compute_gen_fock_block_df = 0

          return

        end function compute_gen_fock_block_df


    end subroutine compute_gen_fock

    integer function copy_semicanonical_mos()

      implicit none

      integer :: offset(nirrep_),dims(nirrep_)
      integer :: error,i_sym,i,j

!      real(wp), allocatable :: tmp(:,:)

      copy_semicanonical_mos = 1

      do i_sym = 1 , nirrep_

        trans_%u_irrep_block(i_sym)%val = 0.0_wp

      end do

      offset = 0
      dims   = ndocpi_
      if (sum(dims) > 0 ) error=copy_semicanonical_mos_block(gen_f_%doc)
      if ( error /= 0 ) call abort_print(531)

      offset = ndocpi_
      dims   = nactpi_
      if ( sum(dims) > 0 ) error=copy_semicanonical_mos_block(gen_f_%act)
      if ( error /= 0 ) call abort_print(532)

      offset = ndocpi_ + nactpi_
      dims   = nextpi_
      if ( sum(dims) > 0 ) error=copy_semicanonical_mos_block(gen_f_%ext)
      if ( error /= 0 ) call abort_print(533)

!      do i_sym = 1 , nirrep_
!
!        allocate(tmp(trans_%nmopi(i_sym),trans_%nmopi(i_sym)))
!        tmp = matmul(trans_%u_irrep_block(i_sym)%val,transpose(trans_%u_irrep_block(i_sym)%val))
!
!        write(*,*)'new symmetry'
!
!        do i=1,trans_%nmopi(i_sym)
!
!          do j=1,trans_%nmopi(i_sym)
!
!            write(*,*)i,j,trans_%u_irrep_block(i_sym)%val(i,j)
!
!          end do
!
!        end do
!
!        deallocate(tmp)
!
!      end do

      copy_semicanonical_mos = 0

      return

      contains

        integer function copy_semicanonical_mos_block(mo_block)

          implicit none

          type(matrix_block), intent(in) :: mo_block(nirrep_)

          integer :: i_sym,nmo,i_f,i_l

          copy_semicanonical_mos_block = 1 

          do i_sym = 1 , nirrep_

            nmo = dims(i_sym)

            if ( nmo == 0 ) cycle

            i_f = offset(i_sym) + 1
            i_l = i_f + dims(i_sym) - 1

            trans_%u_irrep_block(i_sym)%val(i_f:i_l,i_f:i_l) = mo_block(i_sym)%val

          end do

          copy_semicanonical_mos_block = 0
 
          return

        end function copy_semicanonical_mos_block

    end function copy_semicanonical_mos

    integer function diagonalize_gen_fock()

      ! the diagonalblocks of the Fock matrix are destroyed // replaced with eigenvectors

      implicit none

      integer :: error

      diagonalize_gen_fock = 1

      ! inactive block
      error = diagonalize_gen_fock_block(gen_f_%doc,gen_f_%doc_e,ndocpi_)
      if ( error /= 0 ) call abort_print(541)

      ! active block
      error = diagonalize_gen_fock_block(gen_f_%act,gen_f_%act_e,nactpi_)
      if ( error /= 0 ) call abort_print(542)

      ! active block
      error = diagonalize_gen_fock_block(gen_f_%ext,gen_f_%ext_e,nextpi_)
      if ( error /= 0 ) call abort_print(543)

      diagonalize_gen_fock = 0

      return

      contains

        integer function diagonalize_gen_fock_block(f_block,e_block,dims)

          implicit none

          type(matrix_block)  :: f_block(nirrep_)
          type(vector_block)  :: e_block(nirrep_)
          integer, intent(in) :: dims(nirrep_)

          integer :: p_sym,nmo

          do p_sym = 1 , nirrep_

            nmo = dims(p_sym)

            if ( nmo == 0 ) cycle

            error = diagonalize_opdm_block(e_block(p_sym)%val,f_block(p_sym)%val,nmo)
            if ( error /= 0 ) call abort_print(544)
 
          end do

          if ( log_print_ == 1 ) then
            error = print_orbital_energies(e_block)
            if (error /= 0) call abort_print(545)
          end if

          diagonalize_gen_fock_block = error

          return

        end function diagonalize_gen_fock_block

        integer function print_orbital_energies(e_block)

          implicit none

          type(vector_block), intent(in) :: e_block(nirrep_) 

          integer :: p_sym,nmo

          write(fid_,*)

          do p_sym = 1 , nirrep_

            if ( trans_%nmopi(p_sym) == 0 ) cycle

            if ( .not. allocated( e_block(p_sym)%val ) ) cycle
 
            nmo = size(e_block(p_sym)%val)

            write(fid_,'(a,1x,i4,1x,a,1x,i1)')'orbital energies for ',nmo,'orbitals with symmetry',p_sym

            write(fid_,'(8(f9.4,1x))')e_block(p_sym)%val

          end do

          print_orbital_energies = 0

          return

        end function print_orbital_energies
 
    end function diagonalize_gen_fock

    subroutine allocate_semicanonical_initial(ndocpi,nactpi,nextpi,nirrep,int2_nnz)

      implicit none

      integer, intent(in)     :: nirrep 
      integer(ip), intent(in) :: int2_nnz
      integer, intent(in)     :: ndocpi(nirrep)
      integer, intent(in)     :: nactpi(nirrep)
      integer, intent(in)     :: nextpi(nirrep)

      integer :: error
      integer :: nfzcpi(nirrep)

      ! calculate the total number of orbitals in space
      ndoc_tot_ = sum(ndocpi)
      nact_tot_ = sum(nactpi)
      next_tot_ = sum(nextpi)
      nmo_tot_  = ndoc_tot_+nact_tot_+next_tot_

      ! in this case, there are no frozen doubly occupied orbitals
      nfzc_tot_ = 0
      nfzcpi_   = 0

      ! allocate indexing arrays
      call allocate_indexing_arrays(nirrep)

      ! determine integral/density addressing arrays
      call setup_indexing_arrays(nfzcpi,ndocpi,nactpi,nextpi)

      ! allocate transformation matrices
      call allocate_transformation_matrices()

      ! determine indexing arrays (needed for sorts in the integral transformation step)
      call determine_transformation_maps()

      ! allocate blocks of generalized Fock matrix
      call allocate_generalized_fock_matrix()

      ! set up df mapping arrays if density-fitted 2-e integrals are used
      error = 0
      if ( df_vars_%use_df_teints == 1 ) error = df_map_setup(int2_nnz)
      if ( error /= 0 ) call abort_print(20)

      return

    end subroutine allocate_semicanonical_initial

    subroutine deallocate_semicanonical_final()

      implicit none 

      ! deallocate indexing arrays
      call deallocate_indexing_arrays()
       
      ! deallocate transformation matrices
      call deallocate_transformation_matrices()

      ! deallocate generalized fock matrices
      call deallocate_generalized_fock_matrix()
      
      if (allocated(df_vars_%class_to_df_map)) deallocate(df_vars_%class_to_df_map)
!      if (allocated(df_vars_%noccgempi))       deallocate(df_vars_%noccgempi)

      return

    end subroutine deallocate_semicanonical_final

    subroutine allocate_generalized_fock_matrix()

      implicit none

      integer :: i_sym,nmo
  
      allocate(gen_f_%doc(nirrep_))
      allocate(gen_f_%act(nirrep_))
      allocate(gen_f_%ext(nirrep_))
      allocate(gen_f_%doc_e(nirrep_))
      allocate(gen_f_%act_e(nirrep_))
      allocate(gen_f_%ext_e(nirrep_))
 
      ! inactive block
      do i_sym = 1 , nirrep_

        nmo = ndocpi_(i_sym)
 
        if ( nmo == 0 ) cycle

        allocate(gen_f_%doc(i_sym)%val(nmo,nmo))
       
        gen_f_%doc(i_sym)%val = 0.0_wp

      end do

      do i_sym = 1 , nirrep_

        nmo = ndocpi_(i_sym)

        if ( nmo == 0 ) cycle

        allocate(gen_f_%doc_e(i_sym)%val(nmo))

        gen_f_%doc_e(i_sym)%val = 0.0_wp

      end do

      ! active block
      do i_sym = 1 , nirrep_

        nmo = nactpi_(i_sym)

        if ( nmo == 0 ) cycle

        allocate(gen_f_%act(i_sym)%val(nmo,nmo))

        gen_f_%act(i_sym)%val = 0.0_wp

      end do

      do i_sym = 1 , nirrep_

        nmo = nactpi_(i_sym)

        if ( nmo == 0 ) cycle

        allocate(gen_f_%act_e(i_sym)%val(nmo))

        gen_f_%act_e(i_sym)%val = 0.0_wp

      end do

      ! external block
      do i_sym = 1 , nirrep_

        nmo = nextpi_(i_sym)

        if ( nmo == 0 ) cycle

        allocate(gen_f_%ext(i_sym)%val(nmo,nmo))

        gen_f_%ext(i_sym)%val = 0.0_wp

      end do

      do i_sym = 1 , nirrep_

        nmo = nextpi_(i_sym)

        if ( nmo == 0 ) cycle

        allocate(gen_f_%ext_e(i_sym)%val(nmo))

        gen_f_%ext_e(i_sym)%val = 0.0_wp

      end do

      return

    end subroutine allocate_generalized_fock_matrix

    subroutine deallocate_generalized_fock_matrix

      implicit none

      integer :: i_sym

      do i_sym = 1 , nirrep_

        if ( allocated(gen_f_%doc(i_sym)%val) )   deallocate(gen_f_%doc(i_sym)%val)
        if ( allocated(gen_f_%act(i_sym)%val) )   deallocate(gen_f_%act(i_sym)%val)
        if ( allocated(gen_f_%ext(i_sym)%val) )   deallocate(gen_f_%ext(i_sym)%val) 
        if ( allocated(gen_f_%doc_e(i_sym)%val) ) deallocate(gen_f_%doc_e(i_sym)%val)
        if ( allocated(gen_f_%act_e(i_sym)%val) ) deallocate(gen_f_%act_e(i_sym)%val)
        if ( allocated(gen_f_%ext_e(i_sym)%val) ) deallocate(gen_f_%ext_e(i_sym)%val)

      end do
 
      deallocate(gen_f_%doc,gen_f_%act,gen_f_%ext)
      deallocate(gen_f_%doc_e,gen_f_%act_e,gen_f_%ext_e)

      return
 
    end subroutine deallocate_generalized_fock_matrix

    subroutine check_gen_fock()

      implicit none

      integer :: nmo,p_sym
      integer :: p,q

      if ( log_print_ == 0 ) return
       
      write(fid_,'(a)')'printing nonzero elements of Fock matrix'

      do p_sym = 1 , nirrep_

        nmo = ndocpi_(p_sym)

        if ( nmo == 0 ) cycle

        do p = 1 , nmo

          do q = 1 , p

            if ( abs( gen_f_%doc(p_sym)%val(p,q) ) < 1.0e-15_wp ) cycle

            write(fid_,'(2(i3,1x),f9.4)')p,q,gen_f_%doc(p_sym)%val(p,q)

          end do
 
        end do

      end do

      do p_sym = 1 , nirrep_

        nmo = nactpi_(p_sym)

        if ( nmo == 0 ) cycle

        do p = 1 , nmo

          do q = 1 , p

            if ( abs( gen_f_%act(p_sym)%val(p,q) ) < 1.0e-15_wp ) cycle

            write(fid_,'(2(i3,1x),f9.4)')p,q,gen_f_%act(p_sym)%val(p,q)

          end do

        end do

      end do

      do p_sym = 1 , nirrep_

        nmo = nextpi_(p_sym)

        if ( nmo == 0 ) cycle

        do p = 1 , nmo

          do q = 1 , p

            if ( abs( gen_f_%ext(p_sym)%val(p,q) ) < 1.0e-15_wp ) cycle

            write(fid_,'(2(i3,1x),f9.4)')p,q,gen_f_%ext(p_sym)%val(p,q)

          end do

        end do

      end do

    end subroutine check_gen_fock

end module focas_semicanonical
