module focas_gradient_hessian

  use focas_data
  use focas_driver
  use focas_gradient
  use focas_hessian
  use focas_exponential

  integer, allocatable :: rot_pair_map(:,:)

  contains

  subroutine return_gradient_hessian(int1,nnz_int1,int2,nnz_int2,den1,nnz_den1,den2,nnz_den2,    &
                           & nfzcpi,ndocpi,nactpi,nextpi,nirrep,orbopt_data,ret_arr,ret_arr_dim)

    ! Subroutine to compute the orbital gradient and the diagonal elements of the orbital hessian
    ! Only unique (i>j) and nonzero (sym_i = sym_j; e-e, d-d always excluded) elements are returned
    ! Provided that dim(ret_arr) >= 2*nrot, where nrot is the number of unique rotation pairs
    !
    !         ret_arr(1:nrot)        = orbital gradient
    !         ret_arr(nrot+1:2*nrot) = orbital hessian
    !
    ! As an example, consider a system with 9 orbitals = 3 doc, 4 act, and 2 ext with 2 possible IRREPs:
    ! 
    !   orbital index 1 2 3 4 5 6 7 8 9
    !   orbital class d d d a a a a e e
    !   orbital IRREP 1 1 2 1 1 1 2 1 2
    !
    ! With a-a rotations included, there are a total of 17 orbital pairs, and the gradient/hessian elements are ordered as
    !    1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17   
    !   1,4 1,5 1,6 1,8 2,4 2,5 2,6 2,8 3,7 3,9 4,5 4,6 4,8 5,6 5,8 6,8 7,9
    ! With a-a rotations excluded, there are a total of 14 orbital pairs, and the gradient/hessian elements are ordered as
    !    1   2   3   4   5   6   7   8   9  10  11  12  13  14 
    !   1,4 1,5 1,6 1,8 2,4 2,5 2,6 2,8 3,7 3,9 4,8 5,8 6,8 7,9 

    ! integer input
    integer, intent(in)     :: nirrep          ! number of irreps in point group
    integer, intent(in)     :: nnz_int1        ! total number of nonzero 1-e integrals
    integer(ip), intent(in) :: nnz_int2        ! total number of nonzero 2-e integrals
    integer, intent(in)     :: nnz_den1        ! total number of nonzero 1-e density elements
    integer, intent(in)     :: nnz_den2        ! total number of nonzero 2-e density elements
    integer, intent(in)     :: ret_arr_dim     ! dimension of the ret_arr array
    integer, intent(in)     :: nfzcpi(nirrep)  ! number of frozen core orbitals
    integer, intent(in)     :: ndocpi(nirrep)  ! number of doubly occupied orbitals per irrep (uses frozen doubly occupied orbitals)
    integer, intent(in)     :: nactpi(nirrep)  ! number of active orbitals per irrep
    integer, intent(in)     :: nextpi(nirrep)  ! number of virtual orbitals per irrep (excluding forzen virtual orbitals) 
    ! real input
    real(wp), intent(inout) :: ret_arr(ret_arr_dim) ! output array with gradient and hessian elements
    real(wp), intent(in)    :: orbopt_data(15)      ! input/output array
    real(wp), intent(in)    :: int1(nnz_int1)       ! nonzero 1-e integral matrix elements
    real(wp), intent(in)    :: int2(nnz_int2)       ! nonzero 2-e integral matrix elements 
    real(wp), intent(in)    :: den1(nnz_den1)       ! nonzero 1-e density matrix elements
    real(wp), intent(in)    :: den2(nnz_den2)       ! nonzero 2-e density matrix elements
    ! local variables
    integer :: error

    ! set up variables based on input orbopt_data
    nthread_use_                = int(orbopt_data(1))
    include_aa_rot_             = int(orbopt_data(2))
    use_exact_hessian_diagonal_ = int(orbopt_data(7))
    df_vars_%use_df_teints      = int(orbopt_data(10))

    ! calculate the total number of orbitals in space
    nfzc_tot_ = sum(nfzcpi)
    ndoc_tot_ = sum(ndocpi)
    nact_tot_ = sum(nactpi)
    next_tot_ = sum(nextpi)
    nmo_tot_  = ndoc_tot_+nact_tot_+next_tot_
    ngem_tot_ = nmo_tot_ * ( nmo_tot_ + 1 ) /2 

    ! allocate indexing arrays
    call allocate_indexing_arrays(nirrep)

    ! determine integral/density addressing arrays
    call setup_indexing_arrays(nfzcpi,ndocpi,nactpi,nextpi)

    ! allocate trans_ but exclude the actual transformation coefficient matrices
    call allocate_transformation_matrices()

    ! set up transformation maps
    call determine_transformation_maps()
 
    ! determine valid orbital rotation pairs (orbital_gradient allocated upon return)
    call setup_rotation_indeces()

    ! allocate temporary Fock matrices
    call allocate_temporary_fock_matrices()

    ! allocate orbital gradient vector and Hessian
    call allocate_return_gradient_hessian() 

    ! set up df mapping arrays if 3-index integrals are used
    error = 1
    if ( df_vars_%use_df_teints == 1 ) then

      error = df_map_setup(nnz_int2)

      ! allocate intermediate matrices for DF integrals

      call allocate_qint()

    end if

    ! save gradient (1:nrot) and hessian (nrot+1 : 2*nrot)

    if ( 2 * rot_pair_%n_tot <= ret_arr_dim ) then

      ! compute orbital gradient
      call orbital_gradient(int1,int2,den1,den2)

      ! compute orbital hessian
      call diagonal_hessian(q_,z_,int2,den1,den2)

! debug
!      call precondition_step(kappa_)
!      orbital_gradient_ = kappa_
! end debug

      ! set up mapping array for sorting array elements
      call setup_rot_pair_map()

      ! save gradient and hessian
      call save_gradient_hessian(ret_arr)

! debug
!      call compute_exponential(kappa_)
!      call print_U()
! end debug

    end if

! **** debug ****
!
!    write(*,*)'new gradient'
!    call orbital_gradient(int1,int2,den1,den2)
!    call print_vector(orbital_gradient_)
!
!    write(*,*)'new hessian'
!    call diagonal_hessian(q_,z_,int2,den1,den2)
!    call print_vector(orbital_hessian_)
!
!    stop
!
! **** end ****

    ! deallocate variables

    if (df_vars_%use_df_teints == 1 ) call deallocate_qint()

    call deallocate_return_gradient_hessian()

    call deallocate_temporary_fock_matrices() 

    call deallocate_trans_partial()

    call deallocate_indexing_arrays()

    return

  end subroutine return_gradient_hessian

  subroutine print_U()

    implicit none

    integer :: isym,i,j,ni

    do isym = 1 , nirrep_

      ni = trans_%nmopi(isym)
      write(*,'(a,i2,1x,a,1x,i2)')' U for irrep: ',isym,' nmo: ',ni

      do i = 1 , ni

        do j = 1 , ni

          write(*,'(2(i2,1x),es24.16)')j,i,trans_%u_irrep_block(isym)%val(j,i)

        end do

      end do

    end do

    return

  end subroutine print_U

  subroutine save_gradient_hessian(array)

    implicit none

    real(wp) :: array(:)

    integer :: i , j , grad_ind, ind

    ind = 0 

    do i = 1 , nmo_tot_

      do j = i + 1 , nmo_tot_

        grad_ind = rot_pair_map(i,j)

        if ( grad_ind == 0 ) cycle

        ind = ind + 1
 
        array(ind) = orbital_gradient_(grad_ind)

      end do

    end do

    do i = 1 , nmo_tot_

      do j = i + 1 , nmo_tot_

        grad_ind = rot_pair_map(i,j)

        if ( grad_ind == 0 ) cycle

        ind = ind + 1

        array(ind) = orbital_hessian_(grad_ind)

      end do

    end do

    return

  end subroutine save_gradient_hessian

  subroutine setup_rot_pair_map()

    implicit none

    integer :: t_sym , a_sym
    integer :: a , i , t , u
    integer :: grad_ind

    rot_pair_map = 0

    ! ext - doc pairs

    do a_sym = 1 , nirrep_

      grad_ind = rot_pair_%pair_offset(a_sym,rot_pair_%ext_doc_type)

      do i = first_index_(a_sym,1) , last_index_(a_sym,1)

        do a = first_index_(a_sym,3) , last_index_(a_sym,3)

          grad_ind = grad_ind + 1

          rot_pair_map(i,a) = grad_ind

        end do

      end do

    end do
   
    ! act - doc pairs

    do t_sym = 1 , nirrep_

      grad_ind = rot_pair_%pair_offset(t_sym,rot_pair_%act_doc_type)

      do i = first_index_(t_sym,1) , last_index_(t_sym,1)

        do t = first_index_(t_sym,2) , last_index_(t_sym,2)

          grad_ind = grad_ind + 1

          rot_pair_map(i,t) = grad_ind

        end do

      end do

    end do

    ! ext - act pairs
 
    do a_sym = 1 , nirrep_

      grad_ind = rot_pair_%pair_offset(a_sym,rot_pair_%ext_act_type)

      do t = first_index_(a_sym,2) , last_index_(a_sym,2)

        do a = first_index_(a_sym,3) , last_index_(a_sym,3)

          grad_ind = grad_ind + 1

          rot_pair_map(t,a) = grad_ind

        end do

      end do

    end do

    ! act - act pairs

    if ( include_aa_rot_ == 1 ) then

      do t_sym = 1 , nirrep_

        grad_ind = rot_pair_%pair_offset(t_sym,rot_pair_%act_act_type)

        do u = first_index_(t_sym,2) , last_index_(t_sym,2)

          do t = u + 1 , last_index_(t_sym,2)

            grad_ind = grad_ind + 1

            rot_pair_map(u,t) = grad_ind

          end do

        end do

      end do

    end if

    return

  end subroutine setup_rot_pair_map

  subroutine allocate_return_gradient_hessian()

    implicit none

    if ( allocated(orbital_gradient_) ) deallocate(orbital_gradient_)
    
    allocate( orbital_gradient_(rot_pair_%n_tot) )

    if ( allocated(orbital_hessian_) ) deallocate(orbital_hessian_)

    allocate( orbital_hessian_(rot_pair_%n_tot) )

    if ( allocated( rot_pair_map ) ) deallocate(rot_pair_map)

    allocate( rot_pair_map( nmo_tot_ , nmo_tot_ ) )

    rot_pair_map = 0

    if ( allocated(kappa_) ) deallocate(kappa_)

    allocate( kappa_(rot_pair_%n_tot) )

    return

  end subroutine allocate_return_gradient_hessian

  subroutine deallocate_return_gradient_hessian()

    implicit none

    if (allocated(rot_pair_%pair_offset))    deallocate(rot_pair_%pair_offset)

    if (allocated(df_vars_%class_to_df_map)) deallocate(df_vars_%class_to_df_map)

    if ( allocated(orbital_gradient_) ) deallocate(orbital_gradient_)

    if ( allocated(orbital_hessian_) ) deallocate(orbital_hessian_)

    if ( allocated(kappa_) ) deallocate(kappa_)

    return

  end subroutine deallocate_return_gradient_hessian

  subroutine allocate_trans_partial()

      implicit none

      integer :: nmo_i_sym,i_sym,i_class

      ! figure out the total number of mos per irrep

      allocate(trans_%npairpi(nirrep_))
      allocate(trans_%nmopi(nirrep_))

      ! figure out number of orbitals per irrep

      trans_%nmopi = 0

      ! loop over orbital symmetries

      do i_sym = 1 , nirrep_

        ! loop over orbital classes

        do i_class = 1 , 3

          trans_%nmopi(i_sym) = trans_%nmopi(i_sym) + ( last_index_(i_sym,i_class) - first_index_(i_sym,i_class) ) + 1

        end do ! end i_class loop

      end do ! end i_sym loop

      ! allocate the remaining indexing arrays

      allocate(trans_%offset(nirrep_))

      allocate(trans_%irrep_to_class_map(nmo_tot_))

      allocate(trans_%class_to_irrep_map(nmo_tot_))
 
      return

  end subroutine allocate_trans_partial

  subroutine deallocate_trans_partial()

      implicit none

      if (allocated(trans_%npairpi))            deallocate(trans_%npairpi)

      if (allocated(trans_%nmopi))              deallocate(trans_%nmopi)

      if (allocated(trans_%offset))             deallocate(trans_%offset)

      if (allocated(trans_%irrep_to_class_map)) deallocate(trans_%irrep_to_class_map)

      if (allocated(trans_%class_to_irrep_map)) deallocate(trans_%class_to_irrep_map)

  end subroutine deallocate_trans_partial

end module focas_gradient_hessian
