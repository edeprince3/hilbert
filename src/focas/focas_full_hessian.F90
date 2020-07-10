module focas_full_hessian

  use focas_data
  use focas_driver
  use focas_gradient
  use focas_hessian
  use focas_exponential

  contains

  subroutine return_full_hessian(int1,nnz_int1,int2,nnz_int2,den1,nnz_den1,den2,nnz_den2,    &
                           & nfzcpi,ndocpi,nactpi,nextpi,nirrep,orbopt_data,ret_arr,ret_arr_dim)

    ! Subroutine to return all elements of the active-active orbital hessian
    ! Only lower-triangular elements (pq>=rs) of the Hessian are returned
    ! Composite Hessian indeces pq are such that p<q and active orbitals are ordered accoring to irrep
    ! For example, with 2 irreps and 3 orbitals in the 1st irrep (indexed 1- 3) and 4 in the 2nd irrep 
    ! (indexed 4 - 7), the composite hessian indeces pq would be:
    !  p:  1  1  2  4  4  4  5  5  6 
    !  q:  2  3  3  5  6  7  6  7  7
    ! pq:  1  2  3  4  5  6  7  8  9

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
    integer :: dim_hessian

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
    call allocate_return_full_hessian() 

    ! set up df mapping arrays if 3-index integrals are used
    error = 1
    if ( df_vars_%use_df_teints == 1 ) then

      error = df_map_setup(nnz_int2)

      ! allocate intermediate matrices for DF integrals

      call allocate_qint()

    end if

    ! save gradient (1:nrot) and hessian (nrot+1 : 2*nrot)

    dim_hessian = rot_pair_%n_aa

    if ( ( dim_hessian > 0 ) .and. & 
       & ( dim_hessian * ( dim_hessian + 1 ) / 2 <=ret_arr_dim ) ) then

      ! calculate orbital gradient (only the F_i, q, and z matrices are needed)
      call orbital_gradient(int1,int2,den1,den2) 
 
      ! calculate hessian matrix
      error=full_hessian_aa(fock_i_%occ,q_,z_,int2,den1,den2)

!      ! compute orbital hessian
!      call diagonal_hessian(q_,z_,int2,den1,den2)

      ! save gradient and hessian
      call save_full_hessian(ret_arr)

    end if

    ! deallocate variables

    if (df_vars_%use_df_teints == 1 ) call deallocate_qint()

    call deallocate_return_full_hessian()

    call deallocate_temporary_fock_matrices() 

    call deallocate_trans_partial_loc()

    call deallocate_indexing_arrays()

    return

  end subroutine return_full_hessian

  subroutine save_full_hessian(array)

    implicit none

    real(wp) :: array(:)

    integer :: p,q,r,s,pq,rs,r_min,s_min
    integer  :: h_p,h_r    

    integer  :: array_index

    array_index = 0 

    do h_p = 1 , nirrep_

      do p = first_index_(h_p,2) , last_index_(h_p,2)

        do q = p + 1 , last_index_(h_p,2) 

          pq = full_orbital_hessian_%index_map(p,q)

          do h_r = h_p , nirrep_

            r_min = first_index_(h_r,2)
            if ( h_r == h_p ) r_min = p

            do r = r_min , last_index_(h_r,2)

              s_min = r + 1
              if ( r == p ) s_min = q

              do s = s_min , last_index_(h_r,2)

                rs = full_orbital_hessian_%index_map(r,s)
                
                array_index = pq_index(pq,rs)

                array(array_index) = full_orbital_hessian_%aa(rs,pq)

              end do

            end do
 
           end do

        end do

      end do

    end do

    return

  end subroutine save_full_hessian

  subroutine allocate_return_full_hessian()

    implicit none

    if ( allocated(orbital_gradient_) ) deallocate(orbital_gradient_)
    
    allocate( orbital_gradient_(rot_pair_%n_tot) )

    if ( allocated(orbital_hessian_) ) deallocate(orbital_hessian_)

    allocate( orbital_hessian_(rot_pair_%n_tot) )

    if ( allocated(kappa_) ) deallocate(kappa_)

    allocate( kappa_(rot_pair_%n_tot) )

    return

  end subroutine allocate_return_full_hessian

  subroutine deallocate_return_full_hessian()

    implicit none

    integer :: error

    if (allocated(rot_pair_%pair_offset))    deallocate(rot_pair_%pair_offset)

    if (allocated(df_vars_%class_to_df_map)) deallocate(df_vars_%class_to_df_map)

    if ( allocated(orbital_gradient_) ) deallocate(orbital_gradient_)

    if ( allocated(orbital_hessian_) ) deallocate(orbital_hessian_)

    error = deallocate_full_hessian_data()
 
    if ( allocated(kappa_) ) deallocate(kappa_)

    return

  end subroutine deallocate_return_full_hessian

  subroutine allocate_trans_partial_loc()

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

  end subroutine allocate_trans_partial_loc

  subroutine deallocate_trans_partial_loc()

      implicit none

      if (allocated(trans_%npairpi))            deallocate(trans_%npairpi)

      if (allocated(trans_%nmopi))              deallocate(trans_%nmopi)

      if (allocated(trans_%offset))             deallocate(trans_%offset)

      if (allocated(trans_%irrep_to_class_map)) deallocate(trans_%irrep_to_class_map)

      if (allocated(trans_%class_to_irrep_map)) deallocate(trans_%class_to_irrep_map)

  end subroutine deallocate_trans_partial_loc

end module focas_full_hessian
