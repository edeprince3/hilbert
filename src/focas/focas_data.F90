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

module focas_data
  implicit none

  integer, parameter :: test_full_hessian_aa_ = 0 
  ! *** parameters
  
  integer, parameter :: max_regular_int_=2**30                     ! maximum value for signed integer 
  integer, parameter :: wp = selected_real_kind(10)                ! general working precision kind value
  integer, parameter :: ip = selected_int_kind(16)                 ! 64-bit integers (integral addressing)
  integer, parameter :: fid_ = 12345                               ! file identifier for output file
  integer, parameter :: max_nirrep_=8                              ! maximum number of irreps
  character(3), parameter :: g_element_type_(4) =(/'a-i','e-i','a-a','e-a'/)
  integer, parameter :: group_mult_tab_(max_nirrep_,max_nirrep_) & ! irrep multiplication table
      & = reshape( (/    &  
      & 1,2,3,4,5,6,7,8, &
      & 2,1,4,3,6,5,8,7, &
      & 3,4,1,2,7,8,5,6, &
      & 4,3,2,1,8,7,6,5, &
      & 5,6,7,8,1,2,3,4, &
      & 6,5,8,7,2,1,4,3, &
      & 7,8,5,6,3,4,1,2, &
      & 8,7,6,5,4,3,2,1  /), (/8,8/) )

  ! *** allocatable real arrays

  real(wp), allocatable :: q_(:,:)                                 ! auxiliary matrix (asymmetric, nact*nmo storage) 
  real(wp), allocatable :: z_(:,:)                                 ! auxiliary matrix that contains cotractions of the fock_i with den1 ( nact*nmo storage)
!  real(wp), allocatable :: fock_gen_(:,:)                          ! generalized Fock matrix (nmo*nmo storage)
  real(wp), allocatable :: orbital_gradient_(:)                    ! orbital gradient
  real(wp), allocatable :: dk_(:), sk_(:), tk_(:), rk_(:)          ! nonlinear CG direction, step, and temporary vector
  real(wp), allocatable :: gk_(:), gkp1_(:), yk_(:)                ! nonlinear CG gradient and difference
  real(wp), allocatable :: orbital_hessian_(:)                     ! diagonal elements of the orbital hessian
  real(wp), allocatable :: kappa_(:)                               ! orbital rotation parameters (lt elements of skew-symmetric matrix, npair_ storage)
  real(wp), allocatable :: fa_scr_4index_(:,:)                     ! 4-index integrals for active Fock matrix construction (n^2 * nA^2 storage)
  real(wp), allocatable :: fi_scr_4index_(:,:)                     ! 4-index integrals for active Fock matrix construction (n^2 * nD^2 storage)

  ! *** symmetry data for integrals and densities

  type sym_info
    integer, allocatable     :: ngempi(:)                          ! number of geminals per irrep
    integer, allocatable     :: nnzpi(:)                           ! number of nnz matrix elements
    integer, allocatable     :: offset(:)                          ! offset for first matrix element in this irrep
    integer, allocatable     :: gemind(:,:)                        ! symmetry-reduced index of a geminal 
  end type sym_info

  ! full hessian related data

  type full_hessian_data

    integer, allocatable :: index_map(:,:),npair_sym(:)
    real(wp), allocatable :: aa(:,:),HinvXg_aa(:)

  end type full_hessian_data

  ! symmetry data for transformation

  type matrix_block
    real(wp), allocatable :: val(:,:)                       
  end type matrix_block

  type vector_block
    real(wp), allocatable :: val(:)
  end type vector_block

  type trans_info
    integer, allocatable :: U_eq_I(:)                              ! flag for the type of U matrix (==1 U==I and ==0 U/=I)                       
    integer, allocatable :: npairpi(:)                             ! number of orbital rotation pairs per irrep
    integer, allocatable :: nmopi(:)                               ! number of orbitals per irrep
    integer, allocatable :: offset(:)                              ! first index of orbital (with a given symmetry) in the irrep_to_class_map array
    integer, allocatable :: irrep_to_class_map(:)                  ! mapping array to map symmetry-reduced index to class-index
    integer, allocatable :: class_to_irrep_map(:)                  ! mapping array to map class-index to symmetry-reduced index
    type(matrix_block), allocatable :: u_irrep_block(:)            ! transformation matrix for a symmetry block
  end type trans_info

  type rot_info
    integer :: act_doc_type                                        ! index for active-doubly occupied orbital pair
    integer :: ext_doc_type                                        ! index for external-doubly occupied orbital pair
    integer :: act_act_type                                        ! index for active-active orbital pair
    integer :: ext_act_type                                        ! index for external-active orbital pair 
    integer :: n_tot                                               ! total number of rotation pairs
    integer :: n_ad                                                ! total number of active-doubly occupied rotation pairs
    integer :: n_aa                                                ! total number of active-active pairs
    integer :: n_ed                                                ! total number of external-doubly occupied rotation pairs
    integer :: n_ea                                                ! total number of external-active rotation pairs
    integer, allocatable :: pair_offset(:,:)
  end type rot_info

  type df_info 
    integer :: Qstride                                             ! stride of auxiliary index Q ( == 1 --> (ij,Q=1 ... nQ) // == nQ -> (ij=1...n*(n+1)/2,Q) ) 
    integer :: nQ                                                  !  number of auxiliary function for density-fitted integrals
    integer :: use_df_teints                                       ! flag to use density-fitted 2-e integrals
    integer, allocatable :: class_to_df_map(:)                     ! mapping array to map orbital indeces from class order to df order
!    integer, allocatable :: occgemind(:,:)                         ! symmetry reduced geminal indeces for occupied oritals
!    integer, allocatable :: noccgempi(:)                           ! number of symmetry reduced geminals per irrep
  end type df_info

  type diis_info
    integer :: do_diis                                             ! flag for performing DIIS updates (set internally based on max_num_diis)
    integer :: error                                               ! return code from dgesv()
    integer :: update                                              ! internal flag for performing update
    integer :: max_num_diis                                        ! maximum number of diis vectors stored
    integer :: current_index                                       ! current diis index
    real(wp), allocatable :: B(:,:)                                ! DIIS B matrix (max_num_diis+1,max_num_diis+1)
    real(wp), allocatable :: c(:)                                  ! coefficent vector for DIIS interpolation
    integer, allocatable  :: ip(:)                                 ! temporary matrix used during solution of A * x = c
    real(wp), allocatable :: dP(:,:)                               ! dP vectors (npair,max_num_diis)
    real(wp), allocatable :: P(:,:)
  end type diis_info

  type fock_info
    type(matrix_block), allocatable :: occ(:)
    type(vector_block), allocatable :: ext(:)
  end type fock_info

  type qint_info
    type(matrix_block), allocatable :: tuQ(:) 
  end type qint_info

  type gen_f_info
    type(matrix_block), allocatable :: doc(:)
    type(matrix_block), allocatable :: act(:)  
    type(matrix_block), allocatable :: ext(:)
    type(vector_block), allocatable :: doc_e(:)
    type(vector_block), allocatable :: act_e(:)
    type(vector_block), allocatable :: ext_e(:)
  end type gen_f_info
 
  type p_sym_info
    type(matrix_block), allocatable :: t_sym(:)
  end type p_sym_info

  type fock_scr_info
    type(p_sym_info), allocatable :: p_sym(:)
  end type fock_scr_info

  ! *** full hessian data

  type(full_hessian_data) :: full_orbital_hessian_

  ! *** allocatable derived types

  type(sym_info)      :: dens_                                        ! density symmetry data
  type(sym_info)      :: ints_                                        ! integral symmetry data
  type(trans_info)    :: trans_
  type(diis_info)     :: diis_
  type(fock_info)     :: fock_i_ 
  type(fock_info)     :: fock_a_
  type(qint_info)     :: qint_
  type(gen_f_info)    :: gen_f_
  type(fock_scr_info) :: fa_scr_(3)
  type(fock_scr_info) :: fi_scr_(3)

  ! indexing derived types
  
  type(rot_info)   :: rot_pair_                                    ! info for rotation pair indexing
  type(df_info)    :: df_vars_

  ! *** allocatable orbital index arrays

  integer, allocatable :: nfzcpi_(:)                               ! number of frozen (not optimized) doubly occupied orbitals per irrep
  integer, allocatable :: ndocpi_(:)                               ! number of doubly occupied orbitals (including frozen) per irrep
  integer, allocatable :: nactpi_(:)                               ! number of active orbitals per irrep
  integer, allocatable :: nextpi_(:)                               ! number of external orbitals per irrep  
  integer, allocatable :: first_index_(:,:)                        ! index of first orbital in this class and irrep (nirrep,3)
  integer, allocatable :: last_index_(:,:)                         ! index of last orbital in this class and irrep (nirrep,3)
  integer, allocatable :: orb_sym_scr_(:)                          ! scratch array to store the symmetries of orbitals

  ! *** integers

  integer :: nirrep_                                               ! number of irreps in point group
  integer :: nfzc_tot_                                             ! number of frozen doubly occupied orbitals
  integer :: ndoc_tot_                                             ! total number of doubly occupied (including frozen) orbitals
  integer :: nact_tot_                                             ! total number of active orbitals
  integer :: next_tot_                                             ! total number of external orbitals
  integer :: nmo_tot_                                              ! total number of orbitals
  integer :: ngem_tot_                                             ! total number of geminals (ij) with i <= j
  integer :: include_aa_rot_                                       ! 1/0 = include/do not include rotations between active-active orbtials
  integer :: nthread_use_                                          ! number of threads to use in parallel parts of the code (this is the actuaal number of threads used)
  integer :: log_print_                                            ! 1/0 = flag for printing iteration/info for orbtial optimization
  integer :: num_negative_diagonal_hessian_                        ! number of negative diagonal Hessian matrix elements
  integer :: use_exact_hessian_diagonal_                           ! flag to use exact expressions for the diagonal elements of the Hessian
  integer :: num_diis_vectors_
 
  ! *** doubles
  real(wp) :: e1_c_                                                ! core contribution to 1-e energy
  real(wp) :: e2_cc_                                               ! core contribution to 2-e energy all indeces in g(pq|rs) in \D
  real(wp) :: e1_a_                                                ! active contribution to 1-e energy 
  real(wp) :: e2_aa_                                               ! active-active contribution to 2-e energy all indeces in g(pq|rs) in \A
  real(wp) :: e2_ca_                                               ! core-active contribution to 2-e energy only 2 indeces in g(pq|rs) in \A
  real(wp) :: e_nuc_rep_                                           ! nuclear repulsion energy
  real(wp) :: e_frozen_core_                                       ! frozen core energy
  real(wp) :: e_total_                                             ! total energy
  real(wp) :: e1_total_                                            ! total 1-e energy
  real(wp) :: e2_total_                                            ! total 2-e energy
  real(wp) :: e_active_                                            ! active space energy
  real(wp) :: grad_norm_                                           ! norm of the gradient ddot(g,g)
  real(wp) :: min_diag_hessian_                                    ! smallest diagonal Hessian element

  real(wp) :: max_grad_val_                                        ! largest gradient element
  real(wp) :: norm_grad_large_                                     ! total norm of large gradient elements 
  integer :: max_grad_ind_(2)                                      ! orbitalindeces for largest gradient element
  integer :: max_grad_sym_                                         ! orbital pair symmetry
  integer :: max_grad_typ_                                         ! type of orbital rotation
  integer :: n_grad_large_                                         ! number of large gradient elements (val <= +/- 0.75_*max_grad_val) 

  contains

    pure function pq_index(i,j)
! this function computes the two-electron index index (lower triangular reference)
! index = ii*(ii-1)/2+jj where ii=max(i,j) and jj=min(i,j)
! the ishft(k,-1) divides the value of the integer k by 2 and seems to be somewhat
! faster than the regular human-readable expression
      implicit none
      integer, intent(in) ::i,j
      integer(ip) :: pq_index
      if (i.ge.j) then
        pq_index=ishft(i*(i-1),-1)+j
        return
      else
        pq_index=ishft(j*(j-1),-1)+i
        return
      end if
    end function pq_index

    pure function df_aa_index(g,a,a_sym)
! function to return the column index of df(:,ga) where 
! both g & a are an active orbitals (LT storage)
      integer, intent(in)  :: g,a,a_sym
      integer  :: a_i,g_i,df_aa_index

      ! adjust for the number of doubly-ococcupied orbital in this irrep
      a_i = trans_%class_to_irrep_map(a)-ndocpi_(a_sym)

      ! orbital index within irrep
      g_i = trans_%class_to_irrep_map(g)-ndocpi_(a_sym)

      if (a_i.ge.g_i) then
        df_aa_index=ishft(a_i*(a_i-1),-1)+g_i
        return
      else
        df_aa_index=ishft(g_i*(g_i-1),-1)+a_i
        return
      end if

      return
    end function df_aa_index

    pure function df_ga_index(g,a,a_sym)
! function to return the column index of df(:,ga) where 
! g is a general index and a is an active index
! assumes that for each general index g, all the a indeces are stored in contiguous order
      integer, intent(in)  :: g,a,a_sym
      integer  :: a_i,g_i,df_ga_index

      ! adjust for the number of doubly-ococcupied orbital in this irrep
      a_i = trans_%class_to_irrep_map(a)-ndocpi_(a_sym)

      ! orbital index within irrep
      g_i = trans_%class_to_irrep_map(g)


      df_ga_index = ( g_i - 1 ) * nactpi_(a_sym) + a_i

      return
    end function df_ga_index

    pure function df_gd_index(g,d,d_sym)
! function to return the column index of df(:,ga) where 
! g is a general index and d is an doubly-occupied index
! assumes that for each general index g, all the d indeces are stored in contiguous order
      integer, intent(in)  :: g,d,d_sym
      integer  :: d_i,g_i,df_gd_index

      ! adjust for the number of doubly-ococcupied orbital in this irrep
      d_i = trans_%class_to_irrep_map(d)

      ! orbital index within irrep
      g_i = trans_%class_to_irrep_map(g)


      df_gd_index = ( g_i - 1 ) * ndocpi_(d_sym) + d_i

      return
    end function df_gd_index

    pure function df_pq_index(i,j)
! this function computes the two-electron index index (lower triangular reference)
! index = ii*(ii+1)/2+jj where ii=max(i,j) and jj=min(i,j)
! the ishft(k,-1) divides the value of the integer k by 2 and seems to be somewhat
! faster than the regular human-readable expression
      implicit none
      integer, intent(in) ::i,j
      integer(ip) :: i_long,j_long
      integer(ip) :: df_pq_index
      i_long = int(i,kind=ip)
      j_long = int(j,kind=ip)
      if (i_long.ge.j_long) then
        df_pq_index = i_long * ( i_long + 1 ) / 2 + j_long
        if ( df_vars_%Qstride == 1 ) df_pq_index = df_pq_index * int(df_vars_%nQ,kind=ip)
        return
      else
        df_pq_index= j_long * ( j_long + 1 ) / 2 + i_long
        if ( df_vars_%Qstride == 1 ) df_pq_index = df_pq_index * int(df_vars_%nQ,kind=ip)
        return
      end if
    end function df_pq_index

    function timer()
      real(wp) :: omp_get_wtime
      real(wp) :: timer(2)
      ! actual time
#ifdef OMP 
      timer(1) = omp_get_wtime()
#else
      call cpu_time(timer(1))
#endif
      ! total CPU time
      call cpu_time(timer(2))

    end function timer

    function my_ddot(n,vec_1,stride_1,vec_2,stride_2)

      ! simple wrapper for LAPACKs DDOT to avoid integer overflow
      real(wp) :: my_ddot

      integer, intent(in)  :: n,stride_1,stride_2
      real(wp), intent(in) :: vec_1(:),vec_2(:)

      integer(ip)          :: s_1,s_2,e_1,e_2,inc_1,inc_2
      integer              :: n_pass,n_have,n_need,ddot_pass
      real(wp)             :: ddot

!      my_ddot = ddot(n,vec_1,stride_1,vec_2,stride_2)
!      return

      ! number of values accumulated in one pass   
      if ( stride_1 > stride_2 ) then

        n_pass = max_regular_int_ / stride_1

        if ( mod( max_regular_int_ , stride_1) == 0 ) n_pass = n_pass + 1

      else

        n_pass = max_regular_int_ / stride_2
        
        if ( mod( max_regular_int_ , stride_2) == 0 ) n_pass = n_pass + 1

      end if 
      ! increment in starting index
      inc_1   = n_pass * stride_1 ; inc_2   = n_pass * stride_2

      ! initialize
      n_have  = 0 ; my_ddot = 0.0_wp ; s_1 = 1 ; s_2 = 1 

      ! repeated calls to ddot
      do ddot_pass = 1 , n / n_pass

        n_have = n_have + n_pass
 
        e_1 = s_1 + inc_1 - 1 ; e_2 = s_2 + inc_2 - 1

        my_ddot = my_ddot + ddot(n_pass,vec_1(s_1:e_1),stride_1,vec_2(s_2:e_2),stride_2)

        s_1 = s_1 + inc_1 ; s_2 = s_2 + inc_2

      end do

      n_need = n - n_have

      if ( n_need == 0 ) return

      e_1 = s_1 + ( n_need - 1 ) * stride_1 ; e_2 = s_2 + ( n_need -1 ) * stride_2

      my_ddot = my_ddot + ddot(n-n_have,vec_1(s_1:e_1),stride_1,vec_2(s_2:e_2),stride_2)

      return

    end function my_ddot

    subroutine my_dcopy(n,vec_x,stride_x,vec_y,stride_y)

      ! simple wrapper for LAPACKs DDOT to avoid integer overflow
      ! assumes that stride_x >= stride_y

      integer, intent(in)  :: n,stride_x,stride_y
      real(wp), intent(in) :: vec_x(:),vec_y(:)

      integer(ip)          :: s_x,s_y,e_x,e_y,inc_x,inc_y
      integer              :: n_pass,n_have,n_need,dcopy_pass

!      call dcopy(n,vec_x,stride_x,vec_y,stride_y)
!      return

      ! number of values accumulated in one pass   
      if ( stride_x > stride_y ) then

        n_pass = max_regular_int_ / stride_x

        if ( mod( max_regular_int_ , stride_x) == 0 ) n_pass = n_pass + 1

      else

        n_pass = max_regular_int_ / stride_y
        
        if ( mod( max_regular_int_ , stride_y) == 0 ) n_pass = n_pass + 1

      end if

      ! increment in starting index
      inc_x   = n_pass * stride_x ; inc_y   = n_pass * stride_y

      ! initialize
      n_have  = 0 ; s_x = 1 ; s_y = 1

      ! repeated calls to ddot
      do dcopy_pass = 1 , n / n_pass

        n_have = n_have + n_pass

        e_x = s_x + inc_x - 1 ; e_y = s_y + inc_y - 1

        call dcopy(n_pass,vec_x(s_x:e_x),stride_x,vec_y(s_y:e_y),stride_y)

        s_x = s_x + inc_x ; s_y = s_y + inc_y

      end do

      n_need = n - n_have

      if ( n_need == 0 ) return

      e_x = s_x + ( n_need - 1 ) * stride_x ; e_y = s_y + ( n_need - 1 ) * stride_y

      call dcopy(n-n_have,vec_x(s_x:e_x),stride_x,vec_y(s_y:e_y),stride_y)

      return

    end subroutine my_dcopy

    subroutine my_daxpy(n,scale_x,vec_x,stride_x,vec_y,stride_y)

      ! simple wrapper for LAPACKs DDOT to avoid integer overflow
      ! assumes that stride_x >= stride_y

      integer, intent(in)  :: n,stride_x,stride_y
      real(wp), intent(in) :: vec_x(:),vec_y(:)
      real(wp), intent(in) :: scale_x

      integer(ip)          :: s_x,s_y,e_x,e_y,inc_x,inc_y
      integer              :: n_pass,n_have,n_need,dcopy_pass
  
!      call daxpy(n,scale_x,vec_x,stride_x,vec_y,stride_y)
!      return

      ! number of values accumulated in one pass   
      if ( stride_x > stride_y ) then

        n_pass = max_regular_int_ / stride_x

        if ( mod( max_regular_int_ , stride_x) == 0 ) n_pass = n_pass + 1

      else

        n_pass = max_regular_int_ / stride_y

        if ( mod( max_regular_int_ , stride_y) == 0 ) n_pass = n_pass + 1

      end if

      ! increment in starting index
      inc_x   = n_pass * stride_x ; inc_y   = n_pass * stride_y

      ! initialize
      n_have  = 0 ; s_x = 1 ; s_y = 1

      ! repeated calls to ddot
      do dcopy_pass = 1 , n / n_pass

        n_have = n_have + n_pass

        e_x = s_x + inc_x - 1 ; e_y = s_y + inc_y - 1

        call daxpy(n_pass,scale_x,vec_x(s_x:e_x),stride_x,vec_y(s_y:e_y),stride_y)

        s_x = s_x + inc_x ; s_y = s_y + inc_y

      end do

      n_need = n - n_have

      if ( n_need == 0 ) return

      e_x = s_x + ( n_need - 1 ) * stride_x ; e_y = s_y + ( n_need - 1 ) * stride_y

      call daxpy(n_need,scale_x,vec_x(s_x:e_x),stride_x,vec_y(s_y:e_y),stride_y)

      return

    end subroutine my_daxpy

    subroutine my_dscal(n,scale_fac,vec_x,stride_x)

      ! simple wrapper for LAPACKs DDOT to avoid integer overflow
      ! assumes that stride_x >= stride_y

      integer, intent(in)  :: n,stride_x
      real(wp), intent(in) :: vec_x(:)
      real(wp), intent(in) :: scale_fac

      integer(ip)          :: s_x,e_x,inc_x
      integer              :: n_pass,n_have,n_need,dscal_pass

!      call dscal(n,scale_fac,vec_tmp,stride_x)
!      return

      ! number of values accumulated in one pass   

      n_pass = max_regular_int_ / stride_x

      if ( mod( max_regular_int_ , stride_x) == 0 ) n_pass = n_pass + 1

      ! increment in starting index
      inc_x   = n_pass * stride_x 

      ! initialize
      n_have  = 0 ; s_x = 1 

      ! repeated calls to ddot
      do dscal_pass = 1 , n / n_pass

        n_have = n_have + n_pass

        e_x = s_x + inc_x - 1

        call dscal(n_pass,scale_fac,vec_x(s_x:e_x),stride_x)

        s_x = s_x + inc_x 

      end do

      n_need = n - n_have

      if ( n_need == 0 ) return

      e_x = s_x + ( n_need - 1 ) * stride_x

      call dscal(n_need,scale_fac,vec_x(s_x:e_x),stride_x)

      return

    end subroutine my_dscal

    subroutine abort_print(error_code)

      implicit none

      integer, intent(in) :: error_code

      if (error_code == 10) write(*,'(a)')'error encountered in function gather_kappa_block()'

      if (error_code == 11) write(*,'(a)')'error encountered in function compute_block_exponential()'

      if (error_code == 20) write(*,'(a)')'error encountered in function df_map_setup()'

      if (error_code == 30) write(*,'(a)')'error encountered in function transform_oeints()'

      if (error_code == 31) write(*,'(a)')'error encountered in function transform_teints()'

      if (error_code == 32) write(*,'(a)')'error encountered in function transform_mocoeff()'

      if (error_code == 311) write(*,'(a)')'error encountered in function transform_teints_irrep_block()'

      if (error_code == 312) write(*,'(a)')'error encountered in function allocate_transform_scr()'

      if (error_code == 313) write(*,'(a)')'error encountered in function transform_teints_g0_block()'

      if (error_code == 314) write(*,'(a)')'error encountered in function allocate_transform_scr_g0()' 

      if (error_code == 40) write(*,'(a)')'error encountered in function diagonalize_opdm()'

      if (error_code == 50) write(*,'(a)')'error encountered in function transform_mocoeff()'

      if (error_code == 510) write(*,'(a)')'error encountered in function precompute_coulomb()'

      if (error_code == 511) write(*,'(a)')'error encountered in function compute_gen_fock_block_df() for inactive'

      if (error_code == 512) write(*,'(a)')'error encountered in function compute_gen_fock_block_df() for active'

      if (error_code == 513) write(*,'(a)')'error encountered in function compute_gen_fock_block_df() for external'

      if (error_code == 521) write(*,'(a)')'error encountered in function compute_gen_fock_block() for inactive'

      if (error_code == 522) write(*,'(a)')'error encountered in function compute_gen_fock_block() for active'

      if (error_code == 523) write(*,'(a)')'error encountered in function compute_gen_fock_block() for external'

      if (error_code == 531) write(*,'(a)')'error encountered in function copy_semicanonical_mos_block() for inactive'

      if (error_code == 532) write(*,'(a)')'error encountered in function copy_semicanonical_mos_block() for active'

      if (error_code == 533) write(*,'(a)')'error encountered in function copy_semicanonical_mos_block() for external'

      if (error_code == 541) write(*,'(a)')'error encountered in function diagonalize_gen_fock_block() for inactive'

      if (error_code == 542) write(*,'(a)')'error encountered in function diagonalize_gen_fock_block() for active'

      if (error_code == 543) write(*,'(a)')'error encountered in function diagonalize_gen_fock_block() for external'

      if (error_code == 544) write(*,'(a)')'error encountered in function diagonalize_opdm_block()'

      if (error_code == 545) write(*,'(a)')'error encountered in function print_orbital_energies()'

      stop

      return

    end subroutine abort_print

end module focas_data
