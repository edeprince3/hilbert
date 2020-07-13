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

module focas_driver
  use focas_data
  use focas_energy
  use focas_gradient
  use focas_hessian
  use focas_transform_driver
  use focas_exponential
  use focas_redundant
  use focas_diis

  implicit none

  contains

  subroutine focas_optimize(mo_coeff,int1,nnz_int1,int2,nnz_int2,den1,nnz_den1,den2,nnz_den2,    &
                           & nfzcpi,ndocpi,nactpi,nextpi,nirrep,orbopt_data,fname)
 
    ! integer input
    integer, intent(in)     :: nirrep          ! number of irreps in point group
    integer, intent(in)     :: nnz_int1        ! total number of nonzero 1-e integrals
    integer(ip), intent(in) :: nnz_int2        ! total number of nonzero 2-e integrals
    integer, intent(in)     :: nnz_den1        ! total number of nonzero 1-e density elements
    integer, intent(in)     :: nnz_den2        ! total number of nonzero 2-e density elements
    integer, intent(in)     :: nfzcpi(nirrep)  ! number of frozen (not optimized) doubly occupied orbitals per irrep
    integer, intent(in)     :: ndocpi(nirrep)  ! number of doubly occupied orbitals per irrep (includes frozen doubly occupied orbitals)
    integer, intent(in)     :: nactpi(nirrep)  ! number of active orbitals per irrep
    integer, intent(in)     :: nextpi(nirrep)  ! number of virtual orbitals per irrep (excluding forzen virtual orbitals) 
    ! real input
    real(wp), intent(inout) :: orbopt_data(15) ! input/output array
    real(wp), intent(inout) :: mo_coeff(:,:)   ! mo coefficient matrix
    real(wp), intent(in)    :: int1(nnz_int1)  ! nonzero 1-e integral matrix elements
    real(wp), intent(in)    :: int2(nnz_int2)  ! nonzero 2-e integral matrix elements 
    real(wp), intent(in)    :: den1(nnz_den1)  ! nonzero 1-e density matrix elements
    real(wp), intent(in)    :: den2(nnz_den2)  ! nonzero 2-e density matrix elements
    ! character
    character(*), intent(in) :: fname          ! name of file to print output
    ! step size control parameters
    real(wp), parameter     :: r_increase_tol=0.75_wp  ! dE ratio above which step size is increased
    real(wp), parameter     :: r_decrease_tol=0.25_wp  ! dE ratio below which step size is reduced
    real(wp), parameter     :: r_increase_fac=1.20_wp  ! factor by which to increase step size
    real(wp), parameter     :: r_decrease_fac=0.70_wp  ! factor by which to reduce the step size 

    ! quadratic model variables
    real(wp) :: coeff_a,coeff_b,coeff_c

    ! nonlinear CG variables
    real(wp) :: alpha, beta, gamma, delta, tau, theta

    ! timing variables
    real(wp) :: t0(2),t1(2),t_ene,t_gh,t_exp,t_wall_trans,t_cpu_trans,t_wall_aux,t_cpu_aux

    ! iteration variables
    real(wp) :: current_energy,last_energy,delta_energy,gradient_norm_tolerance,delta_energy_tolerance
    real(wp) :: initial_energy,delta_energy_approximate
    integer  :: i,iter,max_iter,error,converged

    ! variables for trust radius
    integer :: evaluate_gradient,reject
    real(wp) :: step_size,step_size_update,step_size_factor,e_new,e_init,de_ratio,e_tmp
    character :: reject_char(1)

    ! other variables
    logical :: fexist
    integer :: alg

    ! set up variables based on input orbopt_data
    nthread_use_                = int(orbopt_data(1))
    include_aa_rot_             = int(orbopt_data(2))
    gradient_norm_tolerance     = orbopt_data(4)
    delta_energy_tolerance      = orbopt_data(5)  
    log_print_                  = int(orbopt_data(6))
    use_exact_hessian_diagonal_ = int(orbopt_data(7))
    diis_%max_num_diis          = int(orbopt_data(8))
    max_iter                    = int(orbopt_data(9)) 
    df_vars_%use_df_teints      = int(orbopt_data(10))
    alg                         = int(orbopt_data(15))

    if ( log_print_ == 1 ) then
      inquire(file=fname,exist=fexist)
      if (fexist) then
        open(fid_,file=fname,status='old',position='append')
      else
        open(fid_,file=fname,status='new')
      endif
    endif

    ! calculate the total number of orbitals in space
    nfzc_tot_ = sum(nfzcpi)
    ndoc_tot_ = sum(ndocpi)
    nact_tot_ = sum(nactpi)
    next_tot_ = sum(nextpi)
    nmo_tot_  = ndoc_tot_+nact_tot_+next_tot_
    ngem_tot_ = nmo_tot_ * ( nmo_tot_ + 1 ) / 2

    ! allocate indexing arrays
    call allocate_indexing_arrays(nirrep)

    ! determine integral/density addressing arrays
    call setup_indexing_arrays(nfzcpi,ndocpi,nactpi,nextpi)

    ! allocate transformation matrices
    call allocate_transformation_matrices() 
   
    ! determine valid orbital rotation pairs (orbital_gradient allocated upon return)
    call setup_rotation_indeces()

    ! allocate temporary Fock matrices
    call allocate_temporary_fock_matrices()

    ! allocate matrices for DIIS extrapolation
    call allocate_diis_data()

    ! allocate remaining arrays/matrices
    call allocate_initial()

    ! determine indexing arrays (needed for sorts in the integral transformation step)
    call determine_transformation_maps()

    ! check for numerically doubly-occupied or empty orbitals
    call compute_opdm_nos(den1)

    ! set up df mapping arrays if density-fitted 2-e integrals are used
    error = 0
    if ( df_vars_%use_df_teints == 1 ) error = df_map_setup(nnz_int2)
    if ( error /= 0 ) call abort_print(20)

    ! allocate intermediate matrices for DF integrals
    if ( df_vars_%use_df_teints == 1 ) call allocate_qint()

    ! **********************************************************************************
    ! *** at this point, everything is allocated and we are ready to do the optimization
    ! **********************************************************************************

    ! initialize trust radius variables
    evaluate_gradient = 1
    step_size_factor  = 1.0_wp

    if ( log_print_ == 1 ) then

      write(fid_,*)

      write(fid_,'((a4,1x),(a16,1x),(a10),(2x),2(a10,1x),1(a4,1x),(a3,1x),(a9,1x),(a8,1x),2(a11,1x),2(a11,1x))') &
                & 'iter','E(k)','dE','||g||','max|g|','type','sym','(i,j)','R_step',     &
                & 'twall_trans','tcpu_trans','twall_aux','tcpu_other'

      write(fid_,'(a)')('-----------------------------------------------------------------&
                       & -----------------------------------------------------------------')

    end if

    ! calculate current energy and gradient
    ! frozen doubly occupied orbitals zeroed in gradient calculation

    call compute_energy(int1,int2,den1,den2)
    call orbital_gradient(int1,int2,den1,den2)
    call diagonal_hessian(q_,z_,int2,den1,den2)

    initial_energy = e_total_
    last_energy    = e_total_
    converged      = 0
    iter           = 0

    gk_ =  orbital_gradient_
    dk_ = -orbital_gradient_
    dk_ = dk_ / abs(orbital_hessian_)

    do
      ! store energy at initial point
      e_init = e_total_
       
      ! check for descent and save terms for quadratic model
      coeff_c = e_total_
      coeff_b = my_ddot(rot_pair_%n_tot,dk_,1,gk_,1)
      if (coeff_b >= 0.0_wp) then
        
        ! not a descent direction; use preconditioned steepest descent
        ! write(*,*) 'Not descent; reset'

        dk_ = -gk_ / abs(orbital_hessian_)
        coeff_b = my_ddot(rot_pair_%n_tot,dk_,1,gk_,1)

        if (coeff_b >= 0.0_wp) then
          ! still not a descent direction; use steepest descent
          ! write(*,*) 'Preconditioner indefinite; gradient'

          orbital_hessian_ = 1.0_wp
          dk_ = -gk_
          coeff_b = my_ddot(rot_pair_%n_tot,dk_,1,gk_,1)
        end if 
      end if

      !!!!! Determine the steplength !!!!!

      ! compute the step and determine the steplength

      alpha = 1.0_wp
      sk_ = alpha*dk_

      ! compute transformation matrix
      call compute_exponential(sk_)

      ! transform the integrals; xk = xk+alpha*dk
      call transform_driver(int1,int2,mo_coeff)
      
      ! calculate the current energy

      call compute_energy(int1,int2,den1,den2)

      ! finish quadratic model

      coeff_a = 2.0_wp*((e_total_ - coeff_c) - coeff_b*alpha) / (alpha*alpha)

      ! If the current step satisfies the Armijo conditions, then accept the step
      ! Otherwise, if the model is a strictly convex quadratic, go to the minimizer
      ! Otherwise accept the full step

      if (e_total_ <= coeff_c + 1e-4*alpha*coeff_b) then

        ! The step satisfies the Armijo condition

        call orbital_gradient(int1,int2,den1,den2)
	! write(*,*), 'Decrse: ', coeff_a, coeff_b, coeff_c, alpha, e_total_, grad_norm_

      else if (coeff_a > 0.0_wp) then

        ! The quadratic is strictly convex; the minimizer is -coeff_b / coeff_a
        ! To evaluate the function we need to step back alpha and step forward to
        ! the new value; the following does this and then fixes up the remaining entries
        
        alpha = -coeff_b / coeff_a - alpha
        sk_ = alpha*sk_

        ! compute transformation matrix
        call compute_exponential(sk_)

        ! transform the integrals; xk = xk+sk_
        call transform_driver(int1,int2,mo_coeff)

        ! calculate the current energy

        call compute_energy(int1,int2,den1,den2)
        call orbital_gradient(int1,int2,den1,den2)

        alpha = -coeff_b / coeff_a
        sk_ = alpha*dk_
        if (e_total_ <= coeff_c + 1e-4*alpha*coeff_b) then
          ! write(*,*), 'CnvexG: ', coeff_a, coeff_b, coeff_c, alpha, e_total_, grad_norm_
        else 
          ! write(*,*), 'CnvexB: ', coeff_a, coeff_b, coeff_c, alpha, e_total_, grad_norm_
        end if
      else
        call orbital_gradient(int1,int2,den1,den2)
	! write(*,*), 'Cncave: ', coeff_a, coeff_b, coeff_c, alpha, e_total_, grad_norm_
      end if

      e_new = e_total_
      delta_energy = e_new - e_init
      
      iter = iter + 1
      if (( abs(delta_energy) <= delta_energy_tolerance ) .and. (grad_norm_ <= gradient_norm_tolerance)) then
        converged = 1
        exit
      else if (iter >= max_iter) then
        converged = 0
        exit
      end if
           
      !!!!! We have accepted a step and have the gradient and can proceed with the
      !!!!! nonlinear conjugate gradient step computation; sk and dk are correct

      !!!!! Update preconditioner; better results without updating
      !! call diagonal_hessian(q_,z_,int2,den1,den2)

      if (alg == 0) then
        !!!!! Preconditioned Steepest Descent
        gk_ =  orbital_gradient_
        dk_ = -orbital_gradient_ / abs(orbital_hessian_)

      else if (alg == 1) then
        !!!!! Preconditioned Hestenes-Stiefel 1952 Method
        gkp1_ = orbital_gradient_
        yk_   = gkp1_ - gk_
        sk_   = orbital_gradient_ / abs(orbital_hessian_)

        tau = my_ddot(rot_pair_%n_tot,dk_,1,yk_,1)
        gamma = my_ddot(rot_pair_%n_tot,yk_,1,sk_,1) 
        beta = gamma / tau

        gk_ = gkp1_
        dk_ = -sk_ + beta*dk_

      else if (alg == 2) then
        !!!!! Preconditioned Dai-Yuan 1999 Method 
        gkp1_ = orbital_gradient_
        yk_   = gkp1_ - gk_
        sk_   = orbital_gradient_ / abs(orbital_hessian_)

        tau = my_ddot(rot_pair_%n_tot,dk_,1,yk_,1)
        gamma = my_ddot(rot_pair_%n_tot,gkp1_,1,sk_,1) 
        beta = gamma / tau
   
        gk_ = gkp1_
        dk_ = -sk_ + beta*dk_

      else if (alg == 3) then
        !!!!! Preconditioned Hager-Zhang 2005 Method; Theta=1.0 per Dai-Kou 2013
        theta = 1.0_wp
        gkp1_ = orbital_gradient_
        yk_   = gkp1_ - gk_
        sk_   = orbital_gradient_ / abs(orbital_hessian_)
        tk_   = yk_ / abs(orbital_hessian_)
      
        tau = my_ddot(rot_pair_%n_tot,dk_,1,yk_,1)
        gamma = my_ddot(rot_pair_%n_tot,yk_,1,sk_,1) 
        beta = gamma / tau
      
        delta = my_ddot(rot_pair_%n_tot,dk_,1,gkp1_,1)
        gamma = delta * my_ddot(rot_pair_%n_tot,yk_,1,tk_,1) 
        beta = beta - theta * gamma / (tau*tau)
      
        gk_ = gkp1_
        dk_ = -sk_ + beta*dk_

      else
        !!!!! Preconditioned Kou-Dai 2015 Method; Theta=1.0
        theta = 1.0_wp
        gkp1_ = orbital_gradient_
        yk_   = gkp1_ - gk_
        sk_   = orbital_gradient_ / abs(orbital_hessian_)
        tk_   = yk_ / abs(orbital_hessian_)
        rk_   = dk_ * abs(orbital_hessian_)

        tau = my_ddot(rot_pair_%n_tot,dk_,1,yk_,1)
        gamma = my_ddot(rot_pair_%n_tot,yk_,1,sk_,1) 
        beta = gamma / tau

        delta = my_ddot(rot_pair_%n_tot,dk_,1,gkp1_,1)
        gamma = delta * my_ddot(rot_pair_%n_tot,yk_,1,tk_,1) 
        beta = beta - theta * gamma / (tau*tau)

        gamma = my_ddot(rot_pair_%n_tot,dk_,1,rk_,1)
        beta = beta - delta / gamma

        if (beta < 0.1_wp * delta / gamma) then
          beta = 0.1_wp * delta / gamma
          gamma = 0.0_wp
        else
          gamma = 0.5_wp * delta / tau
        end if 

        gk_ = gkp1_
        dk_ = -sk_ + beta*dk_ + gamma*tk_
      end if
      
    end do

    last_energy = e_total_

    orbopt_data(11) = real(iter,kind=wp)
    orbopt_data(12) = grad_norm_
    orbopt_data(13) = last_energy - initial_energy
    orbopt_data(14) = real(converged,kind=wp)

    ! deallocate indexing arrays
    call deallocate_indexing_arrays()

    ! allocate matrices for DIIS extrapolation
    call deallocate_diis_data()

    ! final deallocation
    call deallocate_final()

    if ( log_print_ == 1 ) close(fid_)

    return

  end subroutine focas_optimize

  subroutine deallocate_final()
    implicit none
    call deallocate_temporary_fock_matrices()
    if ( df_vars_%use_df_teints == 1 )       call deallocate_qint()
    if (allocated(orbital_gradient_))        deallocate(orbital_gradient_)

    if (allocated(dk_))                      deallocate(dk_)
    if (allocated(sk_))                      deallocate(sk_)
    if (allocated(tk_))                      deallocate(tk_)
    if (allocated(rk_))                      deallocate(rk_)
    if (allocated(gk_))                      deallocate(gk_)
    if (allocated(gkp1_))                    deallocate(gkp1_)
    if (allocated(yk_))                      deallocate(yk_)

    if (allocated(kappa_))                   deallocate(kappa_)
    if (allocated(rot_pair_%pair_offset))    deallocate(rot_pair_%pair_offset)
    if (allocated(df_vars_%class_to_df_map)) deallocate(df_vars_%class_to_df_map)
!    if (allocated(df_vars_%noccgempi))       deallocate(df_vars_%noccgempi)
    call deallocate_transformation_matrices()
    call deallocate_hessian_data()
    return
  end subroutine deallocate_final

  subroutine allocate_initial()
    implicit none

    ! allocate orbital gradient/kappa
    allocate(orbital_gradient_(rot_pair_%n_tot))
    allocate(kappa_(rot_pair_%n_tot))
    orbital_gradient_ = 0.0_wp
    kappa_            = 0.0_wp
 
    ! allocate nonlinear conjugate gradient vectors
    allocate(dk_(rot_pair_%n_tot))
    allocate(sk_(rot_pair_%n_tot))
    allocate(tk_(rot_pair_%n_tot))
    allocate(rk_(rot_pair_%n_tot))
    allocate(gk_(rot_pair_%n_tot))
    allocate(gkp1_(rot_pair_%n_tot))
    allocate(yk_(rot_pair_%n_tot))

    ! hessian data
    call allocate_hessian_data()

    return
  end subroutine allocate_initial

  integer function df_map_setup(nnz_int2)
    implicit none
    integer(ip), intent(in) :: nnz_int2
    integer(ip) :: num
    integer :: npair,i_sym,i_class,ic,idf,j_sym,j_class,i,j,ij_sym,j_max

    df_map_setup = 1

    ! check to make sure that input integral array is of reasonable size
    npair        = nmo_tot_* ( nmo_tot_ + 1 ) / 2

    num = nnz_int2/int(npair,kind=ip)
    num = num * int(npair,kind=ip)

    if ( num /= nnz_int2 ) then
      if ( log_print_ == 1) write(fid_,'(a)')'mod(nnz_int2,nmo_tot_*(nmo_tot_+1)/2) /= 0'
      return 
    endif
 
    ! determine the number of auxiliary function 
    df_vars_%nQ  = nnz_int2/npair

    ! stride of Q
    df_vars_%Qstride = ngem_tot_ 

    ! allocate and determine mapping array from class order to df order
    allocate(df_vars_%class_to_df_map(nmo_tot_))

    idf = 0

    do i_sym = 1 , nirrep_

      do i_class = 1 , 3

        do ic = first_index_(i_sym,i_class) , last_index_(i_sym,i_class)

          df_vars_%class_to_df_map(ic) = idf
          idf = idf + 1          

        end do

      end do

    end do

    if ( minval(df_vars_%class_to_df_map) < 0 ) then

      if ( log_print_ == 1 ) write(fid_,'(a)')'error ... min(class_to_df_map(:)) < 0 )'
      return

    end if

!    allocate(df_vars_%noccgempi(nirrep_))
!    df_vars_%noccgempi=dens_%ngempi

    df_map_setup = 0

    return
  end function df_map_setup

  function compute_approximate_de()

    implicit none

    real(wp) :: compute_approximate_de

    integer  :: i

    compute_approximate_de= 2.0_wp * my_ddot(rot_pair_%n_tot,orbital_gradient_,1,kappa_,1)

    do i = 1 , rot_pair_%n_tot

      compute_approximate_de = compute_approximate_de + kappa_(i) * orbital_hessian_(i) * kappa_(i)
 
    end do

    compute_approximate_de = 0.5_wp * compute_approximate_de 

    return

  end function compute_approximate_de

  subroutine precondition_step(step)

    implicit none

    real(wp) :: step(:)

    integer  :: i
    real(wp) :: h_val

    do i = 1 , rot_pair_%n_tot

      h_val = orbital_hessian_(i)

      if ( h_val < 0.0_wp ) then
        h_val = - h_val
        num_negative_diagonal_hessian_ = num_negative_diagonal_hessian_ + 1
      endif

      step(i) = - orbital_gradient_(i) / h_val

    end do
 
    if ( test_full_hessian_aa_ == 1 ) &
           & i = collect_gH(step,full_orbital_hessian_%HinvXg_aa,'H2g')

    return

  end subroutine precondition_step

  subroutine setup_rotation_indeces()
    implicit none
    ! subroutine to determine nonredundant orbital pairs
    integer :: i_sym,i,j,i_class,j_class,j_class_start,j_start,npair,npair_type(4),pair_ind

    ! initialize orbital pair type indeices
    rot_pair_%act_doc_type = 1
    rot_pair_%ext_doc_type = 2
    rot_pair_%act_act_type = 3
    rot_pair_%ext_act_type = 4

    ! initialize number of rotation paors per irrep
    trans_%npairpi = 0 

    ! initialize counters for each pair type
    npair_type     = 0

    ! initialize counter for total pairs
    npair          = 0

    ! allocate rotation index offset matrix
    allocate(rot_pair_%pair_offset(nirrep_,4))

    ! loop over symmetries for i

    do i_sym = 1 , nirrep_

      ! loop over orbital classes for i

      do i_class = 1 , 3

        ! determine first class for j

        j_class_start = i_class + 1
        if ( ( include_aa_rot_ == 1 ) .and. ( i_class == 2 ) ) j_class_start = i_class

        ! loop over classes for j

        do j_class = j_class_start , 3

          ! figure out the type of pair counter

          if ( i_class == 1 ) then

            pair_ind = rot_pair_%act_doc_type
  
            if ( j_class == 3 ) pair_ind = rot_pair_%ext_doc_type

          else

            pair_ind = rot_pair_%act_act_type

            if ( j_class == 3 ) pair_ind = rot_pair_%ext_act_type

          endif

          ! save offset for this rotation type

          rot_pair_%pair_offset(i_sym,pair_ind) = npair

          ! loop over i indeces

          do i = first_index_(i_sym,i_class) , last_index_(i_sym,i_class)

            ! determine first index for j
 
            j_start = first_index_(i_sym,j_class)
            if ( i_class == j_class ) j_start = i + 1

            ! loop over j indeces

            do j = j_start , last_index_(i_sym,j_class)            
  
              ! update rotation pair count

              npair                 = npair + 1 

              ! update pair count for this symmetry

              trans_%npairpi(i_sym) = trans_%npairpi(i_sym) + 1

              ! update pair counter for this pair

              npair_type(pair_ind)  = npair_type(pair_ind) + 1

            end do ! end j loop

          end do ! end i loop
        
        end do ! end j_class loop 

      end do ! end i_class loop 

    end do ! end i_sym loop

    ! save number of rotation pairs
    rot_pair_%n_tot   = sum(trans_%npairpi)

    ! save type of pair counters
    rot_pair_%n_ad    = npair_type(rot_pair_%act_doc_type)
    rot_pair_%n_ed    = npair_type(rot_pair_%ext_doc_type)
    rot_pair_%n_aa    = npair_type(rot_pair_%act_act_type)
    rot_pair_%n_ea    = npair_type(rot_pair_%ext_act_type)

    return
  end subroutine setup_rotation_indeces

  subroutine allocate_indexing_arrays(nirrep)
    implicit none
    integer, intent(in) :: nirrep
    nirrep_ = nirrep
    call allocate_indexing_array_help(dens_)
    call allocate_indexing_array_help(ints_)
    allocate(first_index_(nirrep_,3))
    allocate(last_index_(nirrep_,3))
    allocate(orb_sym_scr_(nmo_tot_))
    allocate(nfzcpi_(nirrep_))
    allocate(ndocpi_(nirrep_))
    allocate(nactpi_(nirrep_))
    allocate(nextpi_(nirrep_))
    first_index_ = 0
    last_index_  = 0
    orb_sym_scr_ = 0
    nfzcpi_      = 0
    ndocpi_      = 0
    nactpi_      = 0
    nextpi_      = 0
    return 
  end subroutine allocate_indexing_arrays

  subroutine allocate_indexing_array_help(styp)
    implicit none
    type(sym_info) :: styp
    allocate(styp%ngempi(nirrep_),styp%nnzpi(nirrep_),styp%offset(nirrep_),styp%gemind(nmo_tot_,nmo_tot_))
    styp%ngempi = 0
    styp%nnzpi  = 0
    styp%offset = 0
    styp%gemind = 0
    return
  end subroutine allocate_indexing_array_help

  subroutine deallocate_indexing_arrays()
    implicit none
    if (allocated(first_index_)) deallocate(first_index_)
    if (allocated(last_index_))  deallocate(last_index_)
    if (allocated(orb_sym_scr_)) deallocate(orb_sym_scr_)
    if (allocated(nfzcpi_))      deallocate(nfzcpi_)
    if (allocated(ndocpi_))      deallocate(ndocpi_)
    if (allocated(nactpi_))      deallocate(nactpi_)
    if (allocated(nextpi_))      deallocate(nextpi_)    
    call deallocate_indexing_arrays_help(dens_)
    call deallocate_indexing_arrays_help(ints_)
    return
  end subroutine deallocate_indexing_arrays

  subroutine deallocate_indexing_arrays_help(styp)
    implicit none
    type(sym_info) :: styp
    if ( allocated(styp%ngempi) ) deallocate(styp%ngempi)
    if ( allocated(styp%nnzpi) ) deallocate(styp%nnzpi)
    if ( allocated(styp%offset) ) deallocate(styp%offset)
    if ( allocated(styp%gemind) ) deallocate(styp%gemind)
    return
  end subroutine deallocate_indexing_arrays_help

  subroutine setup_indexing_arrays(nfzcpi,ndocpi,nactpi,nextpi)
    implicit none
    integer, intent(in) :: nfzcpi(nirrep_),ndocpi(nirrep_),nactpi(nirrep_),nextpi(nirrep_)
    integer :: irrep,oclass,i,j,i_sym,j_sym,ij_sym
    ! ** Figure out index of the first orbital in each class and each irrep
    ! doubly occupied orbitals
    first_index_(1,1)=1
    do irrep=2,nirrep_
      first_index_(irrep,1) = first_index_(irrep-1,1) + ndocpi(irrep-1)
    end do
    ! active orbitals
    first_index_(1,2) = first_index_(nirrep_,1) + ndocpi(nirrep_)
    do irrep=2,nirrep_
      first_index_(irrep,2) = first_index_(irrep-1,2) + nactpi(irrep-1)
    end do
    ! external orbitals
    first_index_(1,3) = first_index_(nirrep_,2) + nactpi(nirrep_)
    do irrep=2,nirrep_
      first_index_(irrep,3) = first_index_(irrep-1,3) + nextpi(irrep-1)
    end do
    ! ** Figure out indexof the last orbital in each class and each irrep
    ! doubly occupied orbitals
    do irrep=1,nirrep_
      last_index_(irrep,1) = first_index_(irrep,1) + ndocpi(irrep)-1
    end do
    ! active orbitals
    do irrep=1,nirrep_
      last_index_(irrep,2) = first_index_(irrep,2) + nactpi(irrep)-1
    end do
    ! active orbitals
    do irrep=1,nirrep_
      last_index_(irrep,3) = first_index_(irrep,3) + nextpi(irrep)-1
    end do
    ! determine number of orbitals per irrep for each class
    nfzcpi_ = nfzcpi
    ndocpi_ = last_index_(:,1) - first_index_(:,1) + 1
    nactpi_ = last_index_(:,2) - first_index_(:,2) + 1
    nextpi_ = last_index_(:,3) - first_index_(:,3) + 1
    ! ** save orbital symmetries so that we can set up geminal indices **
    do oclass=1,3
      do irrep=1,nirrep_
        do i=first_index_(irrep,oclass),last_index_(irrep,oclass)
          orb_sym_scr_(i) = irrep
        end do
      end do
    end do
    ! figure out reduced geminal indeces for integral addressing
    do i=1,nmo_tot_
      i_sym = orb_sym_scr_(i)
      do j=1,i
        j_sym = orb_sym_scr_(j)
        ij_sym = group_mult_tab_(i_sym,j_sym)
        ints_%ngempi(ij_sym) = ints_%ngempi(ij_sym) + 1
        ints_%gemind(i,j) = ints_%ngempi(ij_sym)
        ints_%gemind(j,i) = ints_%gemind(i,j)
      end do
    end do
    do irrep=1,nirrep_
      ints_%nnzpi(irrep) = ints_%ngempi(irrep) * ( ints_%ngempi(irrep) + 1 ) / 2
    end do
    do irrep=2,nirrep_
      ints_%offset(irrep) = ints_%offset(irrep-1) + ints_%nnzpi(irrep-1)
    end do
    ! figure out reduced geminal indeces for integral addressing
    do i=ndoc_tot_+1,ndoc_tot_+nact_tot_
      i_sym = orb_sym_scr_(i)
      do j=ndoc_tot_+1,i
        j_sym = orb_sym_scr_(j)
        ij_sym = group_mult_tab_(i_sym,j_sym)
        dens_%ngempi(ij_sym) = dens_%ngempi(ij_sym) + 1
        dens_%gemind(i,j) = dens_%ngempi(ij_sym)
        dens_%gemind(j,i) = dens_%gemind(i,j)
      end do
    end do
    do irrep=1,nirrep_
      dens_%nnzpi(irrep) = dens_%ngempi(irrep) * ( dens_%ngempi(irrep) + 1 ) / 2
    end do
    do irrep=2,nirrep_
      dens_%offset(irrep) = dens_%offset(irrep-1) + dens_%nnzpi(irrep-1)
    end do
    deallocate(orb_sym_scr_)
    return
  end subroutine setup_indexing_arrays

  subroutine print_info()
    integer :: irrep,i
    ! print the information gathered so far
    write(fid_,'(a5,2x,3(3(a4,1x),5x))')'irrep','d_f','d_l','n_d','a_f','a_l','n_a','e_f','e_l','n_e'
    do irrep=1,nirrep_
      write(fid_,'(i5,2x,3(3(i4,1x),5x))')irrep,(first_index_(irrep,i),last_index_(irrep,i),&
               & last_index_(irrep,i)-first_index_(irrep,i)+1,i=1,3)
    end do
    write(fid_,'(a)')'density information'
    write(fid_,'(a,8(i9,1x))')'ngempi(:)=',dens_%ngempi
    write(fid_,'(a,8(i9,1x))')' nnzpi(:)=',dens_%nnzpi
    write(fid_,'(a,8(i9,1x))')'offset(:)=',dens_%offset
    write(fid_,'(a)')'integral information'
    write(fid_,'(a,8(i9,1x))')'ngempi(:)=',ints_%ngempi
    if ( df_vars_%use_df_teints == 0 ) then
      write(fid_,'(a,8(i12,1x))')' nnzpi(:)=',ints_%nnzpi
      write(fid_,'(a,8(i12,1x))')'offset(:)=',ints_%offset
    else
      write(fid_,'(a,8(i12,1x))')' nnzpi(:)=',int(ints_%ngempi,kind=ip)*int(df_vars_%nQ,kind=ip)
    endif
    write(fid_,'(a)')'rotation pair information:'
    write(fid_,'(4(a,1x,i6,2x),a,1x,i9)')'act-doc pairs:',rot_pair_%n_ad,'ext-doc pairs:',rot_pair_%n_ed,&
                             &'act-act pairs:',rot_pair_%n_aa,'ext-act pairs:',rot_pair_%n_ea,&
                             &'total orbital pairs:',rot_pair_%n_tot
    write(fid_,*)
!
    return
  end subroutine print_info

end module focas_driver
