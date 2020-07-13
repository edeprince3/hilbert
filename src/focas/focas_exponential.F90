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

module focas_exponential

  use focas_data

  implicit none

  real(wp), parameter :: max_error_tolerance = 1.0e-10_wp

  contains

    subroutine compute_exponential(kappa_in)
      implicit none
      ! subroutine to compute matrix exponential of a skew-symmetric real matrix K
      ! U = exp(K)
      ! Since the matrix K is block-diagonal, U is as well. Thus, the matrix
      ! exponential is computed for each block separately
      
      integer :: i_sym,max_nmopi,error
      real(wp), intent(in) :: kappa_in(:)
      real(wp), allocatable :: k_block(:,:)

      ! maximum size of temporary matrix

      max_nmopi = maxval(trans_%nmopi)
 
      ! allocate temporary matrix

      allocate(k_block(max_nmopi,max_nmopi))

      do i_sym=1,nirrep_

        error = gather_kappa_block(kappa_in,k_block,i_sym)
        if (error /= 0) call abort_print(10) 

        error = 0
        if ( trans_%nmopi(i_sym) > 0 ) error = compute_block_exponential(k_block,i_sym,max_nmopi)
        if (error /= 0) call abort_print(11)
         
      end do

      deallocate(k_block)

      return
 
    end subroutine compute_exponential

    integer function compute_block_exponential(K,block_sym,max_dim)
      implicit none
      ! function to compute the matrix exponential of a matrix according to
      ! U = exp(K) = X * cos(d) * X^T + K * X * d^(-1) * sin(d) * X^T
      ! where X and d are solutions of the eigenvalue equation K^2 * X = lambda * X
      ! and d = sqrt(-lambda)
      integer, intent(in)  :: block_sym,max_dim
      real(wp), intent(in) :: K(max_dim,max_dim)

      real(wp), allocatable :: K2(:,:),X(:,:),tmp_mat_1(:,:),tmp_mat_2(:,:)
      real(wp), allocatable :: d(:),work(:)
      integer, allocatable :: isuppz(:),iwork(:)
      integer :: iwork_tmp(1)
      real(wp) :: work_tmp(1,1),val

      integer :: block_dim,i,j, il,iu,neig_found,lwork,liwork,success,nfzc,nmo
      real(wp) :: vl,vu,diag_tol,max_normalization_error,max_orthogonality_error

      ! initialize return value
      compute_block_exponential = 1

      ! check to see if this block is going to be equal to the I

      if ( trans_%U_eq_I(block_sym) == 1 ) then

        trans_%u_irrep_block(block_sym)%val = 0.0_wp

        do i = 1 , trans_%nmopi(block_sym)

          trans_%u_irrep_block(block_sym)%val(i,i) = 1.0_wp

        end do

        compute_block_exponential = 0

        return

      end if

      ! number of frozen doubly occupied orbitals
     
      nfzc = nfzcpi_(block_sym)

      ! save total number of orbitals

      nmo = trans_%nmopi(block_sym)

      ! dimension of active block

      block_dim = nmo - nfzc

      ! ***************************
      ! allocate temporary matrices
      ! ***************************

      allocate(K2(block_dim,block_dim))
      allocate(d(block_dim))
      allocate(X(block_dim,block_dim))

      ! ***************
      ! compute and K^2
      ! ***************

      do i = 1 , block_dim
        call my_dcopy(block_dim,K(nfzc+1:block_dim+nfzc,i+nfzc),1,X(:,i),1)
      end do

      call dgemm('n','n',block_dim,block_dim,block_dim,1.0_wp,K(nfzc+1:nmo,nfzc+1:nmo),& 
      & block_dim,X,block_dim,0.0_wp,K2,block_dim)

      ! ***************
      ! diagonalize K^2
      ! ***************
      
      X = K2
      call dsyev('v','u',block_dim,X,block_dim,d,work_tmp,-1,success)
      lwork=int(work_tmp(1,1))
      if ( success /= 0 ) then 
        deallocate(K2,d,X)
        return
      end if
      allocate(work(lwork))
      call dsyev('v','u',block_dim,X,block_dim,d, work,lwork,success)
      if ( success /= 0 ) then 
        deallocate(work)
        return
      endif

      ! *************************
      ! initialize unitary matrix
      ! *************************

      trans_%u_irrep_block(block_sym)%val = 0.0_wp

      do i = 1 , nfzc

        trans_%u_irrep_block(block_sym)%val(i,i) = 1.0_wp

      end do

      ! ***********************************************************************************
      ! compute exponential U = exp(K) = X * cos(d) * X^(T) + K * X * d^(-1) * sin(d) * X^T
      ! first compute the terms involving sin(d) followed by the terms involving cos(d)
      ! ***********************************************************************************

      ! scale eigenvalues
      do i = 1 , block_dim
        if (d(i) < 0.0_wp ) then 
          d(i) = sqrt(-d(i))
        else
          d(i) = sqrt(d(i))
        end if
      end do

      allocate(tmp_mat_1(block_dim,block_dim))
      allocate(tmp_mat_2(block_dim,block_dim))

      ! compute d^(-1) * sin(d) * X^T
      ! since both d and sin(d) are diagonal matrices, d^(-1) * sin(d) is also diagonal with diagonal elements sin(d(i))/d(i)
      ! since d^(-1) * sin(d) is diagonal, we can perform the matrix product C = D * M efficiently by recognizing that the ith
      ! row of C is just a scaled version of the ith row of M ... C(:,i) = D(i,i) * M(:,i) 
      ! note to self :: below, we are storing the transpose of d^(-1) * sin(d) * X^T since we directly copy rows --> rows

      do i = 1 , block_dim
        val = 1.0_wp
        if ( d(i) /= 0.0_wp ) val = sin(d(i)) / d(i) 

        call my_dcopy(block_dim,X(:,i),1,tmp_mat_1(:,i),1)
        call my_dscal(block_dim,val,tmp_mat_1(:,i),1)

      end do

      call dgemm('n','t',block_dim,block_dim,block_dim,1.0_wp,X,block_dim,tmp_mat_1, &
                & block_dim,0.0_wp,tmp_mat_2,block_dim)

      call dgemm('n','n',block_dim,block_dim,block_dim,1.0_wp,K(nfzc+1:nmo,nfzc+1:nmo),block_dim,tmp_mat_2, &
                & block_dim,0.0_wp,trans_%u_irrep_block(block_sym)%val(nfzc+1:nmo,nfzc+1:nmo),block_dim)
     
      ! ****************************
      ! *** COMPUTE X * cos(d) * X^T
      ! ****************************

      ! tmp_mat_1 =  cos(d) * X

      do i = 1 , block_dim
        val = cos(d(i))

        call my_dcopy(block_dim,X(:,i),1,tmp_mat_1(:,i),1)
        call my_dscal(block_dim,val,tmp_mat_1(:,i),1)        

      end do 

      call dgemm('n','t',block_dim,block_dim,block_dim,1.0_wp,X,block_dim,tmp_mat_1, & 
                & block_dim,1.0_wp,trans_%u_irrep_block(block_sym)%val(nfzc+1:nmo,nfzc+1:nmo),block_dim)

      max_orthogonality_error = 0.0_wp
      max_normalization_error = 0.0_wp

      do i = 1 , nmo

        ! compute norm of vector

        val = my_ddot(nmo,trans_%u_irrep_block(block_sym)%val(:,i),1,&
              trans_%u_irrep_block(block_sym)%val(:,i),1)

        if ( abs( 1.0_wp - val ) > max_normalization_error ) max_normalization_error = abs ( 1.0_wp - val )

        do j = 1 , i - 1

          val = my_ddot(nmo,trans_%u_irrep_block(block_sym)%val(:,i),1,&
              trans_%u_irrep_block(block_sym)%val(:,j),1)

          if ( abs( val ) > max_orthogonality_error ) max_orthogonality_error = abs ( val )
           
        end do      
 
      end do

      if ( log_print_ == 1 ) then
        if ( ( max_normalization_error > max_error_tolerance ) .or. &
           & ( max_orthogonality_error > max_error_tolerance ) ) then
          write(fid_,'(a,1x,i1,5x,a,1x,i3,5x,a,1x,es10.3,5x,a,1x,es10.3)')'irrep:',block_sym,'nmo:',block_dim,&
               & 'max(normalization_error):',max_normalization_error,'max(orthogonality_error):',max_orthogonality_error
        endif 
      endif

      ! ******************************
      ! deallocate tempporary matrices
      ! ******************************
     
      deallocate(tmp_mat_1,tmp_mat_2)
      deallocate(K2,d,X)

      compute_block_exponential = 0

      return
    end function compute_block_exponential

    integer function gather_kappa_block(kappa_in,block,block_sym)
      implicit none
      real(wp) :: kappa_in(:)
      real(wp) :: block(:,:)
      integer, intent(in) :: block_sym
      integer :: i,j,n_ij,ij,i_irrep,j_irrep,i_class,j_class,j_class_start,j_start

      gather_kappa_block = 1

      ! number of orbital pairs in this block
      n_ij = trans_%npairpi(block_sym)

      if ( n_ij == 0 ) then
        gather_kappa_block = 0
        return
      end if

      ! figure out first/last orbital pair index for this block 

      ij = 0
      if ( block_sym > 1 ) ij = sum(trans_%npairpi(1:block_sym-1))
      
      ! initialize matrix block
      block = 0.0_wp

      ! loop over i/j pairs in this block

      ! the loop structure below cycles through the possible rotation pairs
      ! rotation pairs are sorted according to symmetry and for each irrep,
      ! the rotation pairs are sorted according to orbital classes: ad,ed,aa,ea
      ! for each pair j>i; for more details, see subroutine setup_rotation_indeces in focas_main.F90

 
      do i_class = 1 , 3

        j_class_start = i_class + 1

        if ( ( include_aa_rot_ == 1 ) .and. ( i_class == 2 ) ) j_class_start = i_class

        do j_class = j_class_start , 3

          do i = first_index_(block_sym,i_class) , last_index_(block_sym,i_class)

            j_start = first_index_(block_sym,j_class)

            if ( i_class == j_class ) j_start = i + 1

            do j = j_start , last_index_(block_sym,j_class)

              ! figure out symmetry_reduced indeces

              i_irrep                = trans_%class_to_irrep_map(i)
              j_irrep                = trans_%class_to_irrep_map(j)

              ! update gradient index

              ij                     = ij + 1

              ! copy into matrix

              block(j_irrep,i_irrep) =  kappa_in(ij)
              block(i_irrep,j_irrep) = -kappa_in(ij)

            end do

          end do

        end do    

      end do ! end ij loop

      gather_kappa_block = 0

    end function gather_kappa_block

end module focas_exponential
