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

module focas_transform_oeints

  use focas_data

  implicit none

  contains

    integer function transform_mocoeff(mo_coeff)

      implicit none
      real(wp) :: mo_coeff(:,:)
      integer :: i_sym,nmo_L,i_offset,max_nmopi,first_i,last_i
      real(wp), allocatable :: block_tmp(:,:)

      transform_mocoeff = 1

      if ( size(mo_coeff,dim = 1) /= nmo_tot_ ) then
        if ( log_print_ == 1 ) then
          write(fid_,'(a)')'dimension of mo_coeff(:,:) and nmo_tot_ do not match'
        end if
        return
      end if

      max_nmopi = maxval(trans_%nmopi)
   
      allocate(block_tmp(max_nmopi,max_nmopi))

      do i_sym = 1 , nirrep_

        nmo_L    = trans_%nmopi(i_sym)

        if ( nmo_L == 0 ) cycle

        i_offset = trans_%offset(i_sym)

        first_i  = i_offset + 1
        last_i   = i_offset + nmo_L


        call dgemm('N','N',nmo_L,nmo_L,nmo_L,1.0_wp,mo_coeff(first_i:last_i,first_i:last_i),nmo_L,&
                     trans_%u_irrep_block(i_sym)%val,nmo_L,0.0_wp,block_tmp,max_nmopi)

        mo_coeff(first_i:last_i,first_i:last_i) = block_tmp(1:nmo_L,1:nmo_L)

      end do

      deallocate(block_tmp)

      transform_mocoeff = 0

      return

    end function transform_mocoeff

    integer function transform_oeints(int1)
      implicit none
  
      real(wp) :: int1(:)
      real(wp), allocatable :: irrep_block_1(:,:),irrep_block_2(:,:)
 
      integer :: max_nmopi,ij_class,i_offset,i_class,j_class,i_irrep,j_irrep,sym_L
      integer :: nfz_L,nfz_R,nac_L,nac_R,nmo_L,nmo_R

      ! initialize return value

      transform_oeints = 1

      ! determine maximum block size and allocate temporary matrix
 
      max_nmopi = maxval(trans_%nmopi)
 
      allocate(irrep_block_1(max_nmopi,max_nmopi),irrep_block_2(max_nmopi,max_nmopi))

      ! loop over irreps for i

      do sym_L = 1 , nirrep_

        ! initialize temporary matrix

        irrep_block_1 = 0.0_wp

        ! number of orbitals in this block as well as offset for indexing

        nmo_L    = trans_%nmopi(sym_L)
        
        ! number of frozen core orbitals

        nfz_L = nfzcpi_(sym_L)
 
        ! no orbitals with this symmetry or all orbitals frozen

        if ( ( nmo_L == 0 ) .or. ( trans_%U_eq_I(sym_L) == 1 ) ) cycle

        i_offset = trans_%offset(sym_L)

        ! **************
        ! *** GATHER ***
        ! ************** 

        ! loop over symmetry-reduced indeces for i 

        do i_irrep = 1 , nmo_L

          ! class-reduced index for i needed for integral addressing

          i_class = trans_%irrep_to_class_map(i_irrep+i_offset)

          ! loop over symmetry-reduced indeces for j

          do j_irrep = 1 , i_irrep

            ! class reduced index for j needed for integral addressing

            j_class  = trans_%irrep_to_class_map(j_irrep+i_offset)

            ! one-electron index

            ij_class = ints_%gemind(i_class,j_class) 

            ! save integral making sure that the current irrep block is symmetric

            irrep_block_1(j_irrep,i_irrep) = int1(ij_class)
            irrep_block_1(i_irrep,j_irrep) = irrep_block_1(j_irrep,i_irrep)

          end do ! end j_irrep loop \  0  L^T /u

        end do ! end i_irrep loop

        ! *****************
        ! *** TRANSFORM ***
        ! *****************

! ********************************************************
! THIS CODE DOES NOT INCORPORATE THE SPARSE STRUCTURE OF U
! ********************************************************

        call dgemm('N','N',nmo_L,nmo_L,nmo_L,1.0_wp,irrep_block_1,max_nmopi,&
                     trans_%u_irrep_block(sym_L)%val,nmo_L,0.0_wp,irrep_block_2,max_nmopi)

        call dgemm('T','N',nmo_L,nmo_L,nmo_L,1.0_wp,trans_%u_irrep_block(sym_L)%val,nmo_L, &
                     irrep_block_2,max_nmopi,0.0_wp,irrep_block_1,max_nmopi)

!
!
! ***********************************************************
! THIS CODE INCORPORATES SPARSE STRUCTURE DUE TO FROZEN CORES
! ***********************************************************
!        if ( nfz_L == 0 ) then
!
!          ! no frozen orbitals so transform entire block
!
!          call dgemm('N','N',nmo_L,nmo_L,nmo_L,1.0_wp,irrep_block_1,max_nmopi,&
!                       trans_%u_irrep_block(sym_L)%val,nmo_L,0.0_wp,irrep_block_2,max_nmopi)
!
!          call dgemm('T','N',nmo_L,nmo_L,nmo_L,1.0_wp,trans_%u_irrep_block(sym_L)%val,nmo_L, &
!                       irrep_block_2,max_nmopi,0.0_wp,irrep_block_1,max_nmopi)
!
!        else
!
!          ! some frozen orbitals
!          ! 
!          !               / I_L  0  \   / h11  h12 \   / I_R  0  \   / I_L  0  \   / h11 * I_R  h12 * R \
!          ! L^T * h * R = |         | x |          | x |         | = |         | x |                    | 
!          !               \  0  L^T /   \ h21  h22 /   \  0   R  /   \  0  L^T /   \ h21 * I_R  h22 * R /
!          !
!          !               / I_L * h11 * I_R  I_L * h12 * R \   /     h11         h12 * R    \
!          !             = |                                | = |                            |
!          !               \ L^T * h21 * I_R  L^T * h22 * R /   \  L^T * h21   L^T * h22 * R /
!          !
!          ! matrix block dimensions
!          !
!          ! L   ==> nac_L x nac_L          R   ==> nac_R x nac_R
!          ! h11 ==> nfz_L x nfz_R          h12 ==> nfz_L x nac_R
!          ! h21 ==> nac_L x nfz_R          h22 ==> nac_L x nac_R
!
!          nfz_R = nfz_L
!          nmo_R = nmo_L
!          nac_L = nmo_L - nfz_L
!          nac_R = nmo_R - nfz_R        
!
!          ! h12 block ==> C = h12 * R
!
!          call dgemm('n','n',nfz_L,nac_R,nac_R,1.0_wp, &
!                  &  irrep_block_1(1:nfz_L,nfz_R+1:nmo_R),nfz_L,                         &
!                  &  trans_%u_irrep_block(sym_L)%val(nfz_R+1:nmo_R,nfz_R+1:nmo_R),nac_R, &
!                  &  0.0_wp,irrep_block_2(1:nfz_L,nfz_R+1:nmo_R),nfz_L) 
!
!          ! h21 block == C = L^T * h21 (this is not really needed by include for consistency)
!
!          call dgemm('t','n',nac_L,nfz_R,nac_L,1.0_wp,                                   &
!                   & trans_%u_irrep_block(sym_L)%val(nfz_L+1:nmo_L,nfz_L+1:nmo_L),nac_L, & 
!                   & irrep_block_1(nfz_L+1:nmo_L,1:nfz_R),nac_L,                         &
!                   & 0.0_wp,irrep_block_2(nfz_L+1:nmo_L,1:nfz_L),nac_L) 
!
!          ! h22 block ==> L^T * h22 * R (two half transforms)
!
!          call dgemm('n','n',nac_L,nac_R,nac_R,1.0_wp,                                  &
!                  &  irrep_block_1(nfz_L+1:nmo_L,nfz_R+1:nmo_R),nac_L,                  &
!                  &  trans_%u_irrep_block(sym_L)%val(nfz_R+1:nmo_R,nfz_R+1:nmo_R),nac_R,&
!                  &  0.0_wp,irrep_block_2(nfz_L+1:nmo_L,nfz_R+1:nmo_R),nac_L)
!
!          call dgemm('t','n',nac_L,nac_R,nac_L,1.0_wp,                                   &
!                   & trans_%u_irrep_block(sym_L)%val(nfz_L+1:nmo_L,nfz_L+1:nmo_L),nac_L, &
!                   & irrep_block_2(nfz_L+1:nmo_L,nfz_R+1:nmo_R),nac_L,                   &
!                   & 0.0_wp,irrep_block_1(nfz_L+1:nmo_L,nfz_R+1:nmo_R),nac_L)
!
!
!          irrep_block_1(nfz_L+1:nmo_L,1:nfz_L) = irrep_block_2(nfz_L+1:nmo_L,1:nfz_L) 
!          irrep_block_1(1:nfz_R,nfz_R+1:nmo_R) = irrep_block_2(1:nfz_R,nfz_R+1:nmo_R)
!              
!        end if
!
!
!
        ! ***************
        ! *** SCATTER ***
        ! ***************

        do i_irrep = 1 , nmo_L

          ! class-reduced index for i needed for integral addressing

          i_class = trans_%irrep_to_class_map(i_irrep+i_offset)

          ! loop over symmetry-reduced indeces for j

          do j_irrep = 1, i_irrep

            ! class reduced index for j needed for integral addressing

            j_class  = trans_%irrep_to_class_map(j_irrep+i_offset)

            ! one-electron index

            ij_class = ints_%gemind(i_class,j_class)

            ! save integral

            int1(ij_class) = irrep_block_1(j_irrep,i_irrep)
 
          end do ! end j_irrep loop

        end do ! end i_irrep loop

      end do ! end sym_L loop

      ! deallocate temporary matrices

      deallocate(irrep_block_1,irrep_block_2)

      ! no errors encountered

      transform_oeints = 0

      return
  
    end function transform_oeints

end module focas_transform_oeints
