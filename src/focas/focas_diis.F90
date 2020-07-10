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

module focas_diis

  use focas_data

  implicit none

  contains

    subroutine diis_extrapolate(P,dP)

     implicit none

     real(wp) :: P(:),dP(:)

     integer :: i

     ! update index for current DIIS vector

     if (diis_%current_index == -1 )                 diis_%current_index = 0
     if (diis_%current_index == diis_%max_num_diis ) diis_%current_index = 0

     diis_%current_index = diis_%current_index + 1

     ! decide whether we have enough information to do the update
     
     if ( ( diis_%update == 0 ) .and. ( diis_%current_index == diis_%max_num_diis ) ) diis_%update = 1      

     ! copy current dP/P into the list of stored dP/P

     call my_dcopy(rot_pair_%n_tot,dP,1,diis_%dP(:,diis_%current_index),1)
     call my_dcopy(rot_pair_%n_tot,P,1,diis_%P(:,diis_%current_index),1)

     if ( diis_%update == 0 ) return

     ! calculate the B matrix

     call dgemm('t','n',diis_%max_num_diis,diis_%max_num_diis,rot_pair_%n_tot,&
              & 1.0_wp,diis_%dP,rot_pair_%n_tot,diis_%dP,rot_pair_%n_tot,     &
              & 0.0_wp,diis_%B(1:diis_%max_num_diis,1:diis_%max_num_diis),diis_%max_num_diis)

     diis_%B(diis_%max_num_diis+1,:)                  = -1.0_wp
     diis_%B(:,diis_%max_num_diis+1)                  = -1.0_wp
     diis_%B(diis_%max_num_diis,diis_%max_num_diis+1) =  0.0_wp

     ! set up constraint matrix

     diis_%c                                          = 0.0_wp
     diis_%c(diis_%max_num_diis+1)                    = -1.0_wp

     ! solve system of linear equations A * x = c
     ! on return, c contains the solution x

     call dgesv(diis_%max_num_diis+1,1,diis_%B,diis_%max_num_diis+1,diis_%ip,diis_%c,diis_%max_num_diis+1,diis_%error)     
 
     ! copy solution vector into new kappa

     dP = 0.0_wp
     do i = 1 , diis_%max_num_diis
       call my_daxpy(rot_pair_%n_tot,diis_%c(i),diis_%P(:,i),1,dP,1)
     end do    

     write(*,*)'new dP'
     dP = P + dP
     write(*,*)dP

!
!     call my_dcopy(rot_pair_%n_tot,dP,1,diis_%dP(:,diis_%current_index),1)
!
     return

    end subroutine diis_extrapolate

    subroutine allocate_diis_data()

      implicit none

      diis_%do_diis = 1
      if ( diis_%max_num_diis == 0 ) diis_%do_diis = 0

      if ( diis_%do_diis == 0 ) return

      if (allocated(diis_%B))  deallocate(diis_%B)
      if (allocated(diis_%c))  deallocate(diis_%c)
      if (allocated(diis_%ip)) deallocate(diis_%ip)
      if (allocated(diis_%P))  deallocate(diis_%P)
      if (allocated(diis_%dP)) deallocate(diis_%dP)

      allocate(diis_%B(diis_%max_num_diis+1,diis_%max_num_diis+1))
      allocate(diis_%c(diis_%max_num_diis+1))
      allocate(diis_%ip(diis_%max_num_diis+1))
      allocate(diis_%P(rot_pair_%n_tot,diis_%max_num_diis))
      allocate(diis_%dP(rot_pair_%n_tot,diis_%max_num_diis))
  
      diis_%B                 = 0.0_wp
      diis_%c                 = 0.0_wp 
      diis_%ip                = 0 
      diis_%dP                = 0.0_wp
      diis_%P                 = 0.0_wp
      diis_%current_index     = -1
      diis_%update            = 0

      return

    end subroutine allocate_diis_data    

    subroutine deallocate_diis_data()

      implicit none

      if ( diis_%do_diis == 0 ) return

      if (allocated(diis_%B))  deallocate(diis_%B)
      if (allocated(diis_%c))  deallocate(diis_%c)
      if (allocated(diis_%ip)) deallocate(diis_%ip)
      if (allocated(diis_%P))  deallocate(diis_%P)
      if (allocated(diis_%dP)) deallocate(diis_%dP)      

      return

    end subroutine deallocate_diis_data
 
end module focas_diis
