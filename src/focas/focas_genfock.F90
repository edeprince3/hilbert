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

module focas_genfock

  use focas_data
  use focas_driver
  use focas_gradient

  implicit none

  contains

    subroutine compute_genfock(den1,den2,int1,int2,den1_nnz,den2_nnz,int1_nnz,int2_nnz,    &
                                       & ndocpi,nactpi,nextpi,nirrep,orbopt_data,fname,gen_fock_out, &
                                       & fock_dim)

      implicit none

      integer, intent(in)     :: fock_dim
      integer, intent(in)     :: nirrep
      integer, intent(in)     :: den1_nnz
      integer, intent(in)     :: den2_nnz
      integer, intent(in)     :: int1_nnz
      integer(ip), intent(in) :: int2_nnz

      real(wp), intent(inout) :: gen_fock_out(fock_dim)

      real(wp), intent(in)    :: den1(den1_nnz)
      real(wp), intent(in)    :: den2(den2_nnz)
      real(wp), intent(in)    :: int1(int1_nnz)
      real(wp), intent(in)    :: int2(int2_nnz)
 
      integer, intent(in)     :: ndocpi(nirrep)
      integer, intent(in)     :: nactpi(nirrep)
      integer, intent(in)     :: nextpi(nirrep)

      integer                 :: nfzcpi(nirrep)

      real(wp), intent(inout) :: orbopt_data(14)

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

        write(fid_,'(a)')'computing general Fock matrix'

      endif

      ! set the number of frozen doubly occupied orbitals to zer for each IRREP
      nfzcpi = 0 

      ! allocate variables and set up mapping arrays
      call allocate_genfock_initial(nfzcpi,ndocpi,nactpi,nextpi,nirrep,int2_nnz)

      ! compute generalized Fock matrix 
      call build_entire_gen_fock(int1,int2,den1,den2,gen_fock_out)

      ! deallocate variables
      call deallocate_genfock_final()

      return

    end subroutine compute_genfock

    subroutine build_entire_gen_fock(int1,int2,den1,den2,gen_fock_out)

      real(wp), intent(in) :: int1(:),int2(:),den1(:),den2(:)
      real(wp) :: gen_fock_out(:)
       
      integer :: offset,error,p_class,q_class

      ! calculate inactive Fock matrix
      if ( df_vars_%use_df_teints == 0 ) then
        call compute_f_i(int1,int2)
      else
        call compute_f_i_df_coulomb(int1,int2)
        call compute_f_i_df_exchange(int2)
      endif
      call transpose_matrix(fock_i_)

      ! calculate active Fock matrix
      if ( df_vars_%use_df_teints == 0 ) then
        call compute_f_a(den1,int2)
      else
        call compute_f_a_df_coulomb(den1,int2)
        call compute_f_a_df_exchange(den1,int2)
      endif
      call transpose_matrix(fock_a_)

      ! calculate auxiliary q matrix
      if ( df_vars_%use_df_teints == 0 ) then
        call compute_q(den2,int2)
      else
        call compute_q_df(den2,int2)
      endif

      ! calculate auxiliary z matrix (contraction of fock_i and den1)
      call compute_z(den1)

      offset  = 0
    
      gen_fock_out = 0.0_wp
 
      ! compute nonzero blocks of generalized Fock matrix in the order
      ! dd / da / de / ad / aa / ae 
      do p_class = 1 , 2

        do q_class = 1 , 3

          error = gen_fock_block(p_class,q_class)

        end do

      end do

!      if ( log_print_ == 1 ) then
!
!        ! debug printing
!
!        offset = 0
!
!        do p_class = 1 , 2
!
!          do q_class = 1 , 3
!
!            error = print_gen_fock_block(p_class,q_class)
!
!          end do
!
!        end do
!
!      endif
 
      return

        contains

          integer function gen_fock_block(p_class,q_class)

            integer, intent(in) :: p_class,q_class

            integer  :: p,q
            integer  :: p_i,q_i
            integer  :: p_c,q_c
            integer  :: nmo_p,nmo_q
            integer  :: p_sym

            real(wp) :: val
 
            do p_sym = 1 , nirrep_

              if ( q_class == 1 ) nmo_q = ndocpi_(p_sym)
              if ( q_class == 2 ) nmo_q = nactpi_(p_sym)
              if ( q_class == 3 ) nmo_q = nextpi_(p_sym)

              if ( p_class == 1 ) nmo_p = ndocpi_(p_sym)
              if ( p_class == 2 ) nmo_p = nactpi_(p_sym)
              if ( p_class == 3 ) nmo_p = nextpi_(p_sym)

              if ( nmo_p * nmo_q == 0 ) cycle

              p_c = 0

              do p = first_index_(p_sym,p_class) , last_index_(p_sym,p_class)

                p_i = trans_%class_to_irrep_map(p)

                q_c = 0

                p_c = p_c + 1

                do q = first_index_(p_sym,q_class) , last_index_(p_sym,q_class)

                  q_i = trans_%class_to_irrep_map(q)

                  q_c = q_c + 1

                  if ( p_class == 1 ) then

                    ! d-d  / d-a / d-e terms --> F(p,q) = 2 ( F_i(q,p) + F_a(q,p) ) 
                    ! these terms follow the deifinition in Eq. 12.5.9 in Helgaker

                    val = 2.0_wp * ( fock_i_%occ(p_sym)%val(q_i,p_i) + fock_a_%occ(p_sym)%val(q_i,p_i) ) 

                  elseif ( p_class == 2 ) then

                    ! a-d / a-a / a-e terms --> F(p,q) = z(p,q) - q(p,q)
                    ! the matrix z is the contration of F_i and D1 
                    ! relative to Eq. 12.5.10 in Helgaker, the indexing is reversed

                    val = z_(p - ndoc_tot_ , q) + q_(p - ndoc_tot_ , q)

                  end if

                  if ( gen_fock_out( offset + ( p_c - 1 ) * nmo_q + q_c ) /= 0.0_wp ) write(*,*)'problem'

                  gen_fock_out( offset + ( p_c - 1 ) * nmo_q + q_c ) = val

                end do

              end do

              offset = offset + nmo_p * nmo_q

            end do

            gen_fock_block = 0

            return

          end function gen_fock_block

          integer function print_gen_fock_block(p_class,q_class)

            integer, intent(in) :: p_class,q_class

            integer  :: p,q
            integer  :: p_c,q_c
            integer  :: nmo_p,nmo_q
            integer  :: p_sym,offset_tmp

            real(wp) :: val

            if ( p_class == 1 ) then
              
              if ( q_class == 1 ) then
               
                write(fid_,'(a)')'F(P,Q) --> doc-doc'

              elseif ( q_class == 2 ) then

                write(fid_,'(a)')'F(P,Q) --> doc-act'

              else
 
                write(fid_,'(a)')'F(P,Q) --> doc-ext'

              end if

            else

              if ( q_class == 1 ) then
               
                write(fid_,'(a)')'F(P,Q) --> act-doc'

              elseif ( q_class == 2 ) then

                write(fid_,'(a)')'F(P,Q) --> act-act'

              else

                write(fid_,'(a)')'F(P,Q) --> act-ext'

              end if

            endif

            write(fid_,*)

            do p_sym = 1 , nirrep_

              if ( q_class == 1 ) nmo_q = ndocpi_(p_sym)
              if ( q_class == 2 ) nmo_q = nactpi_(p_sym)
              if ( q_class == 3 ) nmo_q = nextpi_(p_sym)

              if ( p_class == 1 ) nmo_p = ndocpi_(p_sym)
              if ( p_class == 2 ) nmo_p = nactpi_(p_sym)
              if ( p_class == 3 ) nmo_p = nextpi_(p_sym)

              if ( nmo_p * nmo_q == 0 ) cycle

              write(fid_,'(a,1x,i1,1x,2(a,1x,i3,1x))')'p_sym:',p_sym,'nmo_p:',nmo_p,'nmo_q:',nmo_q

              write(fid_,*)

              if ( nmo_p * nmo_q == 0 ) cycle

              p_c = 0
          
              do p = first_index_(p_sym,p_class) , last_index_(p_sym,p_class)

                q_c = 0

                p_c = p_c + 1

                do q = first_index_(p_sym,q_class) , last_index_(p_sym,q_class)

                  q_c = q_c + 1

                  write(fid_,'(2(i3,1x),es17.10)')p_c,q_c, &
                   & gen_fock_out( offset + ( p_c - 1 ) * nmo_q + q_c )

                end do

              end do

              offset = offset + nmo_p * nmo_q

              write(fid_,*)

            end do

            print_gen_fock_block = 0

            return

          end function print_gen_fock_block

    end subroutine build_entire_gen_fock

    subroutine allocate_genfock_initial(nfzcpi,ndocpi,nactpi,nextpi,nirrep,int2_nnz)

      implicit none

      integer, intent(in)     :: nirrep 
      integer(ip), intent(in) :: int2_nnz
      integer, intent(in)     :: ndocpi(nirrep)
      integer, intent(in)     :: nactpi(nirrep)
      integer, intent(in)     :: nextpi(nirrep)
      integer, intent(in)     :: nfzcpi(nirrep)

      integer :: error

      ! calculate the total number of orbitals in space
      ndoc_tot_ = sum(ndocpi)
      nact_tot_ = sum(nactpi)
      next_tot_ = sum(nextpi)
      nmo_tot_  = ndoc_tot_+nact_tot_+next_tot_

      ! zero out total number of frozen doubly occupied orbitals
      nfzc_tot_ = sum(nfzcpi) 

      ! allocate indexing arrays
      call allocate_indexing_arrays(nirrep)

      ! determine integral/density addressing arrays
      call setup_indexing_arrays(nfzcpi,ndocpi,nactpi,nextpi)

      ! allocate transformation matrices
      call allocate_transformation_matrices()

      ! determine indexing arrays (needed for sorts in the integral transformation step)
      call determine_transformation_maps()

      ! set up df mapping arrays if density-fitted 2-e integrals are used
      error = 0
      if ( df_vars_%use_df_teints == 1 ) error = df_map_setup(int2_nnz)
      if ( error /= 0 ) call abort_print(20)

      ! allocate temporary matrices (fock_a_,fock_i_,q_, and z_)
      call allocate_temporary_fock_matrices()

      ! allocate intermediate matrices for DF integrals
      if ( df_vars_%use_df_teints == 1 ) call allocate_qint()

      return

    end subroutine allocate_genfock_initial

    subroutine deallocate_genfock_final()

      implicit none 

      ! deallocate indexing arrays
      call deallocate_indexing_arrays()
       
      ! deallocate transformation matrices
      call deallocate_transformation_matrices()

      if (allocated(df_vars_%class_to_df_map)) deallocate(df_vars_%class_to_df_map)
!      if (allocated(df_vars_%noccgempi))       deallocate(df_vars_%noccgempi)

      ! deallocate temporary fock matrices 
      call deallocate_temporary_fock_matrices()

      ! allocate intermediate matrices for DF integrals
      if ( df_vars_%use_df_teints == 1 ) call deallocate_qint()

      return

    end subroutine deallocate_genfock_final

end module focas_genfock
