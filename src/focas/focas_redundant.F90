module focas_redundant

  use focas_data

  implicit none

  real(wp), parameter :: on_small_tol = 0.01_wp
  real(wp), parameter :: on_large_tol = 2.0_wp -on_small_tol

  contains

    subroutine compute_opdm_nos(opdm)

      implicit none

      real(wp), intent(in) :: opdm(:)

      real(wp), allocatable :: opdm_block(:,:),nos(:)

      integer :: block_sym,nmo,error

      do block_sym = 1 , nirrep_

        nmo = nactpi_(block_sym)

        if ( nmo == 0 ) cycle

        allocate(opdm_block(nmo,nmo),nos(nmo))

        ! gather the opdm elements for this block

        call gather_opdm_block(opdm,opdm_block,block_sym)

        error = diagonalize_opdm_block(nos,opdm_block,nmo)

        call analyze_block_nos(nos,opdm_block,nmo,block_sym)

        if (error /= 0 ) call abort_print(40)

        deallocate(opdm_block,nos)

      end do

      return

    end subroutine compute_opdm_nos

    subroutine analyze_block_nos(nos,opdm_block,nmo,block_sym)

      implicit none

      integer, intent(in)  :: nmo,block_sym
      real(wp), intent(in) :: nos(nmo),opdm_block(nmo,nmo) 

      integer :: n_small,n_large
      integer :: m

      n_small = 0
      n_large = 0

      do m = 1 , nmo

        if ( nos(m) >= on_large_tol ) n_large = n_large + 1
        if ( nos(m) <= on_small_tol ) n_small = n_small + 1

      end do

      if ( log_print_ == 1 ) then

        if ( ( n_large /= 0 ) .and. ( ndocpi_(block_sym) /= 0 ) ) then
        
          write(fid_,'(a,i2,1x,a,1x,es10.3,1x,a,1x,i2)')'warning, there are',n_large,&
               & 'large ONs above',on_large_tol,'for irrep:',block_sym  
          write(fid_,'(10(f6.4,1x))')nos(nmo-n_large+1:nmo)
        end if

        if ( ( n_small /= 0 ) .and. ( nextpi_(block_sym) /= 0 ) ) then

          write(fid_,'(a,i2,1x,a,1x,es10.3,1x,a,1x,i2)')'warning, there are',n_small,&
               & 'small ONs below',on_small_tol,'for irrep:',block_sym
          write(fid_,'(10(f6.4,1x))')nos(1:n_small)

        end if

      end if

      return

    end subroutine analyze_block_nos

    integer function diagonalize_opdm_block(ons,opdm_block,nmo)

      implicit none

      integer  :: nmo
      real(wp) :: ons(nmo),opdm_block(nmo,nmo)

      integer  :: il,liwork,lwork
      real(wp) :: vl,vu,diag_tol

      real(wp), allocatable :: work(:),vecs(:,:)
      integer, allocatable :: isuppz(:),iwork(:)
      integer :: iwork_tmp(1),neig_found
      real(wp) :: work_tmp(1,1),val
 
!      ! determine the values for # of eigenvlaues to be foudn as well as min/max
!      il=1
!      vl=-huge(1.0_wp)
!      vu=0.0_wp
!      diag_tol=2.0_wp*epsilon(1.0_wp)
!
!      ! allocate the first temporary matrix
!      allocate(isuppz(2*nmo))
!
!      ! figure out optimal dimensions for iwork and work
!      call dsyevr('v','a','u',nmo,opdm_block,nmo,vl,vu,il,nmo,diag_tol,neig_found,ons,vecs,&
!                 & nmo,isuppz,work_tmp,-1,iwork_tmp,-1,diagonalize_opdm_block)
!
!      ! something did not go right so return without diagonalization
!      if ( diagonalize_opdm_block /= 0 ) then
!        deallocate(isuppz)
!        return
!      end if
!
!      ! save optimal workspace values
!      liwork=iwork_tmp(1)
!      lwork=int(work_tmp(1,1))
!
!      ! allocate remaining temporary matrices
!      allocate(work(lwork),iwork(liwork),vecs(nmo,nmo))
! 
!      ! diagonalize
!      call dsyevr('v','a','u',nmo,opdm_block,nmo,vl,vu,il,nmo,diag_tol,neig_found,ons,vecs,&
!                 & nmo,isuppz,work,lwork,iwork,liwork,diagonalize_opdm_block)
!
!      ! save eigenvectors
!      opdm_block = vecs
!
!      ! deallocate temporary matrices
!      deallocate(isuppz,work,iwork,vecs)

      call dsyev('v','u',nmo,opdm_block,nmo,ons,work_tmp,-1,diagonalize_opdm_block)
      lwork=int(work_tmp(1,1))
      if ( diagonalize_opdm_block /= 0 ) then
        return
      endif
      allocate(work(lwork))
      call dsyev('v','u',nmo,opdm_block,nmo,ons,work,lwork,diagonalize_opdm_block)
      if ( diagonalize_opdm_block /= 0 ) then 
        deallocate(work)
        return
      endif
    
      deallocate(work)
      return

    end function diagonalize_opdm_block

    subroutine gather_opdm_block(opdm,opdm_block,block_sym)

      implicit none
 
      real(wp), intent(in) :: opdm(:)
      integer, intent(in)  :: block_sym
      real(wp)             :: opdm_block(:,:)
     
      integer :: mn,m,n,m_i,n_i

      do m = first_index_(block_sym,2) , last_index_(block_sym,2)

        m_i = trans_%class_to_irrep_map(m) - ndocpi_(block_sym)

        do n = first_index_(block_sym,2) , m

          n_i = trans_%class_to_irrep_map(n) - ndocpi_(block_sym)

          mn = dens_%gemind(m,n)

          opdm_block(n_i,m_i) = opdm(mn)
          opdm_block(m_i,n_i) = opdm_block(n_i,m_i)

        end do

      end do

      return

    end subroutine gather_opdm_block

end module focas_redundant
