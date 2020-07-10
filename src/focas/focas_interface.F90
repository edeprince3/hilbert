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

subroutine focas_interface(mo_coeff_out,integrals_1,nnz_i1,integrals_2,nnz_i2,density_1,nnz_d1,&
           &density_2,nnz_d2,syms,ncore_in,nact_in,nvirt_in,nirrep_in,orbopt_data_io, &
           &orbopt_log_file,Xcc)

  
  use focas_driver, only        : focas_optimize
  use focas_semicanonical, only : compute_semicanonical_mos
  use focas_genfock, only       : compute_genfock
  use focas_full_hessian, only  : return_full_hessian

  implicit none
  integer, parameter :: fid=99
  integer, parameter :: wp = selected_real_kind(10)
  integer, parameter :: ip = selected_int_kind(16)                  ! 64-bit integers (integral addressing)
  integer, parameter :: mult_tab(8,8) = reshape( (/ &  ! irrep multiplication table
      & 1,2,3,4,5,6,7,8, &
      & 2,1,4,3,6,5,8,7, &
      & 3,4,1,2,7,8,5,6, &
      & 4,3,2,1,8,7,6,5, &
      & 5,6,7,8,1,2,3,4, &
      & 6,5,8,7,2,1,4,3, &
      & 7,8,5,6,3,4,1,2, &
      & 8,7,6,5,4,3,2,1  /), (/8,8/) )

  real(wp) :: orbopt_data_io(15)
  integer :: nirrep_in,ncore_in,nact_in,nvirt_in
  integer :: nnz_d1,nnz_d2,nnz_i1
  integer(ip) :: nnz_i2
  real(wp) :: Xcc((ncore_in+nact_in+nvirt_in)*(ncore_in+nact_in+nvirt_in))
  real(wp) :: integrals_1(nnz_i1),integrals_2(nnz_i2),density_1(nnz_d1),density_2(nnz_d2)
  real(wp) :: mo_coeff_out(ncore_in+nact_in+nvirt_in,ncore_in+nact_in+nvirt_in)
  character(120) :: orbopt_log_file
  integer  :: syms(int(orbopt_data_io(3))+ncore_in+nact_in+nvirt_in)

  integer :: nactpi(nirrep_in),ndocpi(nirrep_in),nextpi(nirrep_in),nfzcpi(nirrep_in)
  integer :: ndoc,nact,next,nmo,nfzc,nirrep
  integer :: nnz_int1,nnz_den1,nnz_den2,df_ints,Xdim
  integer(ip) :: nnz_int2
  integer :: first_index(nirrep_in,3)
  integer :: last_index(nirrep_in,3)
  integer :: nnz_den_psi4(nirrep_in)
  integer :: nnz_den_new(nirrep_in)
  integer(ip) :: nnz_int(nirrep_in)
  integer :: offset_den_psi4(nirrep_in)
  integer :: offset_den_new(nirrep_in)
  integer :: offset_irrep(nirrep_in)
  integer(ip) :: offset_int(nirrep_in)
  integer :: offset_irrep_int1(nirrep_in)
  integer :: offset_irrep_den1(nirrep_in)
  integer, allocatable :: gemind_int(:,:),gemind_act(:,:),gemind_den_new(:,:),gemind_int_new(:,:)
  integer, allocatable :: energy_to_class_map(:),energy_to_irrep_map(:)
  integer, allocatable :: class_to_energy_map(:),class_to_irrep_map(:)
  real(wp), allocatable :: mo_coeff(:,:)

  integer :: irrep

  nfzc=int(orbopt_data_io(3))
  ndoc=ncore_in+nfzc 
  nact=nact_in
  next=nvirt_in
  nmo=ndoc+nact+next
  nirrep=nirrep_in
  ! set density-fitted integral flag
  df_ints = int(orbopt_data_io(10))

  call allocate_mapping_arrays(1)

  call setup_symmetry_arrays(syms)

  call initial_sort()

  nnz_int1 = gemind_int_new(last_index(nirrep,3),last_index(nirrep,3))

  nnz_den1 = gemind_den_new(last_index(nirrep,2),last_index(nirrep,2))

  if (df_ints == 0 ) then

    nnz_int2 = sum(nnz_int)  

  else

    nnz_int2 = nnz_i2 

  endif

  nnz_den2 = sum(nnz_den_new)

  if ( int(orbopt_data_io(9)) > 0 ) then

    call focas_optimize(mo_coeff,integrals_1,nnz_int1,integrals_2,nnz_int2,               &
                      & density_1(1:nnz_den1),nnz_den1,density_2(1:nnz_den2),nnz_den2,    &
                      & nfzcpi,ndocpi,nactpi,nextpi,nirrep,orbopt_data_io,orbopt_log_file)

  elseif ( int(orbopt_data_io(9)) == -1 ) then

    Xdim=(ncore_in+nact_in+nvirt_in)*(ncore_in+nact_in+nvirt_in)

    call compute_genfock(density_1(1:nnz_den1),density_2(1:nnz_den2),integrals_1,&
                      & integrals_2,nnz_den1,nnz_den2,nnz_int1,nnz_int2,ndocpi,   &
                      & nactpi,nextpi,nirrep,orbopt_data_io,orbopt_log_file, &
                      & Xcc,Xdim)

  elseif ( int(orbopt_data_io(9)) == -2 ) then
   
    call compute_semicanonical_mos(density_1(1:nnz_den1),density_2(1:nnz_den2),integrals_1,&
                      & integrals_2,nnz_den1,nnz_den2,nnz_int1,nnz_int2,mo_coeff,ndocpi,   &
                      & nactpi,nextpi,nirrep,orbopt_data_io,orbopt_log_file)

  elseif ( int(orbopt_data_io(9)) == -3 ) then

    Xdim = 0 

    do irrep = 1 , nirrep

      Xdim = Xdim + nactpi(irrep)*(nactpi(irrep)-1)/2

    end do

    Xdim = Xdim*(Xdim + 1)/2

    call return_full_hessian(integrals_1,nnz_int1,integrals_2,nnz_int2,density_1(1:nnz_den1), &
                       & nnz_den1,density_2(1:nnz_den2),nnz_den2,nfzcpi,ndocpi,nactpi,nextpi,  &
                       & nirrep,orbopt_data_io,Xcc(1:Xdim),Xdim)

  end if

  call final_sort()

  call allocate_mapping_arrays(-1)

  contains

    subroutine final_sort()
      implicit none

      integer :: p_sym,q_sym,r_sym,s_sym,pq_sym
      integer :: pq_c,rs_c,pq_i,rs_i,pqrs_c,pqrs_i
      integer :: p_class,q_class,r_class,s_class,q_max,s_max
      integer :: p_c,q_c,r_c,s_c
      integer :: p_i,q_i,r_i,s_i
      integer :: p,q
      integer :: pq_off,max_dim
      real(wp), allocatable :: block(:)

      if ( df_ints == 0 ) then
        max_dim = max(size(integrals_1,dim=1),maxval(nnz_int))
      else
        max_dim = size(integrals_1)
      end if
      allocate(block(max_dim))

      ! copy 1-e integrals

      block=huge(1.0_wp)
      do p_sym = 1 , nirrep
        pq_off = offset_irrep_int1(p_sym)
        do p_class = 1 , 3
          do q_class = 1 , p_class
            do p_c = first_index(p_sym,p_class) , last_index(p_sym,p_class)
              q_max=last_index(p_sym,q_class)
              p_i = class_to_irrep_map(p_c)
              if ( p_class == q_class ) q_max=p_c
              do q_c = first_index(p_sym,q_class),q_max
                q_i = class_to_irrep_map(q_c)
                pq_c = gemind_int_new(p_c,q_c)
                pq_i = pq_ind(p_i,q_i)+pq_off
                block(pq_i)=integrals_1(pq_c)
              end do
            end do
          end do
        end do
      end do
      integrals_1 = block(1:size(integrals_1,dim=1))

      ! copy mo coefficient matrix
      mo_coeff_out=0.0_wp

      do p = nfzc + 1 , nmo
        p_sym=syms(p)
        p_i=energy_to_irrep_map(p) + offset_irrep(p_sym)
        do q = nfzc + 1 , nmo
          q_sym=syms(q)
          if (q_sym/=p_sym) cycle
          q_i=energy_to_irrep_map(q) + offset_irrep(q_sym)
          mo_coeff_out(p-nfzc,q-nfzc)=mo_coeff(p_i,q_i)
        end do
      end do

      ! copy 2-e integrals

      if ( df_ints == 0 ) then
        do pq_sym = 1 , nirrep
          block=huge(1.0_wp)
          pq_off=offset_int(pq_sym)
          do p_class = 1 , 3
            do q_class = 1 , p_class
              do p_sym = 1 , nirrep
                q_sym = mult_tab(pq_sym,p_sym)
                if ( ( q_sym > p_sym ) .and. ( p_class == q_class ) ) cycle
                do r_class = 1 , 3
                  do s_class = 1 , r_class
                    do r_sym = 1 , nirrep
                      s_sym = mult_tab(pq_sym,r_sym)
                      if ( ( s_sym > r_sym ) .and. ( r_class==q_class ) ) cycle
                      do p_c = first_index(p_sym,p_class),last_index(p_sym,p_class)
                        p_i = class_to_energy_map(p_c)
                        q_max = last_index(q_sym,q_class)
                        if ( ( p_class == q_class ) .and. ( p_sym == q_sym ) ) q_max = p_c
                        do q_c = first_index(q_sym,q_class),q_max
                          q_i = class_to_energy_map(q_c)
                          pq_c = gemind_int_new(p_c,q_c)
                          pq_i = gemind_int(p_i,q_i)
                          do r_c=first_index(r_sym,r_class),last_index(r_sym,r_class)
                            r_i = class_to_energy_map(r_c)
                            s_max = last_index(s_sym,s_class)
                            if ( ( s_sym == r_sym ) .and. ( s_class == r_class ) ) s_max = r_c
                            do s_c = first_index(s_sym,s_class),s_max
                              s_i = class_to_energy_map(s_c)
                              rs_c = gemind_int_new(r_c,s_c)
                              rs_i = gemind_int(r_i,s_i)
                              pqrs_c = pq_ind(pq_c,rs_c)
                              pqrs_i = pq_ind(pq_i,rs_i)
                              block(pqrs_i) = integrals_2(pqrs_c+pq_off)
                            end do
                          end do
                        end do
                      end do
                    end do
                  end do
                end do
              end do
            end do
          end do
          integrals_2(pq_off+1:pq_off+nnz_int(pq_sym))=block(1:nnz_int(pq_sym))
        end do
      endif

      deallocate(block)

      return

    end subroutine final_sort

    subroutine initial_sort()

      implicit none

      integer :: p_sym,q_sym,r_sym,s_sym,pq_sym
      integer :: pq_c,rs_c,pq_i,rs_i,pqrs_c,pqrs_i,pp_c,qq_c
      integer :: p_class,q_class,r_class,s_class,q_max,s_max
      integer :: p_c,q_c,r_c,s_c
      integer :: p_i,q_i,r_i,s_i
      integer :: pq_off,max_dim
      integer :: p,q,p_off
      real(wp), allocatable :: block(:)
      real(wp) :: tr_d2

      if ( df_ints == 0 ) then
        max_dim = max(size(integrals_1,dim=1),nnz_d1,maxval(nnz_int))
      else
        max_dim = max(size(integrals_1),nnz_d1,maxval(nnz_den_new))
      endif
      allocate(block(max_dim))

      ! copy 1-e integrals

      block=huge(1.0_wp)
      do p_sym = 1 , nirrep
        pq_off = offset_irrep_int1(p_sym)
        do p_class = 1 , 3
          do q_class = 1 , p_class
            do p_c = first_index(p_sym,p_class) , last_index(p_sym,p_class)
              q_max=last_index(p_sym,q_class)
              p_i = class_to_irrep_map(p_c)
              if ( p_class == q_class ) q_max=p_c
              do q_c = first_index(p_sym,q_class),q_max
                q_i = class_to_irrep_map(q_c)
                pq_c = gemind_int_new(p_c,q_c)
                pq_i = pq_ind(p_i,q_i)+pq_off
                block(pq_c)=integrals_1(pq_i)
              end do
            end do
          end do 
        end do
      end do
      integrals_1 = block(1:size(integrals_1,dim=1))

      ! copy 1-e density
      tr_d2=0.0_wp
      block=huge(1.0_wp)
      do p_sym = 1 , nirrep
        pq_off = offset_irrep_den1(p_sym)
        p_off  = ndocpi(p_sym)
        do p_c = first_index(p_sym,2) , last_index(p_sym,2)
          p_i = class_to_irrep_map(p_c)
          do q_c = first_index(p_sym,2),p_c-1
            q_i = class_to_irrep_map(q_c)
            pq_c = gemind_den_new(p_c,q_c)
            pq_i = pq_ind(p_i-p_off,q_i-p_off)+pq_off
            block(pq_c)=density_1(pq_i)
          end do
          pq_c = gemind_den_new(p_c,p_c)
          pq_i = pq_ind(p_i-p_off,p_i-p_off)+pq_off
          block(pq_c)=density_1(pq_i)
          tr_d2=tr_d2+block(pq_c)
        end do
      end do
      density_1 = block(1:nnz_d1)

      ! copy mo coefficient matrix
      mo_coeff=0.0_wp

      do p = 1 , nfzc

        p_sym = syms(p)
        p_i=energy_to_irrep_map(p) + offset_irrep(p_sym)
        mo_coeff(p_i,p_i) = 1.0_wp

      end do

      do p = nfzc + 1 ,nmo
        p_sym=syms(p)
        p_i=energy_to_irrep_map(p) + offset_irrep(p_sym)
        do q = nfzc + 1 , nmo
          q_sym=syms(q)
          if (q_sym/=p_sym) cycle
          q_i=energy_to_irrep_map(q) + offset_irrep(q_sym)
          mo_coeff(p_i,q_i)=mo_coeff_out(p-nfzc,q-nfzc)
        end do
      end do

      ! copy 2-e integrals
      if ( df_ints == 0 ) then
        do pq_sym = 1 , nirrep
          block=huge(1.0_wp)
          pq_off=offset_int(pq_sym)
          do p_class = 1 , 3
            do q_class = 1 , p_class
              do p_sym = 1 , nirrep
                q_sym = mult_tab(pq_sym,p_sym)
                if ( ( q_sym > p_sym ) .and. ( p_class == q_class ) ) cycle
                do r_class = 1 , 3
                  do s_class = 1 , r_class
                    do r_sym = 1 , nirrep
                      s_sym = mult_tab(pq_sym,r_sym)
                      if ( ( s_sym > r_sym ) .and. ( r_class==q_class ) ) cycle
                      do p_c = first_index(p_sym,p_class),last_index(p_sym,p_class)
                        p_i = class_to_energy_map(p_c)
                        q_max = last_index(q_sym,q_class)
                        if ( ( p_class == q_class ) .and. ( p_sym == q_sym ) ) q_max = p_c
                        do q_c = first_index(q_sym,q_class),q_max
                          q_i = class_to_energy_map(q_c)
                          pq_c = gemind_int_new(p_c,q_c)
                          pq_i = gemind_int(p_i,q_i)
                          do r_c=first_index(r_sym,r_class),last_index(r_sym,r_class)
                            r_i = class_to_energy_map(r_c)
                            s_max = last_index(s_sym,s_class)
                            if ( ( s_sym == r_sym ) .and. ( s_class == r_class ) ) s_max = r_c                        
                            do s_c = first_index(s_sym,s_class),s_max
                              s_i = class_to_energy_map(s_c)
                              rs_c = gemind_int_new(r_c,s_c)
                              rs_i = gemind_int(r_i,s_i)
                              pqrs_c = pq_ind(pq_c,rs_c)
                              pqrs_i = pq_ind(pq_i,rs_i)
                              block(pqrs_c) = integrals_2(pqrs_i+pq_off)
                            end do
                          end do
                        end do
                      end do
                    end do
                  end do
                end do
              end do
            end do
          end do
          integrals_2(pq_off+1:pq_off+nnz_int(pq_sym))=block(1:nnz_int(pq_sym))
        end do
      end if

      ! copy/scale 2-e active density
      do pq_sym = 1 , nirrep
        block=huge(1.0_wp)
        pq_off=offset_den_psi4(pq_sym)
        do p_sym = 1 , nirrep
          q_sym = mult_tab(pq_sym,p_sym)
          if ( q_sym > p_sym ) cycle
          do r_sym = 1 , nirrep
            s_sym = mult_tab(pq_sym,r_sym)
            if ( s_sym > r_sym ) cycle
            do p_c = first_index(p_sym,2),last_index(p_sym,2)
              p_i = class_to_energy_map(p_c)
              q_max = last_index(q_sym,2)
              if ( p_sym == q_sym ) q_max = p_c
              do q_c = first_index(q_sym,2),q_max
                q_i = class_to_energy_map(q_c)
                pq_c = gemind_den_new(p_c,q_c)
                pq_i = gemind_act(p_i,q_i)
                do r_c=first_index(r_sym,2),last_index(r_sym,2)
                  r_i = class_to_energy_map(r_c)
                  s_max = last_index(s_sym,2)
                  if ( s_sym == r_sym ) s_max = r_c
                  do s_c = first_index(s_sym,2),s_max
                    s_i = class_to_energy_map(s_c)
                    rs_c = gemind_den_new(r_c,s_c)
                    rs_i = gemind_act(r_i,s_i)
                    pqrs_c = pq_ind(pq_c,rs_c)
                    pqrs_i = pq_ind(pq_i,rs_i)
                    block(pqrs_c) = den_fac(p_c,q_c,r_c,s_c) * density_2(pqrs_i+pq_off)
                  end do
                end do
              end do
            end do
          end do
        end do
        pq_off = offset_den_new(pq_sym)
        density_2(pq_off+1:pq_off+nnz_den_new(pq_sym)) = block(1:nnz_den_new(pq_sym))
      end do

      density_2(1:sum(nnz_den_new)) = 2.0_wp * density_2(1:sum(nnz_den_new))

      deallocate(block)

      return

    end subroutine initial_sort

    subroutine setup_symmetry_arrays(syms)
      implicit none
      integer :: syms(:)
      integer :: p,q,p_class,p_sym,pq_sym,p_c
      integer :: dims(nirrep),sym_class(nmo)

      ! determine the number of mos per irrep for each class

      nfzcpi = 0
      nactpi = 0
      ndocpi = 0
      nextpi = 0
  
      do p = 1 , nfzc
        nfzcpi(syms(p)) = nfzcpi(syms(p)) + 1
      end do
      do p = 1 , ndoc
        ndocpi(syms(p)) = ndocpi(syms(p)) + 1
      end do
      do p = ndoc + 1 , ndoc + nact
        nactpi(syms(p)) = nactpi(syms(p)) + 1
      end do
      do p = ndoc + nact + 1 , nmo
        nextpi(syms(p)) = nextpi(syms(p)) + 1
      end do
      ! addressing for the density/integrals is the same for psi4
      ! the offset/nnz arrays for new order are also the same

      dims=0
      do p = 1 , nmo
        p_sym = syms(p)
        do q = 1 , p
          pq_sym = mult_tab(syms(q),p_sym)
          dims(pq_sym)=dims(pq_sym)+1
          gemind_int(p,q)=dims(pq_sym)
          gemind_int(q,p)=dims(pq_sym)
        end do
      end do

      nnz_int = 0
      do pq_sym = 1 , nirrep
        nnz_int(pq_sym) = dims(pq_sym)*(dims(pq_sym)+1)/2
      end do

      offset_int = 0
      do pq_sym = 2 , nirrep
        offset_int(pq_sym) = nnz_int(pq_sym-1) + offset_int(pq_sym-1)
      end do

      ! determine energy --> irrep map
      dims=0
      do p=1,nmo
        p_sym = syms(p)
        dims(p_sym) = dims(p_sym) + 1
        energy_to_irrep_map(p) = dims(p_sym)
      end do
 
      offset_irrep = 0
      do p_sym = 2 , nirrep
        offset_irrep(p_sym) = offset_irrep(p_sym-1) + dims(p_sym-1)
      end do

      ! determine dimensions for active density arrays

      dims=0
      do p = ndoc+1 , ndoc+nact
        p_sym = syms(p)
        do q = ndoc+1 , p
          pq_sym = mult_tab(syms(q),p_sym)
          dims(pq_sym)=dims(pq_sym)+1
          gemind_act(p,q) = dims(pq_sym)
          gemind_act(q,p) = dims(pq_sym)
        end do
      end do

      nnz_den_psi4 = 0
      do pq_sym = 1 , nirrep
        nnz_den_psi4(pq_sym) = dims(pq_sym)*(dims(pq_sym)+1)/2
      end do
      offset_den_psi4 = 0
      do pq_sym = 2 , nirrep
        offset_den_psi4(pq_sym) = nnz_den_psi4(pq_sym-1) + offset_den_psi4(pq_sym-1)
      end do

      ! determine irrep->class & class->irrep map arrays

      p_class=0
      do p_sym = 1 , nirrep
        do p = 1 , ndoc
          if ( syms(p) /= p_sym ) cycle
          p_class = p_class + 1
          energy_to_class_map(p) = p_class
          class_to_energy_map(p_class) = p
        end do    
      end do
      do p_sym = 1 , nirrep
        do p = ndoc + 1 , ndoc + nact
          if ( syms(p) /= p_sym ) cycle
          p_class = p_class + 1
          energy_to_class_map(p) = p_class
          class_to_energy_map(p_class) = p
        end do
      end do
      do p_sym = 1 , nirrep
        do p = ndoc + nact + 1 , nmo
          if ( syms(p) /= p_sym ) cycle
          p_class = p_class + 1
          energy_to_class_map(p) = p_class
          class_to_energy_map(p_class) = p
        end do
      end do

      ! set up orbital symmetry arrays in new order

      do p = 1 , nmo
        p_c = energy_to_class_map(p)
        sym_class(p_c) = syms(p)
      end do

      ! set up geminal addressing arrays for the integrals

      dims=0
      do p = 1 , nmo
        p_sym = sym_class(p)
        do q = 1 , p
          pq_sym = mult_tab(sym_class(q),p_sym)
          dims(pq_sym)=dims(pq_sym)+1
          gemind_int_new(p,q)=dims(pq_sym)
          gemind_int_new(q,p)=dims(pq_sym)
        end do
      end do

      ! set up geminal addressing arrays for the densities (active only)

      dims=0
      gemind_den_new = 0
      do p = ndoc + 1 , ndoc + nact
        p_sym = sym_class(p)
        do q = ndoc + 1 , p
          pq_sym = mult_tab(sym_class(q),p_sym)
          dims(pq_sym)=dims(pq_sym)+1
          gemind_den_new(p,q)=dims(pq_sym)
          gemind_den_new(q,p)=dims(pq_sym)
        end do
      end do

      nnz_den_new = 0
      do pq_sym = 1 , nirrep
        nnz_den_new(pq_sym) = dims(pq_sym)*(dims(pq_sym)+1)/2
      end do
      offset_den_new = 0
      do pq_sym = 2 , nirrep
        offset_den_new(pq_sym) = nnz_den_new(pq_sym-1) + offset_den_new(pq_sym-1)
      end do

      first_index = 0
      last_index = 0

      first_index(1,1)=1
      last_index(1,1)=first_index(1,1) + ndocpi(1) - 1
      do p_sym = 2,nirrep
        first_index(p_sym,1) = first_index(p_sym-1,1) + ndocpi(p_sym-1)
        last_index(p_sym,1) = first_index(p_sym,1) + ndocpi(p_sym) - 1
      end do

      first_index(1,2)=first_index(nirrep,1)+ndocpi(nirrep)
      last_index(1,2)=first_index(1,2) + nactpi(1) - 1
      do p_sym = 2,nirrep
        first_index(p_sym,2) = first_index(p_sym-1,2) + nactpi(p_sym-1)
        last_index(p_sym,2) = first_index(p_sym,2) + nactpi(p_sym) - 1
      end do

      first_index(1,3)=first_index(nirrep,2)+nactpi(nirrep)
      last_index(1,3)=first_index(1,3) + nextpi(1) - 1
      do p_sym = 2,nirrep
        first_index(p_sym,3) = first_index(p_sym-1,3) + nextpi(p_sym-1)
        last_index(p_sym,3) = first_index(p_sym,3) + nextpi(p_sym) - 1
      end do

      dims=0
      do p_sym = 1 , nirrep
        dims(p_sym) = 0
        do p_class = 1 , 3
          do p = first_index(p_sym,p_class),last_index(p_sym,p_class)
            dims(p_sym) = dims(p_sym) + 1
            class_to_irrep_map(p) = dims(p_sym)
          end do
        end do
      end do
      offset_irrep_int1 = 0
      do p_sym = 2 , nirrep
        offset_irrep_int1(p_sym) = offset_irrep_int1(p_sym-1) + dims(p_sym-1)*(dims(p_sym-1)+1)/2
      end do

      dims=0
      do p_sym = 1 , nirrep
        do p_class = 2 , 2
          do p = first_index(p_sym,p_class),last_index(p_sym,p_class)
            dims(p_sym) = dims(p_sym) + 1
          end do
        end do
      end do
      offset_irrep_den1 = 0
      do p_sym = 2 , nirrep
        offset_irrep_den1(p_sym) = offset_irrep_den1(p_sym-1) + dims(p_sym-1)*(dims(p_sym-1)+1)/2
      end do
      return
    end subroutine setup_symmetry_arrays

    pure function pq_ind(i,j)
! this function computes the two-electron index index (lower triangular reference)
! index = ii*(ii-1)/2+jj where ii=max(i,j) and jj=min(i,j)
! the ishft(k,-1) divides the value of the integer k by 2 and seems to be somewhat
! faster than the regular human-readable expression
      implicit none
      integer, intent(in) ::i,j
      integer :: pq_ind
      if (i.ge.j) then
        pq_ind=ishft(i*(i-1),-1)+j
        return
      else
        pq_ind=ishft(j*(j-1),-1)+i
        return
      end if
    end function pq_ind

    function den_fac(i,j,k,l)
      integer, intent(in) :: i,j,k,l
      real(wp) :: den_fac
      integer :: fac
      fac=1
      if ( i /= j ) fac = 2
      if ( k /= l ) fac = 2 * fac
      if (pq_ind(i,j) /= pq_ind(k,l)) fac = 2 * fac
      den_fac = 1.0_wp / real(fac,kind=wp)
    end function den_fac

    subroutine allocate_mapping_arrays(alloc)

      integer, intent(in) :: alloc

      if ( alloc > 0 ) then

        allocate(gemind_int(nmo,nmo))
        allocate(gemind_act(nmo,nmo))
        allocate(gemind_int_new(nmo,nmo))
        allocate(gemind_den_new(nmo,nmo))
        allocate(energy_to_irrep_map(nmo))
        allocate(energy_to_class_map(nmo))
        allocate(class_to_energy_map(nmo))
        allocate(class_to_irrep_map(nmo))
        allocate(mo_coeff(nmo,nmo))

      elseif ( alloc < 0 ) then

        deallocate(gemind_int)
        deallocate(gemind_act)
        deallocate(gemind_int_new)
        deallocate(gemind_den_new)
        deallocate(energy_to_irrep_map)
        deallocate(energy_to_class_map)
        deallocate(class_to_energy_map)
        deallocate(class_to_irrep_map)
        deallocate(mo_coeff)
       
      end if

    end subroutine allocate_mapping_arrays

end subroutine focas_interface
