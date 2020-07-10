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

module focas_energy
  use focas_data

  implicit none

  contains

  function compute_energy_tindex(f_i,f_a,q,z,h,d1)

    ! function to compute energy using only two-index quantities: i
    ! F_i, F_a, Q, Z, D, and h
    !
    ! E = 0.5 * [ sum_{p,q} d1(p,q)h(p|q) + sum_{p} F(p|p) ]
    ! 
    ! terms proportional to d1
    !
    ! E = sum_{i} h(i|i) + sum_{u>v} d1(u|v) * h(u|v) + 0.5 * sum_{u} d1(u|u) * h(u|u)  
    !
    ! terms proportional to generalized Fock matrix
    !
    ! E = 0.5 * [ sum_{i} F(i|i) + sum_{u} F(u|u) ]
    !   = sum_{i} ( f_i(i|i) + f_a(i|i) ) + 0.5 * sum_{u} ( z(u|u) + q(u|u) )
    
    implicit none

    type(matrix_block), intent(in) :: f_i(:)
    type(matrix_block), intent(in) :: f_a(:)
    real(wp), intent(in)           :: q(:,:)
    real(wp), intent(in)           :: z(:,:)
    real(wp), intent(in)           :: h(:)
    real(wp), intent(in)           :: d1(:)
    
    integer :: i,u,v,i_i,u_i
    integer :: i_sym,u_sym
    integer :: ii_int,uv_int,uv_den,uu_int,uu_den
   
    real(wp) :: energy,compute_energy_tindex
 
    energy = 0.0_wp

    ! terms involving inactive orbitals
 
    do i_sym = 1 , nirrep_

      do i = first_index_(i_sym,1) , last_index_(i_sym,1)

        ii_int = ints_%gemind(i,i)
        i_i    = trans_%class_to_irrep_map(i)

        energy = energy + h(ii_int) + f_i(i_sym)%val(i_i,i_i) + f_a(i_sym)%val(i_i,i_i)

      end do

    end do 

    do u_sym = 1 , nirrep_

      do u = first_index_(u_sym,2) , last_index_(u_sym,2)

        do v = first_index_(u_sym,2) , u - 1

          uv_den = dens_%gemind(u,v)
          uv_int = ints_%gemind(u,v)

          energy = energy + d1(uv_den) * h(uv_int)

        end do

        uu_den = dens_%gemind(u,u)
        uu_int = ints_%gemind(u,u)

        energy = energy + 0.5_wp * ( d1(uu_den) * h(uu_int) + &
                                   & z(u-ndoc_tot_,u) + q(u-ndoc_tot_,u) ) 
 
      end do

    end do

    compute_energy_tindex=energy

    return

  end function compute_energy_tindex

  subroutine compute_energy(int1,int2,den1,den2)
    implicit none
    real(wp), intent(in) :: int1(:),int2(:),den1(:),den2(:)

    ! initialize values
    e1_c_          = 0.0_wp
    e1_a_          = 0.0_wp
    e2_cc_         = 0.0_wp
    e2_ca_         = 0.0_wp
    e2_aa_         = 0.0_wp   
    e_frozen_core_ = 0.0_wp

    ! compute core contribution to 1-e energy ; 2 indeces \in D
    call compute_core_1e(int1,e1_c_)

    ! compute active contribution to 1-e energy ; 2 indeces \in A
    call compute_active_1e(int1,den1,e1_a_)

    ! compute core-core contribution to 2-e energy ; 4 indeces \in D
    if ( df_vars_%use_df_teints == 0 ) then
      call compute_core_core_2e(int2,e2_cc_)
    else
      call compute_core_core_2e_df(int2,e2_cc_)
    end if

    ! compute core-active contribution to 2-e energy ; 2 indeces \in D && 2 indeces \in A
    if ( df_vars_%use_df_teints == 0 ) then
      call compute_core_active_2e(int2,den1,e2_ca_) 
    else
      call compute_core_active_2e_df(int2,den1,e2_ca_)
    end if

    ! compute active-active contribution to 2-e energy ; 4 indeces \in A
    if ( df_vars_%use_df_teints == 0 ) then
      call compute_active_active_2e(int2,den2,e2_aa_)
    else
      call compute_active_active_2e_df(int2,den2,e2_aa_)
    end if

    ! frozen core energy
    e_frozen_core_ = e1_c_ + e2_cc_

    ! total energy
    e_total_       = e1_c_ + e1_a_ + e2_cc_ + e2_ca_ + e2_aa_

    ! active energy
    e_active_      = e_total_ - e_frozen_core_

    ! total 1-e energy
    e1_total_      = e1_c_ + e1_a_

    ! total 2-e energy
    e2_total_      = e2_cc_ + e2_ca_ + e2_aa_

!    write(fid_,*)e1_c_,e2_cc_
!
!    write(fid_,*)e_frozen_core_,e_active_,e_total_
!    write(fid_,*)e1_total_,e2_total_
!
!    stop

    return
  end subroutine compute_energy

  subroutine compute_core_1e(int1,e_out)
    implicit none
    ! subroutine to accumulate 1-e contribution to the core energy given by
    ! e1_core = 2.0 * SUM(i \in D) h(i|i)
    real(wp), intent(in) :: int1(:)
    real(wp) :: e_out
    integer :: i,i_sym
    integer :: ii

    ! initialize energy

    e_out = 0.0_wp

    ! loop over irreps for i

    do i_sym = 1 , nirrep_

      ! loop over i indeces ( i \in D )

      do i = first_index_(i_sym,1) , last_index_(i_sym,1)

        ! geminal index

        ii = ints_%gemind(i,i)

        ! update core energy

        e_out = e_out + int1(ii)

      end do ! end i loop

    end do ! end i_sym loop

    ! take into account double occupancy of orbital

    e_out = 2.0_wp * e_out

    return
  end subroutine compute_core_1e

  subroutine compute_active_1e(int1,den1,e_out)
    implicit none
    ! subroutine to accumulate the 1-e contribution to the active energy given by
    ! e1_active = sum(i,j \in A) { h(i|j) * d1(i|j) }
    ! because d1 and h are symmetric, we can perform a restricted (i<=j) summation 
    ! and multiply the off-diagonal elements by 2
    real(wp), intent(in) :: int1(:),den1(:)
    real(wp) :: e_out
    integer :: i,j,i_sym
    integer :: ij_den,ij_int
    real(wp) :: e_scr

    ! initialize 1-e active energy

    e_out = 0.0_wp

    ! loop over irreps for i (not that j_sym == i_sym) 

    do i_sym = 1, nirrep_

      ! loop over i indeces ( i \in A )

      do i = first_index_(i_sym,2) , last_index_(i_sym,2)

        ! diagonal element
        ij_int = ints_%gemind(i,i)
        ij_den = dens_%gemind(i,i)
        e_out  = e_out + den1(ij_den) * int1(ij_int)

        ! initialize temporary e value
        e_scr = 0.0_wp

        do j = i + 1 , last_index_(i_sym,2)

          ! off-diagonal element
          ij_int = ints_%gemind(i,j)
          ij_den = dens_%gemind(i,j)
          e_scr = e_scr + den1(ij_den) * int1(ij_int)

        end do ! end j loop

        ! update active 1-e energy (factor of 2 comes from permutational symmetry)

        e_out = e_out + 2.0_wp * e_scr

      end do ! end i loop

    end do ! end i_sym loop
    return
  end subroutine compute_active_1e

  subroutine compute_core_core_2e_df(int2,e_out)
    implicit none
    ! subroutine to accumulate 2-e contribution to the core energy given by
    ! e2_core = sum(i,j \in D) { 2g(ii|jj) - g(ij|ij) }
    ! integrals are stored in density-fitted 3-index format
    ! NOTE: the coulomb terms g(ii|jj) belong to the totally symmetric irrep
    real(wp), intent(in) :: int2(:)
    real(wp) :: e_out,coulomb,exchange
    integer :: j_sym,i_sym,i,j,idf,jdf
    integer(ip) :: ii,jj,ij,nQ
    real(wp) :: int_val
    real(wp) :: v_ii(df_vars_%nQ)

    nQ = int(df_vars_%nQ,kind=ip)

    ! initialize coulomb and exchange energies

    coulomb  = 0.0_wp
    exchange = 0.0_wp

    ! **************************
    ! *** i_sym > j_sym && i > j
    ! **************************

    ! loop over irreps for i

    do i_sym = 1 , nirrep_

      ! loop over irreps for j

      do j_sym = 1 , i_sym - 1

        ! loop over i indeces

        do i = first_index_(i_sym,1) , last_index_(i_sym,1)

          ! orbital index in df order
          idf = df_vars_%class_to_df_map(i)

          ! ii-geminal index
          ii  =  df_pq_index(idf,idf) 
          call my_dcopy(df_vars_%nQ,int2(ii+1:),df_vars_%Qstride,v_ii,1)

          ! loop over j indeces

          do j = first_index_(j_sym,1) , last_index_(j_sym,1)

            ! orbital index in df order
            jdf       = df_vars_%class_to_df_map(j)

            ! geminal indeces
            jj        =  df_pq_index(jdf,jdf)
            ij        =  df_pq_index(idf,jdf)

            ! calculate Coulomb contribution // g(ii,jj) ... i,j \in C)
            ! calculate exchange contribution // g(ij,ij) ... i,j \in C)

            coulomb  = coulomb  + my_ddot(df_vars_%nQ,v_ii,1,int2(jj+1:),df_vars_%Qstride)
            exchange = exchange + my_ddot(df_vars_%nQ,int2(ij+1:),df_vars_%Qstride,int2(ij+1:),df_vars_%Qstride)

          end do ! end j loop

        end do ! end i loop

      end do ! end j_sym loop

      ! ***************************
      ! *** i_sym == j_sym && i > j
      ! ***************************

      ! loop over i indeces

      do i = first_index_(i_sym,1) , last_index_(i_sym,1)

        ! orbital index in df order
        idf = df_vars_%class_to_df_map(i)

        ! ii-geminal index
        ii  = df_pq_index(idf,idf)
        call my_dcopy(df_vars_%nQ,int2(ii+1:),df_vars_%Qstride,v_ii,1)

        ! loop over j indeces

        do j = first_index_(i_sym,1) , i - 1

          ! orbital index in df order
          jdf       = df_vars_%class_to_df_map(j)

          ! geminal indeces
          jj        = df_pq_index(jdf,jdf)
          ij        = df_pq_index(idf,jdf)

          ! calculate Coulomb contribution // g(ii,jj) ... i,j \in C)
          ! calculate exchange contribution // g(ij,ij) ... i,j \in C)

          coulomb  = coulomb  + my_ddot(df_vars_%nQ,v_ii,1,int2(jj+1:),df_vars_%Qstride)
          exchange = exchange + my_ddot(df_vars_%nQ,int2(ij+1:),df_vars_%Qstride,int2(ij+1:),df_vars_%Qstride)

        end do ! end j loop

        ! ****************************
        ! *** i_sym == j_sym && i == j
        ! ****************************

        int_val  = 0.5_wp * my_ddot(df_vars_%nQ,int2(ii+1:),df_vars_%Qstride,int2(ii+1:),df_vars_%Qstride)

        coulomb  = coulomb  + int_val
        exchange = exchange + int_val    

      end do ! end i loop

    end do ! end i_sym loop

    ! subtract the exchange terms from twice the coulomb terms
    e_out = 2.0_wp * ( 2.0_wp * coulomb - exchange )

    return
  end subroutine compute_core_core_2e_df

  subroutine compute_core_core_2e(int2,e_out)
    implicit none
    ! subroutine to accumulate 2-e contribution to the core energy given by
    ! e2_core = sum(i,j \in D) { 2g(ii|jj) - g(ij|ij) }
    ! integrals are stored in full 4-index format
    ! NOTE: the coulomb terms g(ii|jj) belong to the totally symmetric irrep
    real(wp), intent(in) :: int2(:)
    real(wp) :: e_out,coulomb,exchange,ccc
    integer :: j_sym,i_sym,ij_sym,i,j
    integer :: ii,jj,ij
    integer :: int_ind,ij_offset 

    ! initialize coulomb and exchange energies

    coulomb  = 0.0_wp
    exchange = 0.0_wp

    ! loop over irreps for i

    do i_sym = 1 , nirrep_

      ! loop over irreps for j

      do j_sym = 1 , nirrep_

        ! loop over i indeces

        do i = first_index_(i_sym,1) , last_index_(i_sym,1)
 
          ! ii-geminal index
          ii = ints_%gemind(i,i)
  
          ! loop over j indeces

          do j = first_index_(j_sym,1) , last_index_(j_sym,1)             
 
            ! jj geminal index
            jj        = ints_%gemind(j,j)
            ! ij geminal index
            ij        = ints_%gemind(i,j)
            ! ij symmetry
            ij_sym    = group_mult_tab_(i_sym,j_sym)
            ! ij geminal offset
            ij_offset = ints_%offset(ij_sym) 
            ! Coulomb contribution // g(ii,jj) ... i,j \in C)
            int_ind   = pq_index(ii,jj)
            coulomb   = coulomb + int2(int_ind)
            ccc = int2(int_ind)
            ! exchange contribution // g(ij|ij) ... i,j \in C 
            int_ind   = pq_index(ij,ij) + ij_offset
            exchange  = exchange + int2(int_ind)

          end do ! end j loop

        end do ! end i loop

      end do ! end j_sym loop

    end do ! end i_sym loop

    ! subtract the exchange terms from twice the coulomb terms
    e_out = 2.0_wp * coulomb - exchange

    return
  end subroutine compute_core_core_2e

  subroutine compute_core_active_2e_df(int2,den1,e_out)
    implicit none
    ! subroutine to compute the core-active contribution to the 2-e energy
    ! e = 0.5 * sum(ij \in A & k \in D) { 2 * d(i|j) * [ 2g(ij|kk) - g(ik|jk) ] }
    real(wp), intent(in) :: int2(:),den1(:)
    real(wp) :: e_out
    integer :: i,j,k,i_sym,k_sym,ik_sym
    integer :: idf,jdf,kdf,ij_den,ii_den
    integer(ip) :: ii,ij,kk,ik,jk,nQ
    real(wp) :: coulomb,exchange
    real(wp) :: v_ij(df_vars_%nQ),v_ii(df_vars_%nQ)

    nQ = int(df_vars_%nQ,kind=ip)

    ! initialize energy value
    e_out = 0.0_wp

    ! loop over irreps for i

    do i_sym = 1 , nirrep_

      ! loop over irreps for k

      do k_sym = 1 , nirrep_

        ! ik geminal symmetry
        ik_sym    = group_mult_tab_(i_sym,k_sym)

        ! loop over i indeces

        do i = first_index_(i_sym,2) , last_index_(i_sym,2)

          ! orbital index in df order
          idf = df_vars_%class_to_df_map(i)

          ! ***************************
          ! *** i_sym == j_sym && i > j
          ! *************************** 

          ! loop over j indeces 

          do j = first_index_(i_sym,2) , i - 1

            ! orbital index in df order
            jdf = df_vars_%class_to_df_map(j)

            ! ij geminal integral index
            ij  = df_pq_index(idf,jdf)
            call my_dcopy(df_vars_%nQ,int2(ij+1:),df_vars_%Qstride,v_ij,1)

            ! initialize coulomb and exhange contributions
            coulomb  = 0.0_wp
            exchange = 0.0_wp

            ! loop over k indeces

            do k = first_index_(k_sym,1) , last_index_(k_sym,1)

              ! orbital index in df order
              kdf      = df_vars_%class_to_df_map(k)
 
              ! geminal indeces
              kk       = df_pq_index(kdf,kdf)  
              ik       = df_pq_index(idf,kdf) 
              jk       = df_pq_index(jdf,kdf)

              ! calculate Coulomb contribution // g(ij,kk) ... i,j \in A && k \in D
              ! calculate exchange contribution // g(ik,jk) ... i,j \in A && k \in D

              coulomb  = coulomb  + my_ddot(df_vars_%nQ,v_ij,1,int2(kk+1:),df_vars_%Qstride)
              exchange = exchange + my_ddot(df_vars_%nQ,int2(ik+1:),df_vars_%Qstride,int2(jk+1:),df_vars_%Qstride)

            end do ! end k loop

            ! ij geminal density index
            ij_den = dens_%gemind(i,j)
            ! add up terms
            e_out = e_out + den1(ij_den) * ( 2.0_wp * coulomb - exchange )

          end do ! end j loop

          ! ***************************
          ! *** i_sym == j_sym && i > j
          ! ***************************

          ! ij geminal integral index
          ii  = df_pq_index(idf,idf) 
          call my_dcopy(df_vars_%nQ,int2(ii+1:),df_vars_%Qstride,v_ii,1)

          ! initialize coulomb and exhange contributions
          coulomb  = 0.0_wp
          exchange = 0.0_wp

          ! loop over k indeces

          do k = first_index_(k_sym,1) , last_index_(k_sym,1)

            ! orbital index in df order
            kdf      = df_vars_%class_to_df_map(k)

            ! geminal indeces
            kk       = df_pq_index(kdf,kdf) 
            ik       = df_pq_index(idf,kdf) 

            ! calculate Coulomb contribution // g(ii,kk) ... i \in A && k \in D
            ! calculate exchange contribution // g(ik,ik) ... i \in A && k \in D

            coulomb  = coulomb  + my_ddot(df_vars_%nQ,v_ii,1,int2(kk+1:),df_vars_%Qstride)
            exchange = exchange + my_ddot(df_vars_%nQ,int2(ik+1:),df_vars_%Qstride,int2(ik+1:),df_vars_%Qstride)

          end do ! end k loop

          ! ij geminal density index
          ii_den = dens_%gemind(i,i)
          ! add up terms
          e_out = e_out + den1(ii_den) * ( coulomb - 0.5_wp * exchange )          

        end do ! end i loop

      end do ! end k_sym loop 

    end do ! end i_sym loop

   ! account for the j>i elements

    e_out = 2.0_wp * e_out
 
    return
  end subroutine compute_core_active_2e_df

  subroutine compute_core_active_2e(int2,den1,e_out)
    implicit none
    ! subroutine to compute the core-active contribution to the 2-e energy
    ! e = 0.5 * sum(ij \in A & k \in D) { 2 * d(i|j) * [ 2g(ij|kk) - g(ik|jk) ] }
    real(wp), intent(in) :: int2(:),den1(:)
    real(wp) :: e_out
    integer :: i,j,k,i_sym,k_sym,ik_sym
    integer :: ij_int,ij_den,kk,ik,jk
    integer :: ik_offset,ijkk,ikjk
    real(wp) :: coulomb,exchange

    ! initialize energy value
    e_out = 0.0_wp

    ! loop over irreps for i
 
    do i_sym = 1 , nirrep_

      ! loop over irreps for k

      do k_sym = 1 , nirrep_

        ! ik geminal symmetry
        ik_sym    = group_mult_tab_(i_sym,k_sym)
        ! ik geminal offset
        ik_offset = ints_%offset(ik_sym) 

        ! loop over i indeces

        do i = first_index_(i_sym,2) , last_index_(i_sym,2)

          ! loop over j indeces 

          do j = first_index_(i_sym,2) , last_index_(i_sym,2)

            ! ij geminal integral index
            ij_int = ints_%gemind(i,j)

            ! initialize coulomb and exhange contributions
            coulomb  = 0.0_wp
            exchange = 0.0_wp

            ! loop over k indeces

            do k = first_index_(k_sym,1) , last_index_(k_sym,1)

              ! kk geminal index
              kk       = ints_%gemind(k,k)
              ! Coulomb contribution
              ijkk     = pq_index(ij_int,kk)
              coulomb  = coulomb + int2(ijkk)
 
              ! ik/jk geminal indeces
              ik       = ints_%gemind(i,k)
              jk       = ints_%gemind(j,k)
              ! exchange contribution
              ikjk     = pq_index(ik,jk) + ik_offset
              exchange = exchange + int2(ikjk)

            end do ! end k loop

            ! ij geminal density index
            ij_den = dens_%gemind(i,j)
            ! add up terms
            e_out = e_out + den1(ij_den) * ( 2.0_wp * coulomb - exchange )

          end do ! end j loop

        end do ! end i loop
        
      end do ! end k_sym loop 
 
    end do ! end i_sym loop

    return
  end subroutine compute_core_active_2e

  subroutine compute_active_active_2e_df(int2,den2,e_out)
    implicit none
    ! subroutine to compute active-active contribution to 2-electron energy
    ! e = 0.5 * sum(ijkl \in A) { d2(ij|kl) * g(ij|kl) }
    real(wp), intent(in) :: int2(:),den2(:)
    real(wp) :: e_out
    integer :: i,j,k,l,i_sym,j_sym,k_sym,l_sym,ij_sym,idf,jdf,kdf,ldf
    integer :: kl_den,ij_den,ii_den,kk_den
    integer(ip) :: ij,kl,kk,ii,nQ
    integer :: den_ind,ij_den_offset
    real(wp) :: int_val
    real(wp) :: v_ii(df_vars_%nQ),v_ij(df_vars_%nQ)

    nQ = int(df_vars_%nQ,kind=ip)

    ! initialize energy
    e_out = 0.0_wp

    do i_sym = 1 , nirrep_

      do j_sym = 1 , i_sym - 1

        ij_sym = group_mult_tab_(i_sym,j_sym)

        ij_den_offset = dens_%offset(ij_sym)

        do k_sym = 1 , nirrep_
         
          l_sym = group_mult_tab_(ij_sym,k_sym)

          if ( l_sym > k_sym ) cycle

          ! *******************************************************************
          ! *** i_sym > j_sym && k_sym > l_sym --> i > j && k > l (factor of 4)
          ! ******************************************************************* 

          do i = first_index_(i_sym,2) , last_index_(i_sym,2)
 
            idf = df_vars_%class_to_df_map(i)

            do j = first_index_(j_sym,2) , last_index_(j_sym,2)

              jdf    = df_vars_%class_to_df_map(j)

              ij_den = dens_%gemind(i,j)

              ij     = df_pq_index(idf,jdf) 

              do k = first_index_(k_sym,2) , last_index_(k_sym,2)

                kdf = df_vars_%class_to_df_map(k)

                do l = first_index_(l_sym,2) , last_index_(l_sym,2)

                  ldf     = df_vars_%class_to_df_map(l)

                  kl_den  = dens_%gemind(k,l)

                  den_ind = pq_index(ij_den,kl_den) + ij_den_offset

                  kl      = df_pq_index(kdf,ldf) 

                  ! do work here i > j && k > l --> factor of 4
                  
                  int_val  = my_ddot(df_vars_%nQ,int2(ij+1:),df_vars_%Qstride,int2(kl+1:),df_vars_%Qstride)

                  e_out    = e_out + 4.0_wp * int_val * den2(den_ind)  

                end do ! end l loop

              end do ! end k loop

            end do ! end j loop

          end do ! end i loop
 
        end do ! end k_sym loop

      end do ! end j_sym loop

      ! ************************************
      ! *** i_sym == j_sym && k_sym == l_sym
      ! ************************************

      do k_sym = 1 , nirrep_

        do i = first_index_(i_sym,2) , last_index_(i_sym,2)

          idf = df_vars_%class_to_df_map(i)
 
          do j = first_index_(i_sym,2) , i -1

            jdf = df_vars_%class_to_df_map(j) 

            ij_den = dens_%gemind(i,j)

            ij     = df_pq_index(idf,jdf) 
            call my_dcopy(df_vars_%nQ,int2(ij+1:),df_vars_%Qstride,v_ij,1)

            do k = first_index_(k_sym,2) , last_index_(k_sym,2)

              kdf = df_vars_%class_to_df_map(k)

              do l = first_index_(k_sym,2) , k - 1

                ldf     = df_vars_%class_to_df_map(l)

                kl_den  = dens_%gemind(k,l)

                den_ind = pq_index(ij_den,kl_den)

                kl      = df_pq_index(kdf,ldf)  
 
                ! do work here i > j && k > l --> factor of 4

                int_val  = my_ddot(df_vars_%nQ,v_ij,1,int2(kl+1:),df_vars_%Qstride)

                e_out    = e_out + 4.0_wp * int_val * den2(den_ind)                

              end do ! end l loop

              kk_den  = dens_%gemind(k,k)
 
              den_ind = pq_index(ij_den,kk_den)

              kk      = df_pq_index(kdf,kdf) 
 
              ! do work here i > j && k == l --> factor of 2

              int_val  = my_ddot(df_vars_%nQ,v_ij,1,int2(kk+1:),df_vars_%Qstride)

              e_out    = e_out + 2.0_wp * int_val * den2(den_ind)

            end do ! end k loop

          end do ! end j loop

          ii     = df_pq_index(idf,idf) 
          call my_dcopy(df_vars_%nQ,int2(ii+1:),df_vars_%Qstride,v_ii,1)

          ii_den = dens_%gemind(i,i) 
 
          do k = first_index_(k_sym,2) , last_index_(k_sym,2)

            kdf = df_vars_%class_to_df_map(k)
 
            do l = first_index_(k_sym,2) , k - 1

              ldf = df_vars_%class_to_df_map(l)

              kl_den  = dens_%gemind(k,l) 

              den_ind = pq_index(ii_den,kl_den)
 
              kl      = df_pq_index(kdf,ldf)  

              ! do work here i == j && k > l --> factor of 2

              int_val  = my_ddot(df_vars_%nQ,v_ii,1,int2(kl+1:),df_vars_%Qstride)

              e_out    = e_out + 2.0_wp * int_val * den2(den_ind)

            end do ! end l loop

            kk_den  = dens_%gemind(k,k)
 
            den_ind = pq_index(ii_den,kk_den)

            kk      = df_pq_index(kdf,kdf) 

            ! do work here i == j && k == l --> factor of 1

            int_val  = my_ddot(df_vars_%nQ,v_ii,1,int2(kk+1:),df_vars_%Qstride)

            e_out    = e_out + int_val * den2(den_ind)

          end do ! end k loop          

        end do ! end i loop

      end do ! end k_sym loop

    end do ! end i_sym loop

    e_out = 0.5_wp * e_out

    return

  end subroutine compute_active_active_2e_df

  subroutine compute_active_active_2e(int2,den2,e_out)
    implicit none
    ! subroutine to compute active-active contribution to 2-electron energy
    ! e = 0.5 * sum(ijkl \in A) { d2(ij|kl) * g(ij|kl) }
    real(wp), intent(in) :: int2(:),den2(:) 
    real(wp) :: e_out
    integer :: i,j,k,l,i_sym,j_sym,k_sym,l_sym,ij_sym
    integer :: ij_int,kl_int,ij_den,kl_den
    integer :: ij_den_offset,ij_int_offset,den_ind,int_ind
     
    ! initialize energy
    e_out = 0.0_wp

    ! loop over irreps for i
    
    do i_sym = 1 , nirrep_

      ! loop over irreps for j
     
      do j_sym = 1 , nirrep_

        ! ij geminal symmetry    
        ij_sym        = group_mult_tab_(i_sym,j_sym)
        ! ij geminal integral offset
        ij_int_offset = ints_%offset(ij_sym)
        ! ij geminal density offset
        ij_den_offset = dens_%offset(ij_sym)

        ! loop over irreps for k

        do k_sym = 1 , nirrep_

          ! figure out symmetry of l such that ij_sym == kl_sym
          l_sym = group_mult_tab_(ij_sym,k_sym)

          ! loop over i indeces
 
          do i = first_index_(i_sym,2) , last_index_(i_sym,2)

            ! loop over j indeces

            do j = first_index_(j_sym,2) , last_index_(j_sym,2)

              ! ij geminal integral index
              ij_int = ints_%gemind(i,j)
              ! ij geminal density index
              ij_den = dens_%gemind(i,j)

              do k = first_index_(k_sym,2) , last_index_(k_sym,2)

                do l = first_index_(l_sym,2) , last_index_(l_sym,2)

                  ! kl geminal integral index
                  kl_int  = ints_%gemind(k,l)

                  ! ij geminal density index
                  kl_den  = dens_%gemind(k,l)
                  ! integral g(ij|kl) index 
                  int_ind = pq_index(ij_int,kl_int) + ij_int_offset
                  ! density d2(ij|kl) index 
                  den_ind = pq_index(ij_den,kl_den) + ij_den_offset
                  ! update 2-e active contribution
                  e_out   = e_out + int2(int_ind) * den2(den_ind) 

                end do ! end l loop

              end do ! end k loop

            end do ! end j_loop

          end  do ! end i loop

        end do ! end k_sym loop

      end  do ! end j_sym loop

    end do ! end i_sym loop

    e_out = 0.5_wp * e_out 

    return
  end subroutine compute_active_active_2e

end module focas_energy
