subroutine force_spcf
!===================================================================
!
! modified flexible SPC model (SPC/F2) for water:
!    L. Lobaugh and G. A. Voth, J. Chem. Phys. 106, 2400 (1997).
! NOTE:  the atoms must be sorted as O, H, H, O, H, H, ...
!
!===================================================================
  use Parameters, only: &
    r, fr, pi, dipr => dipoler, Natom, Nbead, pot => Eenergy, &
    Ista, Iend
  use utility, only: program_abort

  implicit none
      real(8), parameter:: au_length = 0.529177249d-10
      real(8), parameter:: au_mass   = 9.1093897d-31
      real(8), parameter:: au_energy = 4.3597482d-18
      real(8), parameter:: au_time   = 0.024188843d-15
      real(8), parameter:: au_charge = 1.60217646d-19
      real(8), parameter:: boltz     = 0.316682968d-5   ! KtoAU

  integer ::  i, j, k, m, i_o, i_h, i_g
  real(8) ::  rho_w, d_w, b_oh, b_hh, b_const, c_const, d_const,  &
              sigma, epsilon, q_o, q_h, es6, es12, dvdr

  real(8) ::  rr(3), r1, r2, rinv, rinv2, rinv6, rinv12

  real(8) ::  rr_oh(3), r2_oh, r_oh,   &
              rr_og(3), r2_og, r_og,   &
              rr_hh(3), r2_hh, r_hh

  real(8) :: qi = 0.d0, qj = 0.d0

!-------------------------------------------------------------------
! /*   potential parameters                                       */
!-------------------------------------------------------------------


  rho_w   =   2.361d+10  *au_length
  d_w     =   0.708d-18  /au_energy

  b_oh    =   1.000d-10  /au_length
  b_hh    =   b_oh *sin(pi*108.d0/2.d0/180.d0) *2.d0

  b_const =   1.803d+02  /au_energy *au_length**2
  c_const = - 1.469d+02  /au_energy *au_length**2
  d_const =   0.776d+02  /au_energy *au_length**2

  sigma   =   3.165d-10 /au_length
  epsilon =   78.22d0   *boltz

  q_o     = - 0.82d0
  q_h     = + 0.41d0

  pot(:)    = 0.0d0
  fr(:,:,:) = 0.0d0
  dipr(:,:) = 0.0d0

!-------------------------------------------------------------------
! /*   intra-molecular term                                       */
!-------------------------------------------------------------------

  !do k = 1, nbead
  do k = Ista, Iend
     do i = 1, natom, 3

        i_o = i
        i_h = i + 1
        i_g = i + 2

        rr_oh(:) = r(:,i_o,k) - r(:,i_h,k)
        rr_og(:) = r(:,i_o,k) - r(:,i_g,k)
        rr_hh(:) = r(:,i_h,k) - r(:,i_g,k)

        !r2_oh = rx_oh*rx_oh + ry_oh*ry_oh + rz_oh*rz_oh
        r2_oh = dot_product(rr_oh(:),rr_oh(:))
        r_oh  = sqrt(r2_oh)

        !r2_og = rx_og*rx_og + ry_og*ry_og + rz_og*rz_og
        r2_og = dot_product(rr_og(:),rr_og(:))
        r_og  = sqrt(r2_og)

        !r2_hh = rx_hh*rx_hh + ry_hh*ry_hh + rz_hh*rz_hh
        r2_hh = dot_product(rr_hh(:),rr_hh(:))
        r_hh  = sqrt(r2_hh)

        pot(k) = pot(k)                                 &
           + rho_w*rho_w*d_w*(r_oh-b_oh)*(r_oh-b_oh)    &
           + rho_w*rho_w*d_w*(r_og-b_oh)*(r_og-b_oh)    &
           + b_const*0.5d0*(r_hh-b_hh)*(r_hh-b_hh)      &
           + c_const*(r_oh+r_og-2.d0*b_oh)*(r_hh-b_hh)  &
           + d_const*(r_oh-b_oh)*(r_og-b_oh)

        fr(:,i_o,k) = fr(:,i_o,k)                             &
           - 2.d0*rho_w*rho_w*d_w*rr_oh(:)/r_oh*(r_oh-b_oh)   &
           - 2.d0*rho_w*rho_w*d_w*rr_og(:)/r_og*(r_og-b_oh)   &
           - c_const*rr_oh(:)/r_oh*(r_hh-b_hh)                &
           - c_const*rr_og(:)/r_og*(r_hh-b_hh)                &
           - d_const*rr_oh(:)/r_oh*(r_og-b_oh)                &
           - d_const*rr_og(:)/r_og*(r_oh-b_oh)

        fr(:,i_h,k) = fr(:,i_h,k)                             &
           + 2.d0*rho_w*rho_w*d_w*rr_oh(:)/r_oh*(r_oh-b_oh)   &
           - b_const*rr_hh(:)/r_hh*(r_hh-b_hh)                &
           + c_const*rr_oh(:)/r_oh*(r_hh-b_hh)                &
           - c_const*(r_oh+r_og-2.d0*b_oh)*rr_hh(:)/r_hh      &
           + d_const*rr_oh(:)/r_oh*(r_og-b_oh)

        fr(:,i_g,k) = fr(:,i_g,k)                              &
           + 2.d0*rho_w*rho_w*d_w*rr_og(:)/r_og*(r_og-b_oh)    &
           + b_const*rr_hh(:)/r_hh*(r_hh-b_hh)                 &
           + c_const*rr_og(:)/r_og*(r_hh-b_hh)                 &
           + c_const*(r_oh+r_og-2.d0*b_oh)*rr_hh(:)/r_hh       &
           + d_const*rr_og(:)/r_og*(r_oh-b_oh)

     end do
  end do

!-------------------------------------------------------------------
! /*   inter-molecular Lennard-Jones term                         */
!-------------------------------------------------------------------

  es6  = epsilon*sigma**6
  es12 = epsilon*sigma**12

  !do k = 1, nbead
  do k = Ista, Iend

     do i = 1, natom, 3

        m = (i-1)/3

        do j = (m+1)*3+1, natom, 3

           rr(:) = r(:,i,k) - r(:,j,k)

           !r2     =  rx*rx + ry*ry + rz*rz
           r2     =  dot_product(rr(:),rr(:))
           r1      =  sqrt(r2)
           rinv   =  1.d0/r1
           rinv2  =  rinv*rinv
           rinv6  =  rinv2*rinv2*rinv2
           rinv12 =  rinv6*rinv6

           pot(k) = pot(k) - 4.d0*es6*rinv6 + 4.d0*es12*rinv12

           dvdr = + 24.d0*es6*rinv6*rinv - 48.d0*es12*rinv12*rinv

           fr(:,i,k) = fr(:,i,k) - dvdr*rr(:)*rinv

           fr(:,j,k) = fr(:,j,k) + dvdr*rr(:)*rinv

        end do

     end do

  end do

!-------------------------------------------------------------------
! /*   inter-molecular Coulomb term                               */
!-------------------------------------------------------------------

  !do k = 1, nbead
  do k = Ista, Iend

     do i = 1, natom

        m = (i-1)/3

        do j = (m+1)*3+1, natom

           rr(:) = r(:,i,k) - r(:,j,k)

           if ( mod( i, 3 ) .eq. 1 ) qi = q_o
           if ( mod( i, 3 ) .eq. 2 ) qi = q_h
           if ( mod( i, 3 ) .eq. 0 ) qi = q_h

           if ( mod( j, 3 ) .eq. 1 ) qj = q_o
           if ( mod( j, 3 ) .eq. 2 ) qj = q_h
           if ( mod( j, 3 ) .eq. 0 ) qj = q_h

           !r2    =  rx*rx + ry*ry + rz*rz
           r2     =  dot_product(rr(:),rr(:))
           r1    =  sqrt(r2)
           rinv  =  1.d0/r1
           rinv2 =  rinv*rinv

           pot(k) = pot(k) + qi*qj*rinv

           dvdr = - qi*qj*rinv2

           fr(:,i,k) = fr(:,i,k) - dvdr*rr(:)*rinv

           fr(:,j,k) = fr(:,j,k) + dvdr*rr(:)*rinv

        end do

     end do

  end do

!-------------------------------------------------------------------
! /*   dipole moment                                              */
!-------------------------------------------------------------------

  !do k = 1, nbead
  do k = Ista, Iend

     dipr(:,k) = 0.d0

     do i = 1, natom, 3
        i_o = i
        i_h = i + 1
        i_g = i + 2

        dipr(:,k) = dipr(:,k) + q_o*r(:,i_o,k) + q_h*(r(:,i_h,k)+r(:,i_g,k))
     end do

  end do

  return
end subroutine force_spcf



!subroutine force_water
!  use common_variables, only :  &
!     model_water, iboundary, iounit
!
!  implicit none
!  integer, save :: iset = 0
!
!!-------------------------------------------------------------------
!! /*   select water model                                         */
!!-------------------------------------------------------------------
!
!  if ( iset .eq. 0 ) then
!     call read_int1 ( model_water, '<model_water>', 13, iounit )
!     iset = 1
!  end if
!
!  if ( iboundary .ne. 0 ) then
!     write( 6, '(a)' ) 'Error - Water potentials are available for free boundary only.'
!     write( 6, '(a)' )
!     call error_handling( 1, 'subroutine force_water', 22 )
!  end if
!
!  if ( model_water .eq. 1 ) then
!     call force_spcf
!  else if ( model_water .eq. 2 ) then
!     stop 'Not Updated kk'
!     call force_rwk
!  end if
!
!  return
!end subroutine force_water

