subroutine Setup_time_mass
  use Parameters
  use utility,    only: set_random_seed
  implicit none
  double precision :: rndnumber
  integer :: Inhc, Icolor, Imode

!YK Initiate Random Number Generator
  !Call RandomG(0,rndnumber)
  call set_random_seed(Lfixed_random = .True.)

!     /*   for multiple time step   */
  if ( Isimulation == 10 ) then
    dt_ref = dt
    ista = 1
    iend = 1
  else
    dt_ref = dt/dble(Nref)
  end if
!
!     /*   bath parameters for path integral MD   */
!
! /*   thermostat attached to centroids   */
  select case(Ncent)
    case(1)
      !qmcent11(1) = 3.d0*dble(Natom)/beta/omega2
      qmcent11(1) = dble(Ndim*Natom)/beta/omega2
      do Inhc = 2, Nnhc
        qmcent11(Inhc) = 1.d0/beta/omega2
      enddo
    case(3)
      qmcent31(1) = dble(Ndim*Natom)/beta/omega2
      !qmcent31(1) = 3.d0*dble(Natom)/beta/omega2
      ! qmcent31(1) = 1.0d0/beta/omega2 ??
      do Inhc = 2, Nnhc
        qmcent31(Inhc) = 1.d0/beta/omega2
      enddo
  end select


! /*   thermostat attached to non-centroid modes   */
  if ( Isimulation /= 10 ) then
    qmass(1) = 0.0d0
    do Imode = 2, Nbead
      qmass(Imode) = 1.d0/beta/omega_p2
    enddo
  end if
!
!  For centroid MD, qmass should be scaled by gamma2, since
!  natural frequencies of modes are omega_p**2/gamma**2
!  adiabaticity parameter for centroid MD
  gamma2 = gamma1*gamma1
  If( Isimulation == 2 ) then
    do Imode = 2, Nbead
      qmass(Imode) = gamma2*qmass(Imode)
    enddo
  EndIf

!
!     /*   set parameters for higher order decomposition   *
!      *   of propagator (Yoshida and Suzuki parameters)   */
!
  if (nys==1) then
    ysweight(1) = 1.d0
  else if (nys==3) then
    ysweight(1) = 1.d0/(2.d0 - 2.d0**(1.d0/3.d0))
    ysweight(2) = 1.d0 - 2.d0*ysweight(1)
    ysweight(3) = ysweight(1)
  else if (nys==5) then
    ysweight(1) = 1.d0/(4.d0 - 4.d0**(1.d0/3.d0))
    ysweight(2) = ysweight(1)
    ysweight(3) = 1.d0 - 4.d0*ysweight(1)
    ysweight(4) = ysweight(1)
    ysweight(5) = ysweight(1)
  end if


  return
end subroutine Setup_time_mass
