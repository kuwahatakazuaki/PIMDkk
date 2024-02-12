Subroutine Setup_MPI
  Use Parameters
  Implicit None
  double precision :: rndnumber
  integer :: inhc, icolor, imode

!YK Initiate Random Number Generator
  Call RandomG(0,rndnumber)

!     /*   for multiple time step   */

  dt_ref = dt/dble(nref)
!
!     /*   bath parameters for path integral MD   */
!
! /*   thermostat attached to centroids   */
  select case(Ncent)
    case(1)
      qmcent11(1) = 3.d0*dble(natom)/beta/omega2
      ! qmcent11(1) = 3.d0/beta/omega2 ??
      do inhc=2,nnhc
        qmcent11(inhc) = 1.d0/beta/omega2
      enddo
    case(3)
      qmcent31(1) = 3.d0*dble(natom)/beta/omega2
      ! qmcent31(1) = 1.0d0/beta/omega2 ??
      do inhc=2,nnhc
        qmcent31(inhc) = 1.d0/beta/omega2
      enddo
  end select


!YK different qmass required for centroid MD
! /*   thermostat attached to non-centroid modes   */
  qmass(1) = 0.0d0
  do imode = 2, nbead
     qmass(imode) = 1.d0/beta/omega_p2
  enddo
!
!     /*   For centroid MD, qmass should be scaled by gamma2, since   *
!      *   natural frequencies of modes are omega_p**2/gamma**2       */
!
!YK If CMD
!     /*   adiabaticity parameter for centroid MD   */
  gamma2 = gamma*gamma
  If( Isimulation == 2 ) then
    do imode = 2, nbead
      qmass(imode) = gamma2*qmass(imode)
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


Return
End Subroutine
