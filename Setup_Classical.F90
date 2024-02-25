Subroutine Setup_Classical
  use Parameters
  implicit none
  double precision :: rndnumber
  integer :: inhc, icolor, imode, iatom

!YK Initiate Random Number Generator
  Call RandomG(0,rndnumber)

!     /*   for multiple time step   */
  !dt_ref = dt
  if ( Isimulation == 10 ) then
    dt_ref = dt
  else
    dt_ref = dt/dble(nref)
  end if
  ista = 1
  iend = 1
!
!     /*   bath parameters for path integral MD   */
!
!YK different qmass required for centroid MD
  select case(Ncent)
    case(1)
  !If(NCent==1) Then
      qmcent11(1) = 3.d0*dble(natom)/beta/omega2
      do inhc=2,nnhc
         qmcent11(inhc) = 1.d0/beta/omega2
      enddo
  !EndIf
  !If(NCent==3) Then
    case(3)
      qmcent31(1) = 3.d0*dble(natom)/beta/omega2
      do inhc=2,nnhc
         qmcent31(inhc) = 1.d0/beta/omega_p2
      enddo
  end select
  !EndIf


!    /*   set parameters for higher order decomposition   *
!     *   of propagator (Yoshida and Suzuki parameters)   */

  if(nys==1) then
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

  do iatom = 1, natom
    fictmass(iatom,1) = physmass(iatom)
  enddo

  Return
End Subroutine
