subroutine Init_Mass
  use Parameters
  implicit none
  integer :: imode, iatom
  real(8) :: dp, di, twopi

  twopi=2.d0*pi
  dp = dble(Nbead)

! Kuwahata 2020/12/04  bug fixed
!     dnmmass(:,1) = 0.d0
!     dnmmass(:,Nbead) = 4.d0*physmass(iatom)
!  do iatom = 1, natom
!     do imode = 1, (Nbead-2)/2
!        di = dble(imode)
!        dnmmass(iatom,2*imode) = 2.d0*(1.d0-dcos(twopi*di/dp))*physmass(iatom)
!        dnmmass(iatom,2*imode+1) = dnmmass(iatom,2*imode)
!     enddo
!  enddo


  dnmmass(:,1) = 0.d0
  do iatom = 1, natom
    dnmmass(iatom,Nbead) = 4.d0*dp*physmass(iatom)
    do imode = 1, (Nbead-2)/2
      di = dble(imode)
      dnmmass(iatom,2*imode) = 2.d0*(1.d0-dcos(twopi*di/dp))*dp*physmass(iatom)
      dnmmass(iatom,2*imode+1) = dnmmass(iatom,2*imode)
    enddo
  enddo

! End Kuwahata 2020/12/04  bug fixed

!
!      /*  fictitious mass for centroid MD  */
!

  if ( Isimulation == 1 .or. Isimulation == 10 ) then
! +++ For RPMD or conventional MD simulation +++
    do iatom = 1, Natom
      fictmass(iatom,:) = physmass(iatom)
    enddo
! +++ For RPMD or conventional MD simulation +++
  !else if ( Isimulation == 0 .or. Isimulation == 2 ) then
  else
! +++ For PIMD or CMD simulation +++
    do iatom = 1, Natom
      fictmass(iatom,1) = physmass(iatom)
      do imode = 2, Nbead
        fictmass(iatom,imode) = gamma2*dnmmass(iatom,imode)
      enddo
    enddo
! +++ For CMD simulation +++
  end if
!  dnmmass(:,:) = dnmmass(:,:) * dp

  return
end subroutine Init_Mass
