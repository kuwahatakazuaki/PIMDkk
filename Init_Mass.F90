subroutine Init_Mass
  use Parameters
  implicit none
  integer :: Imode, Iatom
  real(8) :: dp, di, twopi

  twopi=2.d0*pi
  dp = dble(Nbead)

! Kuwahata 2020/12/04  bug fixed
!     dnmmass(:,1) = 0.d0
!     dnmmass(:,Nbead) = 4.d0*physmass(Iatom)
!  do Iatom = 1, Natom
!     do Imode = 1, (Nbead-2)/2
!        di = dble(Imode)
!        dnmmass(Iatom,2*Imode) = 2.d0*(1.d0-dcos(twopi*di/dp))*physmass(Iatom)
!        dnmmass(Iatom,2*Imode+1) = dnmmass(Iatom,2*Imode)
!     enddo
!  enddo


  dnmmass(:,1) = 0.d0
  do Iatom = 1, Natom
    dnmmass(Iatom,Nbead) = 4.d0*dp*physmass(Iatom)
    do Imode = 1, (Nbead-2)/2
      di = dble(Imode)
      dnmmass(Iatom,2*Imode) = 2.d0*(1.d0-dcos(twopi*di/dp))*dp*physmass(Iatom)
      dnmmass(Iatom,2*Imode+1) = dnmmass(Iatom,2*Imode)
    enddo
  enddo

! End Kuwahata 2020/12/04  bug fixed

!
!      /*  fictitious mass for centroid MD  */
!

  if ( Isimulation == 1 .or. Isimulation == 10 ) then
! +++ For RPMD or conventional MD simulation +++
    do Iatom = 1, Natom
      fictmass(Iatom,:) = physmass(Iatom)
    enddo
! +++ For RPMD or conventional MD simulation +++
  !else if ( Isimulation == 0 .or. Isimulation == 2 ) then
  else
! +++ For PIMD or CMD simulation +++
    do Iatom = 1, Natom
      fictmass(Iatom,1) = physmass(Iatom)
      do Imode = 2, Nbead
        fictmass(Iatom,Imode) = gamma2*dnmmass(Iatom,Imode)
      enddo
    enddo
! +++ For CMD simulation +++
  end if
!  dnmmass(:,:) = dnmmass(:,:) * dp

  return
end subroutine Init_Mass
