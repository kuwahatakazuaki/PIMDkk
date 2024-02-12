Subroutine Init_Mass
  Use Parameters
  Implicit None
  Double Precision :: di,twopi
  integer :: imode, iatom

  twopi=2.d0*pi
  dp = dble(nbead)

! Kuwahata 2020/12/04  bug fixed
!     dnmmass(:,1) = 0.d0
!     dnmmass(:,nbead) = 4.d0*physmass(iatom)
!  do iatom = 1, natom
!     do imode = 1, (nbead-2)/2
!        di = dble(imode)
!        dnmmass(iatom,2*imode) = 2.d0*(1.d0-dcos(twopi*di/dp))*physmass(iatom)
!        dnmmass(iatom,2*imode+1) = dnmmass(iatom,2*imode)
!     enddo
!  enddo


  do iatom = 1, natom
     dnmmass(iatom,1) = 0.d0
     dnmmass(iatom,nbead) = 4.d0*dp*physmass(iatom)
     do imode = 1, (nbead-2)/2
        di = dble(imode)
        dnmmass(iatom,2*imode) = 2.d0*(1.d0-dcos(twopi*di/dp))*dp*physmass(iatom)
        dnmmass(iatom,2*imode+1) = dnmmass(iatom,2*imode)
     enddo
  enddo

! End Kuwahata 2020/12/04  bug fixed

!
!      /*  fictitious mass for centroid MD  */
!
  If ( Isimulation == 1) Then

    do iatom = 1, natom
      !fictmass(iatom,1) = physmass(iatom)
      !do imode = 2, nbead
      !  fictmass(iatom,imode) = physmass(iatom)
      !enddo
      fictmass(iatom,:) = physmass(iatom)
    enddo

  Else

    do iatom = 1, natom
      fictmass(iatom,1) = physmass(iatom)
      do imode = 2, nbead
        fictmass(iatom,imode) = gamma2*dnmmass(iatom,imode)
      enddo
    enddo

  ENDIF
!  dnmmass(:,:) = dnmmass(:,:) * dp

  Return
End Subroutine
