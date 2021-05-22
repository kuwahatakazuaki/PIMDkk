subroutine init_position

! +++ gaussian distribution of normal modes +++
! +++ corresponding to free particle        +++

  use global_variable
  implicit none
  integer :: iatom, imode
  real(8) :: gasd


  if (Lrandomc .eqv. .True. ) then
    do imode = 2, Nbead
      do iatom = 1, Natom
!        usigma = dsqrt(1.d0/beta/omega_p2/dnmmass(iatom,imode))
!        call gasdev(gasd)
!        ux(iatom,imode) = usigma*gasd
!        call gasdev(gasd)
!        uy(iatom,imode) = usigma*gasd
!        call gasdev(gasd)
!        uz(iatom,imode) = usigma*gasd
      enddo
    enddo
  else
    u(:,:,2:Nbead) = 0.0d0
  end if


!  do imode = 2, nbead
!     do iatom = 1, natom
!        If(nrandomc==1) Then
!           usigma = dsqrt(1.d0/beta/omega_p2/dnmmass(iatom,imode))
!           call gasdev(gasd)
!           ux(iatom,imode) = usigma*gasd
!           call gasdev(gasd)
!           uy(iatom,imode) = usigma*gasd
!           call gasdev(gasd)
!           uz(iatom,imode) = usigma*gasd
!        Else   ! Default -> nrandomc = 0
!           ux(iatom,imode) = 0.d0
!           uy(iatom,imode) = 0.d0
!           uz(iatom,imode) = 0.d0
!        EndIf
!     enddo
!  enddo
!
!  Return
End Subroutine
