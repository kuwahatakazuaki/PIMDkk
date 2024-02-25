subroutine Uupdate
  use Parameters
  implicit none
  integer :: imode, iatom

  do imode = 1, nbead
    do iatom = 1, natom
      ur(:,iatom,imode) = ur(:,iatom,imode) + dt_ref*vur(:,iatom,imode)
      !ur(:,iatom,imode) = ur(:,iatom,imode) + 0.5 * dt_ref*vur(:,iatom,imode)
    enddo
  enddo
  return
end subroutine

subroutine Vupdate
  use Parameters
  implicit none
  integer :: imode, iatom

  do imode = 1, nbead
    do iatom = 1, natom
      vur(:,iatom,imode) = vur(:,iatom,imode) + 0.5d0*dt*fur(:,iatom,imode)/fictmass(iatom,imode)
    enddo
  enddo
  return
end subroutine
