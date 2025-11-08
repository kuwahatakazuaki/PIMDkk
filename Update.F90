subroutine Uupdate
!subroutine update_pos_nor
  use Parameters
  implicit none
  integer :: imode, iatom

  do imode = 1, Nbead
    do iatom = 1, Natom
      ur(:,iatom,imode) = ur(:,iatom,imode) + dt_ref*vur(:,iatom,imode)
      !ur(:,iatom,imode) = ur(:,iatom,imode) + 0.5 * dt_ref*vur(:,iatom,imode)
    enddo
  enddo
  return
!end subroutine update_pos_nor
end subroutine

subroutine Vupdate
  use Parameters
  implicit none
  integer :: imode, iatom

  do imode = 1, Nbead
    do iatom = 1, Natom
      vur(:,iatom,imode) = vur(:,iatom,imode) + 0.5d0*dt*fur(:,iatom,imode)/fictmass(iatom,imode)
    enddo
  enddo
  return
end subroutine
