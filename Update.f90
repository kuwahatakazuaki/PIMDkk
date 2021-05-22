Subroutine Uupdate
  Use Parameters
  Implicit None

  do imode = 1, nbead
     do iatom = 1, natom
        ux(iatom,imode) = ux(iatom,imode) + dt_ref*vux(iatom,imode)
        uy(iatom,imode) = uy(iatom,imode) + dt_ref*vuy(iatom,imode)
        uz(iatom,imode) = uz(iatom,imode) + dt_ref*vuz(iatom,imode)
     enddo
  enddo
  Return
End Subroutine

Subroutine Vupdate
  Use Parameters
  Implicit None

  do imode = 1, nbead
     do iatom = 1, natom
        vux(iatom,imode) = vux(iatom,imode) + 0.5d0*dt*fux(iatom,imode)/fictmass(iatom,imode)
        vuy(iatom,imode) = vuy(iatom,imode) + 0.5d0*dt*fuy(iatom,imode)/fictmass(iatom,imode)
        vuz(iatom,imode) = vuz(iatom,imode) + 0.5d0*dt*fuz(iatom,imode)/fictmass(iatom,imode)
     enddo
  enddo
  Return
End Subroutine
