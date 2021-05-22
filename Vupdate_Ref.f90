Subroutine Vupdate_Ref

  Use Parameters
  Implicit None

  Integer                            :: i

  do imode = 1, nbead
     do iatom = 1, natom
        vux(iatom,imode) = vux(iatom,imode) + 0.5d0*dt_ref*fux_ref(iatom,imode)/fictmass(iatom,imode)
        vuy(iatom,imode) = vuy(iatom,imode) + 0.5d0*dt_ref*fuy_ref(iatom,imode)/fictmass(iatom,imode)
        vuz(iatom,imode) = vuz(iatom,imode) + 0.5d0*dt_ref*fuz_ref(iatom,imode)/fictmass(iatom,imode)
     enddo
  enddo


End Subroutine
