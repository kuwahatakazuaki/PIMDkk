Subroutine Kinetic_Energy (Tmp)

  Use Parameters
  Implicit None
  Double Precision :: Tmp

  Tmp = 0.d0
  do imode = 1, nbead
     do iatom = 1, natom
        Tmp  = Tmp                                                      &
             + fictmass(iatom,imode)*vux(iatom,imode)*vux(iatom,imode)  &
             + fictmass(iatom,imode)*vuy(iatom,imode)*vuy(iatom,imode)  &
             + fictmass(iatom,imode)*vuz(iatom,imode)*vuz(iatom,imode)
     enddo
  enddo
  Return
End Subroutine
