Subroutine Temp_ctr
  use Parameters
  implicit none
  Double Precision :: temp_scale, tempi
  integer :: imode, iatom
  real(8) :: get_kinetic_ene

  !tempi=0.0D+00
  !call kinetic_energy(tempi)
  !tempi = 0.5d0*tempi
  tempi = get_kinetic_ene()
  tempi = 2.d0*tempi/(3.d0*dble(natom))/boltz
  tempi = tempi/dble(nbead)

  temp_scale = dsqrt(temperature/tempi)
  do imode = 1, nbead
    do iatom = 1, natom
      !vux(iatom,imode) = vux(iatom,imode)*temp_scale
      !vuy(iatom,imode) = vuy(iatom,imode)*temp_scale
      !vuz(iatom,imode) = vuz(iatom,imode)*temp_scale
      vur(:,iatom,imode) = vur(:,iatom,imode)*temp_scale
    enddo
  enddo

Return
End Subroutine
