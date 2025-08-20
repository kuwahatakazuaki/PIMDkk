subroutine Temp_ctr
  use Parameters
  implicit none
  Double Precision :: temp_scale, tempi
  !integer :: imode, iatom
  real(8) :: get_kinetic_ene

  tempi = get_kinetic_ene()
  tempi = 2.d0*tempi/(3.d0*dble(natom))/KtoAU
  tempi = tempi/dble(nbead)

  temp_scale = dsqrt(temperature/tempi)
  vur(:,:,:) = vur(:,:,:) * temp_scale

return
end subroutine Temp_ctr
