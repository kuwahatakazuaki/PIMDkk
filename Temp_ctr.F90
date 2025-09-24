subroutine Temp_ctr
  use Parameters
  implicit none
  real(8) :: temp_scale, tempi
  real(8) :: get_kinetic_ene

  dkinetic = get_kinetic_ene()
  tempi = 2.d0*dkinetic*AU2K/dble(3*Natom*Nbead)
  !tempi = 2.d0*tempi/(3.d0*dble(natom))/K2AU
  !tempi = tempi/dble(nbead)

  temp_scale = dsqrt(temperature/tempi)
  vur(:,:,:) = vur(:,:,:) * temp_scale

return
end subroutine Temp_ctr
