subroutine Temp_ctr
  use Parameters
  implicit none
  real(8) :: temp_scale, tempi
  real(8) :: get_kinetic_ene

  dkinetic = get_kinetic_ene()
  tempi = 2.d0*dkinetic*AU2K/dble(Ndof)
  !tempi = 2.d0*dkinetic*AU2K/dble(Ndim*Natom*Nbead)

  temp_scale = dsqrt(temperature/tempi)
  vur(:,:,:) = vur(:,:,:) * temp_scale

return
end subroutine Temp_ctr
