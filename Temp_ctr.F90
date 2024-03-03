Subroutine Temp_ctr
  use Parameters
  implicit none
  Double Precision :: temp_scale, tempi
  integer :: imode, iatom
  real(8) :: get_kinetic_ene

  tempi = get_kinetic_ene()
  tempi = 2.d0*tempi/(3.d0*dble(natom))/KtoAU
  tempi = tempi/dble(nbead)

  temp_scale = dsqrt(temperature/tempi)
  vur(:,:,:) = vur(:,:,:) * temp_scale
  !do imode = 1, nbead
  !  do iatom = 1, natom
  !    vur(:,iatom,imode) = vur(:,iatom,imode)*temp_scale
  !  enddo
  !enddo

Return
End Subroutine
