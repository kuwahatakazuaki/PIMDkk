! 2019/11/08 : Adding 3N/2beta, change the sign from Force to Potential
subroutine Virial_Estimator
  use Parameters
  implicit none
  integer          :: i, j, k, imode, iatom
  double precision :: e_virial1

  E_Virial=0.0D0
  do imode=1,nbead
    do iatom=1,natom
      !E_Virial = E_Virial + fx(iatom,imode) * (x(iatom,imode)-ux(iatom,1))
      !E_Virial = E_Virial + fy(iatom,imode) * (y(iatom,imode)-uy(iatom,1))
      !E_Virial = E_Virial + fz(iatom,imode) * (z(iatom,imode)-uz(iatom,1))
      E_Virial = E_Virial &
                 + dot_product(fr(:,iatom,imode),(r(:,iatom,imode)-ur(:,iatom,1)))
    enddo
  enddo

  E_Virial=E_Virial/2.0D+00

! Kuwahata 2019/11/08
  e_virial1 = 1.5d0*dble(natom)/beta
  E_Virial = e_virial1 - E_Virial
! End Kuwahata 2019/11/08

  return
end subroutine Virial_Estimator
