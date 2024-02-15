!Subroutine Getfnm
subroutine nmtrans_force_r2ur
  use Parameters
  implicit none
  integer :: i, j, k, imode, jmode, iatom

  fur(:,:,:) = 0.0d0

!    *  fu(i) = fu(i) + sum_j fx(j)*tnm(j,i)  */
  Do iatom = 1, natom
    Do imode = 1, nbead
      Do jmode = 1, nbead
        !fux(iatom,imode) = fux(iatom,imode) + fx(iatom,jmode)*tnm(jmode,imode)
        !fuy(iatom,imode) = fuy(iatom,imode) + fy(iatom,jmode)*tnm(jmode,imode)
        !fuz(iatom,imode) = fuz(iatom,imode) + fz(iatom,jmode)*tnm(jmode,imode)
        fur(:,iatom,imode) = fur(:,iatom,imode) + fr(:,iatom,jmode)*tnm(jmode,imode)
      Enddo
    Enddo
  Enddo

Return
end subroutine nmtrans_force_r2ur
