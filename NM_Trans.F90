subroutine nmtrans_ur2r
!      r(i) = r(i) + sum_j tnm(i,j)*u(j)
  use Parameters
  integer :: i, j, xyz
  !r(:,:,:) = 0.0d0
  do i = 1, Natom
    !do j = 1, Nbead
    !  do xyz = 1, Ndim
    !    r(xyz,i,j) = r(xyz,i,j) + dot_product(tnm(j,:),ur(xyz,i,:))
    !  end do
    !end do
    r(:, i, :) = matmul( ur(:, i, :), transpose(tnm) )
  end do
end subroutine nmtrans_ur2r

subroutine nmtrans_r2ur
!       u(i) = u(i) + sum_j tnminv(i,j)*x(j)
  use Parameters
  integer :: i, j, xyz
  !ur(:,:,:) = 0.0d0
  !do i = 1, Natom
  !  do j = 1, Nbead
  !    do xyz = 1, Ndim
  !      ur(xyz,i,j) = ur(xyz,i,j) + dot_product(tnminv(j,:),r(xyz,i,:))
  !    end do
  !  end do
  !end do
  do i = 1, Natom
    ur(:, i, :) = matmul( r(:, i, :), transpose(tnminv) )
  end do
end subroutine nmtrans_r2ur

subroutine nmtrans_fr2fur
!    *  fu(i) = fu(i) + sum_j fx(j)*tnm(j,i)  */
  use Parameters
  implicit none
  integer :: Iatom

  !fur(:,:,:) = 0.0d0
  !do iatom = 1, natom
  !  do imode = 1, nbead
  !    do jmode = 1, nbead
  !      fur(:,iatom,imode) = fur(:,iatom,imode) + fr(:,iatom,jmode)*tnm(jmode,imode)
  !    end do
  !  end do
  !end do
  do Iatom = 1, Natom
    fur(:, Iatom, :) = matmul( fr(:, Iatom, :), tnm )
  end do

return
end subroutine nmtrans_fr2fur

