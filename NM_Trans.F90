subroutine nmtrans_ur2r
  use Parameters
  integer :: i, j, xyz
  r(:,:,:) = 0.0d0
  do i = 1, Natom
    do j = 1, Nbead
      do xyz = 1, 3
        r(xyz,i,j) = r(xyz,i,j) + dot_product(tnm(j,:),ur(xyz,i,:))
      end do
    end do
  end do
end subroutine nmtrans_ur2r

subroutine nmtrans_r2ur
  use Parameters
  integer :: i, j, xyz
  r(:,:,:) = 0.0d0
  do i = 1, Natom
    do j = 1, Nbead
      do xyz = 1, 3
        ur(xyz,i,j) = ur(xyz,i,j) + dot_product(tnminv(j,:),r(xyz,i,:))
      end do
    end do
  end do
end subroutine nmtrans_r2ur

subroutine nmtrans_fr2fur
  use Parameters
  implicit none
  integer :: i, j, k, imode, jmode, iatom

  fur(:,:,:) = 0.0d0

!    *  fu(i) = fu(i) + sum_j fx(j)*tnm(j,i)  */
  do iatom = 1, natom
    do imode = 1, nbead
      do jmode = 1, nbead
        fur(:,iatom,imode) = fur(:,iatom,imode) + fr(:,iatom,jmode)*tnm(jmode,imode)
      end do
    end do
  end do

return
end subroutine nmtrans_fr2fur

!Subroutine NM_Trans(NTmp)
!  Use Parameters
!  Implicit None
!  integer :: i, j, iatom, imode, jmode
!  integer, intent(in) :: NTmp
!
!  if (NTmp ==  0) then
!!  /*  from normal mode variables to real variables  *
!!   *  x(i) = x(i) + sum_j tnm(i,j)*u(j)             */
!    r(:,:,:) = 0.0d0
!    do iatom = 1, natom
!      do imode = 1, nbead
!        do jmode = 1, nbead
!          !x(iatom,imode) = x(iatom,imode) + tnm(imode,jmode)*ux(iatom,jmode)
!          !y(iatom,imode) = y(iatom,imode) + tnm(imode,jmode)*uy(iatom,jmode)
!          !z(iatom,imode) = z(iatom,imode) + tnm(imode,jmode)*uz(iatom,jmode)
!          r(1,iatom,imode) = r(1,iatom,imode) + tnm(imode,jmode)*ur(1,iatom,jmode)
!          r(2,iatom,imode) = r(2,iatom,imode) + tnm(imode,jmode)*ur(2,iatom,jmode)
!          r(3,iatom,imode) = r(3,iatom,imode) + tnm(imode,jmode)*ur(3,iatom,jmode)
!        enddo
!      enddo
!    enddo
!
!  elseif (NTmp == 1) then  ! We don't use this
!!     /*  from real variables to normal mode variables  *
!!      *  u(i) = u(i) + sum_j tnminv(i,j)*x(j)          */
!!     /*  initialize array  */
!    ur(:,:,:) = 0.0d0
!
!    do iatom = 1, natom
!      do imode = 1, nbead
!        do jmode = 1, nbead
!          !ux(iatom,imode) = ux(iatom,imode) + tnminv(imode,jmode)*x(iatom,jmode)
!          !uy(iatom,imode) = uy(iatom,imode) + tnminv(imode,jmode)*y(iatom,jmode)
!          !uz(iatom,imode) = uz(iatom,imode) + tnminv(imode,jmode)*z(iatom,jmode)
!          ur(1,iatom,imode) = ur(1,iatom,imode) + tnminv(imode,jmode)*r(1,iatom,jmode)
!          ur(2,iatom,imode) = ur(2,iatom,imode) + tnminv(imode,jmode)*r(2,iatom,jmode)
!          ur(3,iatom,imode) = ur(3,iatom,imode) + tnminv(imode,jmode)*r(3,iatom,jmode)
!        enddo
!      enddo
!    enddo
!  endif
!
!  Return
!End Subroutine
