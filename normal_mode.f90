subroutine normal_mode_matrix
  use global_variable, only: Nbead,  pi, &
                             tnm, tnminv, myrank
  implicit none
  integer :: i, j
  real(8) :: signf
  real(8) :: sdbinv, pref
  real(8) :: di, dj
  real(8) :: temp
  real(8) :: sqp, sqpinv
  real(8), allocatable :: unitary(:,:), uni_inv(:,:)
  allocate( unitary(Nbead, Nbead),       source=0.0d0 )
  allocate( uni_inv(Nbead, Nbead),       source=0.0d0 )

  sdbinv = 1.0d0 / dsqrt( dble(Nbead) )
  pref   = dsqrt(2.0d0 / dble(Nbead))

  signf = -1.0d0
  do i = 1, Nbead
    unitary(i,1) = sdbinv
    unitary(i,Nbead) = signf * sdbinv
    signf = signf * (-1.0d0)
  end do

  do i = 1, (Nbead-2) / 2
    di = dble(i)
    do j = 1, Nbead
      dj = dble(j)
      temp = 2.0d0 * pi * di * dj / dble(Nbead)
      unitary(j, 2*i  ) = pref * dcos(temp)
      unitary(j, 2*i+1) = pref * dsin(temp)
    end do
  end do

  do i = 1, Nbead
    do j = 1, Nbead
      uni_inv(j,i) = unitary(i,j)
    end do
  end do

  sqp    = sqrt(dble(Nbead))
  sqpinv = 1.0d0 / dsqrt(dble(Nbead))
  tnm(:,:) = sqp * unitary(:,:)
  tnminv(:,:) = sqpinv * uni_inv(:,:)
!  unitary(:,:) = sqp    * unitary(:,:)
!  uni_inv(:,:) - sqpinv * uni_inv(:,:)

!if (myrank == 0) then
!  do i = 1, Nbead
!    do j = 1, Nbead
!      print '(2I3,F10.5)', i,j, unitary(i,j)/pref
!    end do
!  end do
!  print *, "pref is ", pref
!end if

end subroutine normal_mode_matrix






! --- x(i) = x(i) + sum_j tnm(i,j)*u(j) --- 
subroutine normal_mode_trans
  use global_variable
  use mpi
  implicit none
  integer :: i, j, xyz

  r(:,:,:) = 0.0d0
  do i = 1, Natom
    do j = 1, Nbead
      do xyz = 1, 3
        r(xyz,i,j) = r(xyz,i,j) + dot_product(tnm(j,:),u(xyz,i,:))
      end do
    end do
  end do
end subroutine normal_mode_trans

!  Integer       :: NTmp, jmode
!  integer :: i, j
!
!
!  if (NTmp ==  0) then
!!     /*  from normal mode variables to real variables  *
!!      *  x(i) = x(i) + sum_j tnm(i,j)*u(j)             */
!!     /*  initialize array  */
!     do imode = 1, nbead
!        do iatom = 1, natom
!           x(iatom,imode) = 0.d0
!           y(iatom,imode) = 0.d0
!           z(iatom,imode) = 0.d0
!        enddo
!     enddo
!
!     do iatom = 1, natom
!        do imode = 1, nbead
!           do jmode = 1, nbead
!              x(iatom,imode) = x(iatom,imode) + tnm(imode,jmode)*ux(iatom,jmode)
!              y(iatom,imode) = y(iatom,imode) + tnm(imode,jmode)*uy(iatom,jmode)
!              z(iatom,imode) = z(iatom,imode) + tnm(imode,jmode)*uz(iatom,jmode)
!           enddo
!        enddo
!     enddo
!
!  elseif (NTmp == 1) then  ! We don't use this
!!     /*  from real variables to normal mode variables  *
!!      *  u(i) = u(i) + sum_j tnminv(i,j)*x(j)          */
!!     /*  initialize array  */
!     do imode = 1, nbead
!        do iatom = 1, natom
!           ux(iatom,imode) = 0.d0
!           uy(iatom,imode) = 0.d0
!           uz(iatom,imode) = 0.d0
!        enddo
!     enddo
!
!     do iatom = 1, natom
!        do imode = 1, nbead
!           do jmode = 1, nbead
!              ux(iatom,imode) = ux(iatom,imode) + tnminv(imode,jmode)*x(iatom,jmode)
!              uy(iatom,imode) = uy(iatom,imode) + tnminv(imode,jmode)*y(iatom,jmode)
!              uz(iatom,imode) = uz(iatom,imode) + tnminv(imode,jmode)*z(iatom,jmode)
!           enddo
!        enddo
!     enddo
!  endif


! --- fu(i) = fu(i) + sum_j fx(j)*tnm(j,i) ---
subroutine normal_mode_force
  use global_variable
  implicit none
  integer :: i, j, xyz


  fu(:,:,:) = 0.0d0
  do i = 1, Natom
    do j = 1, Nbead
      do xyz = 1, 3
      fu(xyz,i,j) = fu(xyz,i,j) + dot_product(f(xyz,i,:),tnm(:,j))
      end do
    end do
  end do
end subroutine normal_mode_force


