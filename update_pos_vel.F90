subroutine update_pos_nor
  use Parameters
  implicit none
  !integer :: Imode, Iatom
  !do Imode = 1, Nbead
  !  do Iatom = 1, Natom
  !    ur(:,Iatom,Imode) = ur(:,Iatom,Imode) + dt_ref*vur(:,Iatom,Imode)
  !  enddo
  !enddo
  ur(:,:,:) = ur(:,:,:) + dt_ref*vur(:,:,:)
  return
end subroutine update_pos_nor

subroutine update_pos_car
  use Parameters
  implicit none
  !ur(:,:,:) = ur(:,:,:) + dt_ref*vur(:,:,:)
  r(:,:,:) = r(:,:,:) + dt_ref*vr(:,:,:)
  return
end subroutine update_pos_car


subroutine update_vel_nor
  use Parameters
  implicit none
  integer :: Imode, Iatom

  !do Imode = 1, Nbead
  !  do Iatom = 1, Natom
  !    vur(:,Iatom,Imode) = vur(:,Iatom,Imode) + 0.5d0*dt*fur(:,Iatom,Imode)/fictmass(Iatom,Imode)
  !  enddo
  !enddo
  vur(:,:,:) = vur(:,:,:) + 0.5d0 * dt * fur(:,:,:) / spread(fictmass, dim=1, ncopies=3)
  return
end subroutine update_vel_nor


subroutine update_vel_ref_nor
  use Parameters
  use utility, only: program_abort
  implicit none
  integer:: j, i

  !do j = 2, Nbead
  !  do i = 1, Natom
  !    vur(:,i,j) = vur(:,i,j) + 0.5d0*dt_ref*fur_ref(:,i,j)/fictmass(i,j)
  !  end do
  !end do
  vur(:,:,:) = vur(:,:,:) + 0.5d0 * dt_ref * fur_ref(:,:,:) / spread(fictmass, dim=1, ncopies=3)
end subroutine update_vel_ref_nor

