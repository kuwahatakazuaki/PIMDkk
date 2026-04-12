subroutine update_pos_nor
  use Parameters
  implicit none
  ur(:,:,:) = ur(:,:,:) + dt_ref*vur(:,:,:)
  return
end subroutine update_pos_nor

subroutine update_pos_car
  use Parameters
  implicit none
  r(:,:,:) = r(:,:,:) + dt_ref*vr(:,:,:)
  return
end subroutine update_pos_car


subroutine update_vel_nor
  use Parameters
  implicit none
  !integer :: Imode, Iatom
  !do Imode = 1, Nbead
  !  do Iatom = 1, Natom
  !    vur(:,Iatom,Imode) = vur(:,Iatom,Imode) + 0.5d0*dt*fur(:,Iatom,Imode)/fictmass(Iatom,Imode)
  !  enddo
  !enddo
  vur(:,:,:) = vur(:,:,:) + 0.5d0 * dt * fur(:,:,:) / spread(fictmass, dim=1, ncopies=Ndim)
  return
end subroutine update_vel_nor


subroutine update_vel_ref_nor
  use Parameters
  use utility, only: program_abort
  implicit none
  !integer:: j, i
  !do j = 2, Nbead
  !  do i = 1, Natom
  !    vur(:,i,j) = vur(:,i,j) + 0.5d0*dt_ref*fur_ref(:,i,j)/fictmass(i,j)
  !  end do
  !end do
  vur(:,:,:) = vur(:,:,:) + 0.5d0 * dt_ref * fur_ref(:,:,:) / spread(fictmass, dim=1, ncopies=Ndim)
end subroutine update_vel_ref_nor


subroutine update_pos_vel_analy
  use Parameters
  implicit none
  real(8) :: omega0(Ndim,Natom,2:Nbead)
  real(8) :: cos0(Ndim,Natom,2:Nbead), sin0(Ndim,Natom,2:Nbead)
  real(8) :: ur0(Ndim,Natom,2:Nbead)

  ! for centroid
  ur(:,:,1) = ur(:,:,1) + vur(:,:,1) * dt

  ! for non-centroid
  omega0(:,:,2:Nbead) = spread( sqrt( dnmmass(:,2:Nbead)/fictmass(:,2:Nbead) * omega_p2 ), dim=1, ncopies=Ndim )
  cos0(:,:,2:Nbead) = cos(omega0(:,:,2:Nbead)*dt)
  sin0(:,:,2:Nbead) = sin(omega0(:,:,2:Nbead)*dt)
  ur0(:,:,2:Nbead)  = ur(:,:,2:Nbead)

  ur(:,:,2:Nbead)  = ur0(:,:,2:Nbead)*cos0(:,:,2:Nbead) + vur(:,:,2:Nbead)*sin0(:,:,2:Nbead)/omega0(:,:,2:Nbead)
  vur(:,:,2:Nbead) = vur(:,:,2:Nbead)*cos0(:,:,2:Nbead) - ur0(:,:,2:Nbead)*sin0(:,:,2:Nbead)*omega0(:,:,2:Nbead)

end subroutine update_pos_vel_analy

