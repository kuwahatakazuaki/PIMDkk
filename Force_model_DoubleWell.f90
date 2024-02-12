! Obtaining force(fx) from position(x)
! ++++++++++++++++++++++++++++
! ++++++ Double well 1D ++++++
! ++++++++++++++++++++++++++++
subroutine Force_model_DoubleWell1D
  Use Parameters, &
    !only: x, y, z, fx, fy, fz, Natom, Nbead, Eenergy, potential, &
    only: r, fr, Natom, Nbead, Eenergy, potential, &
          alabel, dp_inv, address, istepsv, &
          Lsave_force, &
          AUtoAng => bohr_inv, KtoAU => boltz, AngtoAU => bohr
  !Use Parameter_tk
  implicit none
  integer :: i, j, imode, iatom
  real(8) :: f_two(3), power
  integer :: Udis, Ucoor, Ufor, Uene

  real(8), parameter :: xwid = 0.5 * AngtoAU
  real(8), parameter :: height   = 1000 * KtoAU

  ! +++ Calculating Forcde which atom (i) feels from atom (j) +++
  fr(:,:,:) = 0.0d0
  do imode = 1, Nbead
    do i = 1, Natom
      fr(1,i,imode) = - 4.0 * height / xwid**4 * r(1,i,imode) * (r(1,i,imode)**2 - xwid**2)
    end do
  end do
  fr(:,:,:) = fr(:,:,:) * dp_inv
  ! +++ End Calculating Forcde which atom (i) feels from atom (j) +++

  ! +++ Calculating enetemp +++
  Eenergy(:) = 0.0d0
  do imode = 1, Nbead
    do i = 1, Natom
      Eenergy(imode) = Eenergy(imode) + height / xwid**4 * (r(1,i,imode)**2 - xwid**2)**2
    end do
  end do
  ! +++ End Calculating enetemp +++

  9998 format(3E23.15)
  9999 format(a2,1x,E15.9,1x,E15.9,1x,E15.9)
  return
end subroutine Force_model_DoubleWell1D
! ++++++++++++++++++++++++++++
! ++++++ Double well 1D ++++++
! ++++++++++++++++++++++++++++


! ++++++++++++++++++++++++++++
! ++++++ Double well 3D ++++++
! ++++++++++++++++++++++++++++
subroutine Force_model_DoubleWell3D
  Use Parameters, &
    !only: x, y, z, fx, fy, fz, Natom, Nbead, Eenergy, potential, &
    only: r, fr, Natom, Nbead, Eenergy, potential, &
          alabel, dp_inv, address, istepsv, &
          Lsave_force, &
          AUtoAng => bohr_inv, KtoAU => boltz, AngtoAU => bohr
  implicit none
  integer :: i, j, imode, iatom
  real(8) :: f_two(3), power
  !real(8) :: r(3,Natom,Nbead), f(3,Natom,Nbead)
  real(8) :: rij(3), dij2, dij1
  integer :: Udis, Ucoor, Ufor, Uene

  real(8), parameter :: xwid = 1.0 * AngtoAU
!   real(8), parameter :: height = 1000 * KtoAU
  real(8), parameter :: height = 0.003

  !! r(xyz,atom,bead)
  !r(1,:,:) =  x(:,:)
  !r(2,:,:) =  y(:,:)
  !r(3,:,:) =  z(:,:)

! open(newunit=Udis,file=trim(address)//'/dis.dat',position='append')
  ! +++ Calculating Forcde which atom (i) feels from atom (j) +++
  fr(:,:,:) = 0.0d0
  Eenergy(:) = 0.0d0
  do imode = 1, Nbead
    do i = 1, Natom
      do j = i+1, Natom
        rij(:) = r(:,i,imode) - r(:,j,imode)
        dij2 = dot_product(rij(:),rij(:))
        dij1 = dsqrt(dij2)
!        dij1 = dsqrt(dij2) - 1.0*AngtoAU    ! Shift the potential

        fr(:,i,imode) = - 4.0 * height / xwid**4 * (dij2 - xwid**2) * rij(:)
!        f(:,i,imode) = - 4.0 * height / xwid**4 * (dij1**2 - xwid**2) * rij(:)
        fr(:,j,imode) = (-1) * fr(:,i,imode)
      ! +++ Calculating enetemp +++
        Eenergy(imode) = Eenergy(imode) + height / xwid**4 * (dij2 - xwid**2)**2
!        Eenergy(imode) = Eenergy(imode) + height / xwid**4 * (dij1**2 - xwid**2)**2

!        write(Udis,*) dij1
      end do
    end do
  end do
!close(Udis)
  fr(:,:,:) = fr(:,:,:) * dp_inv
  ! +++ End Calculating Forcde which atom (i) feels from atom (j) +++


  9998 format(3E23.15)
  9999 format(a2,1x,E15.9,1x,E15.9,1x,E15.9)
  return
end subroutine Force_model_DoubleWell3D
! ++++++++++++++++++++++++++++
! ++++++ Double well 3D ++++++
! ++++++++++++++++++++++++++++
