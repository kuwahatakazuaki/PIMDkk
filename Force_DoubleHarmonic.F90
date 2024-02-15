! Obtaining force(fx) from position(x)
! ++++++++++++++++++++++++++++
! ++++++ Double Harmonic ++++++
! ++++++++++++++++++++++++++++
subroutine Force_DoubleHarmonic
  Use Parameters, &
    only: x, y, z, fx, fy, fz, Natom, Nbead, Eenergy, potential, &
          alabel, dp_inv, address, istepsv, &
          Lsave_force, &
          AUtoAng => bohr_inv, KtoAU => boltz, AngtoAU => bohr
  implicit none
  integer :: i, j, imode, iatom
  real(8) :: r(3,Natom,Nbead), f(3,Natom,Nbead)
  integer :: Upre, Udip, Uchar, Uhfc, Ucoor, Ufor, Uene
!  real(8) :: rij(3), dij2, dij1
  real(8) :: r12(3), r32(3)
  real(8) :: dis12, dis32
  real(8) :: f12(3), f32(3)

  real(8) :: cons = 1.0 ! reading from input
  real(8) :: req = 1.0d0 * AngtoAU

!  real(8), parameter :: xwid = 0.5 * AngtoAU
!  real(8), parameter :: height = 0.003

  ! r(xyz,atom,bead)
  r(1,:,:) =  x(:,:)
  r(2,:,:) =  y(:,:)
  r(3,:,:) =  z(:,:)

  ! +++ Calculating Forcde which atom (i) feels from atom (j) +++
  f(:,:,:) = 0.0d0
  Eenergy(:) = 0.0d0
  do imode = 1, Nbead
    r12(:) = r(:,1,imode) - r(:,2,imode)
    r32(:) = r(:,3,imode) - r(:,2,imode)
    dis12  = dsqrt( dot_product(r12(:),r12(:)) )
    dis32  = dsqrt( dot_product(r32(:),r32(:)) )

    f12(:) = - cons * (dis12-req) * r12(:) / dis12
    f32(:) = - cons * (dis32-req) * r32(:) / dis32

    f(:,1,imode) = f(:,1,imode) + f12(:)
    f(:,2,imode) = f(:,2,imode) - f12(:) - f32(:)
    f(:,3,imode) = f(:,3,imode) + f32(:)

!      ! +++ Calculating enetemp +++
    Eenergy(imode) = 0.5d0 * cons * ( (dis12-req)**2 + (dis32-req)**2 )
!    do i = 1, Natom
!      do j = i+1, Natom
!        rij(:) = r(:,i,imode) - r(:,j,imode)
!        dij2 = dot_product(rij(:),rij(:))
!        dij1 = dsqrt(dij2)
!
!        f(:,i,imode) = - 4.0 * height / xwid**4 * (dij2 - xwid**2) * rij(:)
!        f(:,j,imode) = (-1) * f(:,i,imode)
!      ! +++ Calculating enetemp +++
!        Eenergy(imode) = Eenergy(imode) + height / xwid**4 * (dij2 - xwid**2)**2
!
!      end do
!    end do
  end do
  f(:,:,:) = f(:,:,:) * dp_inv
  ! +++ End Calculating Forcde which atom (i) feels from atom (j) +++

  fx(:,:) = f(1,:,:)
  fy(:,:) = f(2,:,:)
  fz(:,:) = f(3,:,:)

  ! +++ Writting output +++
  open(Ucoor,file=trim(address)//'/coor.xyz',status='unknown',form='formatted',position='append')
    write(Ucoor,'(I5)') natom*nbead
    write(Ucoor,'(I10)') istepsv
    do j = 1, Nbead
      do i = 1, Natom
        write(Ucoor,9999) alabel(i), r(:,i,j)*AUtoAng
      end do
    end do
  close(Ucoor)

  if (Lsave_force .eqv. .True.) then
    open(newunit=Ufor,file=trim(address)//'/force.dat',status='unknown',form='formatted',position='append')
      write(Ufor,'("#",I10)') istepsv
      do imode=1,nbead
        do iatom=1,natom
          write(Ufor,9998) f(:,iatom,imode)
        end do
      end do
    close(Ufor)
  end if

  !open(igete,file=trim(address)//'/ene.dat',status='unknown',form='formatted',position='append')
  !  write(igete,*) istepsv
  !  do imode = 1, Nbead
  !    write(igete,*) Eenergy(imode)
  !  end do
  !close(igete)

  potential = 0.0d0
  do imode = 1, Nbead
    potential = potential + Eenergy(imode)
  end do
  potential = potential * dp_inv
  ! +++ End Writting output +++


  9998 format(3E23.15)
  9999 format(a2,1x,E15.9,1x,E15.9,1x,E15.9)
  return
end subroutine Force_DoubleHarmonic
