! Obtaining force(fx) from position(x)
subroutine Force_model_DoubleWell
  Use Parameters, &
    only: x, y, z, fx, fy, fz, Natom, Nbead, imode, iatom, Eenergy, potential, &
          alabel, dp_inv, address, igete, igetx, igetf, istepsv, &
          Save_force, &
          AUtoAng => bohr_inv, KtoAU => boltz, AngtoAU => bohr
  Use Parameter_tk
implicit none
integer :: i, j!, imode
real(8) :: f_two(3), power, dis !, Epoten
real(8) :: r(3,Natom,Nbead), f(3,Natom,Nbead)
real(8) :: rij(3), temp


real(8), parameter :: xwid = 0.5 * AngtoAU
real(8), parameter :: v0   = 1000 * KtoAU


! r(xyz,atom,bead)
r(1,:,:) =  x(:,:)
r(2,:,:) =  y(:,:)
r(3,:,:) =  z(:,:)

! +++ Calculating Forcde which atom (i) feels from atom (j) +++
f(:,:,:) = 0.0d0
do imode = 1, Nbead
  do i = 1, Natom
    f(1,i,imode) = - 4.0 * v0 / xwid**4 * r(1,i,imode) * (r(1,i,imode)**2 - xwid**2)
  end do
end do
f(:,:,:) = f(:,:,:) * dp_inv
! +++ End Calculating Forcde which atom (i) feels from atom (j) +++

fx(:,:) = f(1,:,:)
fy(:,:) = f(2,:,:)
fz(:,:) = f(3,:,:)

! +++ Calculating enetemp +++
Eenergy(:) = 0.0d0
do imode = 1, Nbead
  do i = 1, Natom
    Eenergy(imode) = Eenergy(imode) + v0 / xwid**4 * (r(1,i,imode)**2 - xwid**2)**2
  end do
end do
! +++ End Calculating enetemp +++


! +++ Writting output +++
open(igetx,file=trim(address)//'/coor.xyz',status='unknown',form='formatted',position='append')
  write(igetx,*) Natom * Nbead
  write(igetx,*) istepsv
  do j = 1, Nbead
    do i = 1, Natom
!      write(igetx,9999) alabel(i), x(i,j)*bohr_inv, y(i,j)*bohr_inv,z(i,j)*bohr_inv
      write(igetx,9999) alabel(i), r(:,i,j)*AUtoAng
    end do
  end do
close(igetx)

if (Save_force .eqv. .True.) then
  open(igetf,file=trim(address)//'/force.dat',status='unknown',form='formatted',position='append')
    write(igetf,*) istepsv
    do imode=1,nbead
      do iatom=1,natom
        write(igetf,9998) fx(iatom,imode),fy(iatom,imode),fz(iatom,imode)
!        write(igetf,9998) f(:,iatom,imode)
      end do
    end do
  close(igetf)
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
end subroutine Force_model_DoubleWell


