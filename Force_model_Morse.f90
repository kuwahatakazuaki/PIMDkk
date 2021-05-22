! Obtaining force(fx) from position(x)
subroutine Force_model_Morse
  Use Parameters, & 
    only: x, y, z, fx, fy, fz, Natom, Nbead, imode, iatom, Eenergy, potential, &
          alabel, dp_inv, address, igete, igetx, igetf, istepsv, &
          Save_force, &
          AUtoAng =>  bohr_inv
  Use Parameter_tk
implicit none
integer :: i, j!, imode
real(8) :: f_two(3), power, dis !, Epoten
real(8) :: r(3,Natom,Nbead), f(3,Natom,Nbead)
real(8) :: rij(3), temp

!! +++ Constants for conversion +++
!real(8), parameter :: eVtoAU    = 1.0d0/27.21162
!real(8), parameter :: AngtoAU = 1/0.529177249d0
!real(8), parameter :: AUtoAng = 0.529177249d0
!
!! +++ Parameters for Morse Potential +++
!real(8), parameter :: De   = 4.519d0  * eVtoAU   ! eV
!real(8), parameter :: r0   = 0.74d0   * AngtoAU  ! angstrom
!real(8), parameter :: width = 1.981d0 * AUtoAng  ! 1/Angstrom
!! print *, "De and r0", De, r0

! Parameters in Atomic unit
! Hydrogen molecule
  real(8), parameter :: r0    = 1.41014d0
  real(8), parameter :: De    = 0.1745d0
  real(8), parameter :: width = 1.0213d0


! r(xyz,atom,bead)
r(1,:,:) =  x(:,:)
r(2,:,:) =  y(:,:)
r(3,:,:) =  z(:,:)

! +++ Calculating Forcde which atom (i) feels from atom (j) +++
f(:,:,:) = 0.0d0
do imode = 1, Nbead
  do i = 1, Natom
    do j = i+1, Natom
      rij(:) = r(:,i,imode)-r(:,j,imode)
      dis = dsqrt( dot_product( rij(:),rij(:) ) )
      power = -width * (dis-r0)
      f_two(:) = 2 * width * De * exp(power) * (exp(power) - 1) * (rij(:))/dis
      f(:,i,imode) = f(:,i,imode) + f_two(:)
      f(:,j,imode) = f(:,j,imode) - f_two(:)
    end do
  end do
end do
f(:,:,:) = f(:,:,:) * dp_inv
! +++ End Calculating Forcde which atom (i) feels from atom (j) +++

fx(:,:) = f(1,:,:)
fy(:,:) = f(2,:,:)
fz(:,:) = f(3,:,:)

! +++ Calculating enetemp +++
open(20,file=trim(address)//'/distance.dat',status='unknown',position='append')
Eenergy(:) = 0.0d0
write(20,*) istepsv
do imode = 1, Nbead
  do i = 1, Natom
    do j = i+1, Natom
      rij(:) = r(:,i,imode)-r(:,j,imode)
      dis = dsqrt( dot_product( rij(:),rij(:) ) )
      temp = 1 - exp(-width * (dis-r0))
      Eenergy(imode) = Eenergy(imode) + De * temp * temp
      write(20,*) dis
    end do
  end do
end do
close(20)
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

open(igete,file=trim(address)//'/ene.dat',status='unknown',form='formatted',position='append')
  write(igete,*) istepsv
  do imode = 1, Nbead
    write(igete,*) Eenergy(imode)
  end do
close(igete)

potential = 0.0d0
do imode = 1, Nbead
  potential = potential + Eenergy(imode)
end do
potential = potential * dp_inv
! +++ End Writting output +++


9998 format(3E23.15)
9999 format(a2,1x,E15.9,1x,E15.9,1x,E15.9)
return
end subroutine Force_model_Morse


