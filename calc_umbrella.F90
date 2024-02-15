! Obtaining force(fx) from position(x)
subroutine calc_umbrella
  Use Parameters, &
    only: r, fr, Natom, Nbead,  Eenergy, potential, myrank, &
    !only: x, y, z, fx, fy, fz, Natom, Nbead,  Eenergy, potential, myrank, &
          alabel, dp_inv, address, istepsv, Iforce, &
          AUtoAng => bohr_inv, KtoAU => boltz, AngtoAU => bohr, &
          atom1 => umbrella_atom1, atom2  => umbrella_atom2, &
          umbrella_width, umbrella_height, ista, iend, &
          umbrella_constant, umbrella_sampling
          ! cons => umbrella_constant, umbrella_sampling
  !use MPI
implicit none
integer :: i, j, imode, iatom
real(8) :: r12(3), f12(3), d12
!real(8) :: r(3,Natom,Nbead), f(3,Natom,Nbead)
real(8) :: cons

! cons = cons / AngtoAU
! cons = umbrella_constant / (AngtoAU*AngtoAU)
cons = umbrella_constant * AUtoAng * AUtoAng
! cons = umbrella_constant

! r(xyz,atom,bead)
!r(1,:,:) =  x(:,:)
!r(2,:,:) =  y(:,:)
!r(3,:,:) =  z(:,:)

!f(1,:,:) = fx(:,:)
!f(2,:,:) = fy(:,:)
!f(3,:,:) = fz(:,:)

! +++ Calculating Force +++
!! +++ Specific in Double well potential +++ 
!if (Iforce == 12) then
!  block
!    real(8) :: a, distance, half_dis
!    real(8) :: width, height
!    width  = umbrella_width  * AngtoAU
!    height = umbrella_height * KtoAU
!    a = height / (width * width)
!    do imode = 1, Nbead
!
!      do i = 1, Natom
!        f12(1) = (-2) * a * r(1,i,imode) * dp_inv
!        f12(2) = 0.0d0
!        f12(3) = 0.0d0
!
!        f(:,i,imode) = f(:,i,imode) + f12(:)
!      end do
!    ! +++ End Specific in Double well potential +++
!    end do
!  end block
!end if

if ( umbrella_sampling == 1 ) then
  do imode = ista, iend
    r12(:) = r(:,atom1,imode) - r(:,atom2,imode)
    f12(:) = (-2) * cons * r12(:) * dp_inv

    fr(:,atom1,imode) = fr(:,atom1,imode) + f12(:)
    fr(:,atom2,imode) = fr(:,atom2,imode) - f12(:)
  end do
end if

!fx(:,:) = f(1,:,:)
!fy(:,:) = f(2,:,:)
!fz(:,:) = f(3,:,:)

Call Start_Send_Recv_MPI_tk


return
end subroutine calc_umbrella



