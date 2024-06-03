! Obtaining force(fx) from position(x)
subroutine calc_umbrella
  use Parameters, &
    only: r, fr, Natom, Nbead,  Eenergy, potential, myrank, &
          alabel, dp_inv, address, istepsv, Iforce, &
          AUtoAng, KtoAU, AngtoAU, &
          atom1 => umbrella_atom1, atom2  => umbrella_atom2, atom3 => umbrella_atom3, &
          Ista, Iend, Iumbrella, cons => umbrella_constant
  implicit none
  integer :: i, j, Imode, Iatom
  real(8) :: rij(3), fij(3), dij

  select case(Iumbrella)
    case(1)
      do imode = Ista, Iend
        rij(:) = r(:,atom1,imode) - r(:,atom2,imode)
        fij(:) = (-2) * cons * rij(:) * dp_inv

        fr(:,atom1,imode) = fr(:,atom1,imode) + fij(:)
        fr(:,atom2,imode) = fr(:,atom2,imode) - fij(:)
      end do
    case(2)
  end select

return
end subroutine calc_umbrella


  ! cons = umbrella_constant * AUtoAng * AUtoAng
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

