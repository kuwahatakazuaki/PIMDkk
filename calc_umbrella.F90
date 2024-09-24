! Obtaining force(fx) from position(x)
subroutine calc_umbrella
  use Parameters, &
    only: r, fr, Natom, Nbead,  Eenergy, potential, myrank, &
          alabel, dp_inv, dir_result, istepsv, Iforce, &
          AUtoAng, KtoAU, AngtoAU, &
          atom1 => umb_atom1, atom2  => umb_atom2, atom3 => umb_atom3, &
          Ista, Iend, Iumb, umb_cons, umb_pot
  implicit none
  integer :: i, j, Imode, Iatom
  real(8) :: rij(3), fij(3), dij
  real(8) :: r1(3), r2(3), r3(3), r_mid(3)

  select case(Iumb)
    case(1)
      do Imode = Ista, Iend
        rij(:) = r(:,atom1,Imode) - r(:,atom2,Imode)
        fij(:) = (-2) * umb_cons * rij(:) * dp_inv

        fr(:,atom1,Imode) = fr(:,atom1,Imode) + fij(:)
        fr(:,atom2,Imode) = fr(:,atom2,Imode) - fij(:)
      end do
    case(2)
      do Imode = Ista, Iend
        r1(:) = r(:,atom1,Imode)
        r2(:) = r(:,atom2,Imode)
        r3(:) = r(:,atom3,Imode)

        r_mid(:) = 0.5d0*(r1(:)+r2(:))
        rij(:)  = r3(:) -  r_mid(:)
        umb_pot = umb_cons * dot_product(rij(:),rij(:))
        Eenergy(Imode) = Eenergy(Imode) + umb_pot

        fij(:) = umb_cons * rij(:) * dp_inv
        fr(:,atom1,Imode) = fr(:,atom1,Imode) + fij(:)
        fr(:,atom2,Imode) = fr(:,atom2,Imode) + fij(:)
        fr(:,atom3,Imode) = fr(:,atom3,Imode) + fij(:) * (-2.0d0)
      end do
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

