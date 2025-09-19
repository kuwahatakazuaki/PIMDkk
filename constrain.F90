subroutine  add_constrain
  use Parameters, &
    only: r, fr, ur, Natom, Nbead,  pot_bead, potential, myrank, &
          alabel, dp_inv, dir_result, istepsv, Iforce, &
          AU2Ang, K2AU, Ang2AU, &
          Icons, cons_strenght, cons_val, &
          atom1 => cons_atom1, atom2  => cons_atom2, atom3 => cons_atom3
  real(8) :: dis2, uij(3)
  uij(:) = ur(:,atom1,1) - ur(:,atom2,1)

end subroutine add_constrain


!subroutine calc_umbrella
!  use Parameters, &
!    only: r, fr, Natom, Nbead,  pot_bead, potential, myrank, &
!          alabel, dp_inv, dir_result, istepsv, Iforce, &
!          AU2Ang, K2AU, Ang2AU, &
!          atom1 => umb_atom1, atom2  => umb_atom2, atom3 => umb_atom3, &
!          Ista, Iend, Iumb, umb_cons, umb_pot
!  implicit none
!  integer :: i, j, Imode, Iatom
!  real(8) :: rij(3), fij(3), dij
!  real(8) :: r1(3), r2(3), r3(3), r_mid(3)
!
!  select case(Iumb)
!    case(1)
!      do Imode = Ista, Iend
!        rij(:) = r(:,atom1,Imode) - r(:,atom2,Imode)
!        fij(:) = (-2) * umb_cons * rij(:) * dp_inv
!
!        fr(:,atom1,Imode) = fr(:,atom1,Imode) + fij(:)
!        fr(:,atom2,Imode) = fr(:,atom2,Imode) - fij(:)
!      end do
!    case(2)
!      do Imode = Ista, Iend
!        r1(:) = r(:,atom1,Imode)
!        r2(:) = r(:,atom2,Imode)
!        r3(:) = r(:,atom3,Imode)
!
!        r_mid(:) = 0.5d0*(r1(:)+r2(:))
!        rij(:)  = r3(:) -  r_mid(:)
!        umb_pot = umb_cons * dot_product(rij(:),rij(:))
!        pot_bead(Imode) = pot_bead(Imode) + umb_pot
!
!        fij(:) = umb_cons * rij(:) * dp_inv
!        fr(:,atom1,Imode) = fr(:,atom1,Imode) + fij(:)
!        fr(:,atom2,Imode) = fr(:,atom2,Imode) + fij(:)
!        fr(:,atom3,Imode) = fr(:,atom3,Imode) + fij(:) * (-2.0d0)
!      end do
!  end select
!
!return
!end subroutine calc_umbrella

