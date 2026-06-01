subroutine add_constrain
  use Parameters, &
    only: r, fur_ref, tnm, Ndim, Natom, Nbead, dp_inv, Ang2AU, &
          Icons, cons_strength, cons_val, potential_cons, cons_cv_ave, &
          dVdcons_cv, atom1 => cons_atom1, atom2 => cons_atom2
  use utility, only: program_abort
  implicit none
  integer :: imode, iatom
  real(8) :: dr, r0_cons, k_cons
  real(8) :: rij(Ndim), fij(Ndim)
  real(8), allocatable :: dist(:)
  real(8), allocatable :: fcart_cons(:,:,:)
  real(8), parameter :: tiny_dist = 1.0d-12

  potential_cons = 0.0d0
  cons_cv_ave = 0.0d0
  dVdcons_cv = 0.0d0
  if ( Icons == 0 ) return

  call nmtrans_ur2r

  select case(Icons)
    case(1)
      allocate( dist(Nbead) )
      allocate( fcart_cons(Ndim,Natom,Nbead) )
      fcart_cons(:,:,:) = 0.0d0
      r0_cons = cons_val * Ang2AU
      k_cons = cons_strength

      do imode = 1, Nbead
        rij(:) = r(:,atom1,imode) - r(:,atom2,imode)
        dist(imode) = norm2(rij(:))
        if ( dist(imode) < tiny_dist ) then
          call program_abort('ERROR!!: constrained atom distance is too small')
        end if
        cons_cv_ave = cons_cv_ave + dist(imode)
      end do
      cons_cv_ave = cons_cv_ave * dp_inv

      dr = cons_cv_ave - r0_cons
      potential_cons = 0.5d0 * k_cons * dr * dr
      dVdcons_cv = k_cons * dr

      do imode = 1, Nbead
        rij(:) = r(:,atom1,imode) - r(:,atom2,imode)
        fij(:) = -dVdcons_cv * dp_inv * rij(:) / dist(imode)
        fcart_cons(:,atom1,imode) = fcart_cons(:,atom1,imode) + fij(:)
        fcart_cons(:,atom2,imode) = fcart_cons(:,atom2,imode) - fij(:)
      end do

      do iatom = 1, Natom
        fur_ref(:,iatom,:) = fur_ref(:,iatom,:) + matmul(fcart_cons(:,iatom,:), tnm)
      end do
      deallocate( dist )
      deallocate( fcart_cons )
    case default
      call program_abort('ERROR!!: unsupported Icons option in add_constrain')
  end select

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
