subroutine calc_internal_pressure
  use Parameters
  use utility, only: program_abort
  implicit none
  integer :: imode, iatom
  real(8) :: K_c, W_cv, E_cv, qkinetic_eff

  if ( .not. allocated(W_pot_bead) ) then
    call program_abort("W_pot_bead is not allocated")
  end if

  if ( vol <= 0.0d0 ) then
    call program_abort("calc_internal_pressure requires positive cell volume")
  end if

  W_pot = sum(W_pot_bead(:)) * dp_inv

  K_c = 0.0d0
  do iatom = 1, Natom
    K_c = K_c + 0.5d0 * physmass(iatom) * dot_product(vur(:,iatom,1), vur(:,iatom,1))
  end do

  W_cv = 0.0d0
  do imode = 1, Nbead
    do iatom = 1, Natom
      W_cv = W_cv + dot_product(fr(:,iatom,imode), r(:,iatom,imode) - ur(:,iatom,1))
    end do
  end do
  E_cv = 1.5d0*dble(Natom)/beta - 0.5d0*W_cv

  if ( Nbead == 1 ) then
    qkinetic_eff = 0.0d0
  else
    qkinetic_eff = qkinetic
  end if

  press_inst = (2.0d0*K_c - W_cv + W_pot) / (3.0d0*vol)
  press_cv   = (2.0d0*E_cv + W_pot) / (3.0d0*vol)
  press_prim = (3.0d0*dble(Natom*Nbead)/beta - 2.0d0*qkinetic_eff + W_pot) / (3.0d0*vol)

  return
end subroutine calc_internal_pressure
