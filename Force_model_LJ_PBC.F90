subroutine Force_model_LJ_PBC
  use Parameters, &
    only: r, fr, Natom, Nbead, pot_bead, W_pot_bead, dp_inv, &
          lattice, Lperiodic, Ndim, Ista, Iend, K2AU, Ang2AU
  use utility, only: program_abort, calc_inv_mat33
  implicit none
  integer :: imode, iatom, jatom
  real(8), parameter :: epsilon_K = 119.8d0
  real(8), parameter :: sigma_Ang = 3.405d0
  real(8), parameter :: min_r2 = 1.0d-20
  real(8) :: epsilon, sigma, cutoff, cutoff2
  real(8) :: cell(3,3), invcell(3,3)
  real(8) :: avec(3), bvec(3), cvec(3)
  real(8) :: rij(3), sij(3), fij(3)
  real(8) :: r2, inv_r2, sr2, sr6, sr12, force_factor

  if ( .not. Lperiodic ) then
    call program_abort('Force_model_LJ_PBC requires Lperiodic = .T.')
  end if
  if ( Ndim /= 3 ) then
    call program_abort('Force_model_LJ_PBC requires Ndim = 3')
  end if

  epsilon = epsilon_K * K2AU
  sigma = sigma_Ang * Ang2AU

  cell(:,1) = lattice(1,:) * Ang2AU
  cell(:,2) = lattice(2,:) * Ang2AU
  cell(:,3) = lattice(3,:) * Ang2AU
  invcell(:,:) = calc_inv_mat33(cell)

  avec(:) = cell(:,1)
  bvec(:) = cell(:,2)
  cvec(:) = cell(:,3)
  cutoff = 0.5d0 * min(sqrt(dot_product(avec,avec)), &
                       sqrt(dot_product(bvec,bvec)), &
                       sqrt(dot_product(cvec,cvec)))
  if ( cutoff <= 0.0d0 ) then
    call program_abort('Force_model_LJ_PBC requires positive lattice vectors')
  end if
  cutoff2 = cutoff * cutoff

  fr(:,:,:) = 0.0d0
  pot_bead(:) = 0.0d0
  W_pot_bead(:) = 0.0d0

  do imode = Ista, Iend
    do iatom = 1, Natom - 1
      do jatom = iatom + 1, Natom
        rij(:) = r(:,iatom,imode) - r(:,jatom,imode)
        sij(:) = matmul(invcell, rij)
        sij(:) = sij(:) - dnint(sij(:))
        rij(:) = matmul(cell, sij)
        r2 = dot_product(rij, rij)

        if ( r2 <= min_r2 ) then
          call program_abort('Zero pair distance in Force_model_LJ_PBC')
        end if

        if ( r2 < cutoff2 ) then
          inv_r2 = 1.0d0 / r2
          sr2 = sigma * sigma * inv_r2
          sr6 = sr2 * sr2 * sr2
          sr12 = sr6 * sr6
          pot_bead(imode) = pot_bead(imode) + 4.0d0 * epsilon * (sr12 - sr6)
          force_factor = 24.0d0 * epsilon * (2.0d0*sr12 - sr6) * inv_r2
          fij(:) = force_factor * rij(:)
          fr(:,iatom,imode) = fr(:,iatom,imode) + fij(:)
          fr(:,jatom,imode) = fr(:,jatom,imode) - fij(:)
          W_pot_bead(imode) = W_pot_bead(imode) + dot_product(fij, rij)
        end if
      end do
    end do
  end do

  fr(:,:,Ista:Iend) = fr(:,:,Ista:Iend) * dp_inv

  return
end subroutine Force_model_LJ_PBC
