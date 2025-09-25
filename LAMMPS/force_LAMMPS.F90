subroutine force_LAMMPS
  use Parameters
  !USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_int64_t, c_int32_t
  USE Struct_,           ONLY: xyz_, force_
  USE LAMMPSCalculator_, ONLY: lammps_calculator_, LAMMPSCalculator!, atom_id
  USE LIBLAMMPS
  TYPE(xyz_),   ALLOCATABLE :: cartesian_coordinates (:)   ! Atomic coordinates (Ang)
  TYPE(force_), ALLOCATABLE :: forces(:)                   ! Forces acting on atoms (eV/Ang)
  TYPE(lammps_calculator_) :: calculator   ! LAMMPS calculator
  integer :: Iatom, Imode
  integer :: i, Uout
  !INTEGER(KIND=c_int64_t), POINTER :: total_atoms
  real(8) :: r_temp(3,Natom)


! Set LAMMPS calculator
  calculator = LAMMPSCalculator( lammps_file = 'lammps.in' )


!! Total number of atoms from LAMMPS
  !total_atoms = calculator%lammps%extract_global( 'natoms' )

  !atom_id = calculator%lammps%extract_atom( 'id' )

  allocate( cartesian_coordinates(Natom) )
  allocate( forces(Natom) )

! Calculation in LAMMPS
  do Imode = Ista, Iend
    r_temp(:,:) = r(:,:,Imode) * AU2Ang
    do Iatom = 1, Natom
      cartesian_coordinates(Iatom)%x = r_temp(1,Iatom)
      cartesian_coordinates(Iatom)%y = r_temp(2,Iatom)
      cartesian_coordinates(Iatom)%z = r_temp(3,Iatom)
    end do
    call calculator%run( cartesian_coordinates = cartesian_coordinates(:) )

    pot_bead(Imode) = calculator%potential_energy
    do Iatom = 1, Natom
      fr(1,Iatom,Imode) = calculator%forces( Iatom )%x
      fr(2,Iatom,Imode) = calculator%forces( Iatom )%y
      fr(3,Iatom,Imode) = calculator%forces( Iatom )%z
    end do
  end do
  !pot_bead(:) = pot_bead(:) * kcalPmol2AU
  !fr(:,:,Ista:Iend) = fr(:,:,Ista:Iend) * dp_inv * kcalPmol2AU / Ang2AU
  pot_bead(:) = pot_bead(:) * eV2AU
  fr(:,:,Ista:Iend) = fr(:,:,Ista:Iend) * dp_inv * (eV2AU/Ang2AU)

  CALL calculator%close()
  DEALLOCATE( cartesian_coordinates )
  DEALLOCATE( forces )

end subroutine force_LAMMPS

