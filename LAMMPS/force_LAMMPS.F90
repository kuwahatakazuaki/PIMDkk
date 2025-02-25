subroutine force_LAMMPS
  use Parameters
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_int64_t
  USE Struct_,           ONLY: xyz_, abc_, force_
  USE LAMMPSCalculator_, ONLY: lammps_calculator_, LAMMPSCalculator
  TYPE(xyz_),   ALLOCATABLE :: cartesian_coordinates (:)   ! Atomic coordinates (Ang)
  TYPE(force_), ALLOCATABLE :: forces(:)                   ! Forces acting on atoms (eV/Ang)
  TYPE(lammps_calculator_) :: calculator   ! LAMMPS calculator
  integer :: Iatom, Imode
  INTEGER(KIND=c_int64_t), POINTER :: total_atoms
  integer, allocatable :: lammps_id(:)
  real(8) :: r_temp(3,Natom)

! Set LAMMPS calculator
  calculator = LAMMPSCalculator( lammps_file = 'lammps.in' )

!! Total number of atoms from LAMMPS
!  total_atoms = calculator%lammps%extract_global( 'natoms' )

allocate( cartesian_coordinates(Natom) )
allocate( forces(Natom) )
allocate( lammps_id(Natom) )

! Calculation in LAMMPS
  !do Imode = 1, Nbead
  do Imode = Ista, Iend
    r_temp(:,:) = r(:,:,Imode)
    r_temp(:,:) = r_temp(:,:) * AUtoAng
    do Iatom = 1, Natom
      cartesian_coordinates(Iatom)%x = r_temp(1,Iatom)
      cartesian_coordinates(Iatom)%y = r_temp(2,Iatom)
      cartesian_coordinates(Iatom)%z = r_temp(3,Iatom)
    end do
    call calculator%run( cartesian_coordinates = cartesian_coordinates(:) )

    Eenergy(Imode) = calculator%potential_energy
    do Iatom = 1, Natom
      fr(1,Iatom,Imode) = calculator%forces(Iatom)%x
      fr(2,Iatom,Imode) = calculator%forces(Iatom)%y
      fr(3,Iatom,Imode) = calculator%forces(Iatom)%z
    end do
  end do
  Eenergy(:) = Eenergy(:) * kcalPmol2AU
  fr(:,:,Ista:Iend) = fr(:,:,Ista:Iend) * dp_inv * kcalPmol2AU / AngtoAU

  lammps_id = calculator%atom_id

  CALL calculator%close()
  DEALLOCATE( cartesian_coordinates )
  DEALLOCATE( forces )

end subroutine force_LAMMPS

