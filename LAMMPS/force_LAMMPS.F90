subroutine force_LAMMPS
  use Parameters
  use utility
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_int64_t
  USE Struct_,           ONLY: xyz_, abc_, force_
  USE LAMMPSCalculator_, ONLY: lammps_calculator_, LAMMPSCalculator
  TYPE(xyz_),   ALLOCATABLE :: cartesian_coordinates (:)   ! Atomic coordinates (Ang)
  TYPE(force_), ALLOCATABLE :: forces(:)                   ! Forces acting on atoms (eV/Ang)
  TYPE(lammps_calculator_) :: calculator   ! LAMMPS calculator
  integer :: Iatom, Imode
  integer :: i, Uout
  INTEGER(KIND=c_int64_t), POINTER :: total_atoms
  integer, allocatable :: lammps_id(:)
  real(8) :: r_temp(3,Natom), lattice(3,3), lat_inv(3,3)
  real(8), allocatable :: s(:,:)


lattice(1,:) = [12.8512848674,    0.0000000000,    0.0000000000]
lattice(2,:) = [ 0.2132032072,   12.6594123484,    0.0000000000]
lattice(3,:) = [ 0.1438686986,   -0.0357082118,   12.7223433492]
call get_inv_mat(lattice,lat_inv,3)

! Set LAMMPS calculator
  calculator = LAMMPSCalculator( lammps_file = 'lammps.in' )

!! Total number of atoms from LAMMPS
!  total_atoms = calculator%lammps%extract_global( 'natoms' )

  allocate( cartesian_coordinates(Natom) )
  allocate( forces(Natom) )
  allocate( lammps_id(Natom) )
allocate( s(3,Natom) )

! Calculation in LAMMPS
  !do Imode = 1, Nbead
  do Imode = Ista, Iend
    r_temp(:,:) = r(:,:,Imode)
    r_temp(:,:) = r_temp(:,:) * AUtoAng
    do Iatom = 1, Natom
s(:,Iatom) = matmul(r_temp(:,Iatom),lat_inv)
s(:,Iatom) = s(:,Iatom) - dble(nint(s(:,Iatom)+0.5d0))+1.0d0
r_temp(:,Iatom) = matmul(s(:,Iatom),lattice)
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

!print *, "# step", istepsv
!lammps_id = calculator%atom_id
!open(newunit=Uout,file='lammps_id.out',position='append')
!  write(Uout,*) "# step", istepsv
!  do i = 1, Natom
!    write(Uout,*) lammps_id(i)
!  end do
!close(Uout)
!open(newunit=Uout,file='lammps_force.out',position='append')
!  write(Uout,*) istepsv
!  do i = 1, Natom
!    write(Uout,*) fr(:,i,1)
!  end do
!close(Uout)

  CALL calculator%close()
  DEALLOCATE( cartesian_coordinates )
  DEALLOCATE( forces )

contains

  subroutine get_inv_mat(mat,inv,n)
    integer :: n
    real(8), intent(in)  :: mat(n,n)
    real(8), intent(out) :: inv(n,n)
    integer :: lwork, lda, info
    real(8), allocatable :: work(:)
    integer, allocatable :: ipiv(:)
    inv(:,:) = mat(:,:)
    lda = n
    lwork = 64*n
    allocate(work(lwork),ipiv(n))
    call dgetrf(N, N, inv, lda, ipiv, info)
    call dgetri(N, inv, lda, ipiv, work, lwork, info)
  end subroutine get_inv_mat

end subroutine force_LAMMPS

