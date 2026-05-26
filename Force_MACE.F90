module mace_force_config
  use iso_c_binding, only: c_long
  use Parameters, only: lattice, Lperiodic
  use python_mace_interface
  use utility, only: lowerchr, program_abort
  implicit none
  private
  public :: get_mace_paths
  public :: initialize_mace_interface
  public :: finalize_mace_interface
  public :: symbol_to_atomic_number

  character(len=256), save :: mace_python_dir = './MACE'
  character(len=256), save :: mace_model_path = './MACE/small_mace_model.model'
  logical, save :: mace_ready = .false.

contains

  subroutine get_mace_paths(python_dir, model_path)
    character(len=*), intent(out) :: python_dir
    character(len=*), intent(out) :: model_path
    integer :: env_len, env_status

    call get_environment_variable('MACE_PYTHON_DIR', mace_python_dir, &
         length=env_len, status=env_status)
    if (env_status /= 0) mace_python_dir = './MACE'

    call get_environment_variable('MACE_MODEL', mace_model_path, &
         length=env_len, status=env_status)
    if (env_status /= 0) mace_model_path = './MACE/small_mace_model.model'

    python_dir = trim(mace_python_dir)
    model_path = trim(mace_model_path)
  end subroutine get_mace_paths

  subroutine initialize_mace_interface()
    character(len=256) :: python_dir, model_path

    if (mace_ready) return

    call get_mace_paths(python_dir, model_path)
    call mace_initialize(trim(python_dir))
    mace_ready = .true.
  end subroutine initialize_mace_interface

  subroutine finalize_mace_interface()
    if (.not. mace_ready) return

    call mace_finalize()
    mace_ready = .false.
  end subroutine finalize_mace_interface

  integer(c_long) function symbol_to_atomic_number(symbol) result(z)
    character(len=*), intent(in) :: symbol
    character(len=len(symbol)) :: lower_symbol

    lower_symbol = lowerchr(trim(symbol))

    select case(trim(lower_symbol))
      case('h', 'd', 'mu')
        z = 1_c_long
      case('he')
        z = 2_c_long
      case('li')
        z = 3_c_long
      case('be')
        z = 4_c_long
      case('b')
        z = 5_c_long
      case('c')
        z = 6_c_long
      case('n')
        z = 7_c_long
      case('o')
        z = 8_c_long
      case('f')
        z = 9_c_long
      case('ne')
        z = 10_c_long
      case('na')
        z = 11_c_long
      case('mg')
        z = 12_c_long
      case('al')
        z = 13_c_long
      case('si')
        z = 14_c_long
      case('p')
        z = 15_c_long
      case('s')
        z = 16_c_long
      case('cl')
        z = 17_c_long
      case('ar')
        z = 18_c_long
      case('k')
        z = 19_c_long
      case('ca')
        z = 20_c_long
      case('sc')
        z = 21_c_long
      case('ti')
        z = 22_c_long
      case('v')
        z = 23_c_long
      case('cr')
        z = 24_c_long
      case('mn')
        z = 25_c_long
      case('fe')
        z = 26_c_long
      case('co')
        z = 27_c_long
      case('ni')
        z = 28_c_long
      case('cu')
        z = 29_c_long
      case('zn')
        z = 30_c_long
      case('ga')
        z = 31_c_long
      case('ge')
        z = 32_c_long
      case('as')
        z = 33_c_long
      case('se')
        z = 34_c_long
      case('br')
        z = 35_c_long
      case('kr')
        z = 36_c_long
      case('rb')
        z = 37_c_long
      case('sr')
        z = 38_c_long
      case('y')
        z = 39_c_long
      case('zr')
        z = 40_c_long
      case('nb')
        z = 41_c_long
      case('mo')
        z = 42_c_long
      case('tc')
        z = 43_c_long
      case('ru')
        z = 44_c_long
      case('rh')
        z = 45_c_long
      case('pd')
        z = 46_c_long
      case('ag')
        z = 47_c_long
      case('cd')
        z = 48_c_long
      case('in')
        z = 49_c_long
      case('sn')
        z = 50_c_long
      case('sb')
        z = 51_c_long
      case('te')
        z = 52_c_long
      case('i')
        z = 53_c_long
      case('xe')
        z = 54_c_long
      case('cs')
        z = 55_c_long
      case('ba')
        z = 56_c_long
      case('la')
        z = 57_c_long
      case('ce')
        z = 58_c_long
      case('pt')
        z = 78_c_long
      case('au')
        z = 79_c_long
      case('hg')
        z = 80_c_long
      case default
        call program_abort('ERROR!!! Unknown element for MACE: '//trim(symbol))
    end select
  end function symbol_to_atomic_number

end module mace_force_config

subroutine Set_MACE
  use Parameters, only: Lperiodic, lattice, Fout, MyRank
  use mace_force_config, only: get_mace_paths, initialize_mace_interface
  use utility, only: program_abort
  implicit none

  character(len=256) :: python_dir, model_path
  character(len=512) :: helper_path
  integer :: unit_out
  integer :: access
  real(8) :: lattice_norm

  call get_mace_paths(python_dir, model_path)
  helper_path = trim(python_dir)//'/mace_function.py'

  if (access(trim(python_dir), ' ') /= 0) then
    call program_abort('ERROR!!! There is no MACE directory: '//trim(python_dir))
  end if
  if (access(trim(helper_path), ' ') /= 0) then
    call program_abort('ERROR!!! There is no MACE Python helper: '//trim(helper_path))
  end if
  if (access(trim(model_path), ' ') /= 0) then
    call program_abort('ERROR!!! There is no MACE model: '//trim(model_path))
  end if

  lattice_norm = sum(abs(lattice(:,:)))
  if (Lperiodic .eqv. .true. .and. lattice_norm < 1.0d-12) then
    call program_abort('ERROR!!! Invalid $lattice for MACE periodic calculation')
  end if

  call initialize_mace_interface()

  if (MyRank == 0) then
    open(newunit=unit_out, file=Fout, status='old', position='append')
      write(unit_out,'(a)') ' +++++ MACE force field +++++'
      write(unit_out,'(a,a)') ' +++++ MACE Python dir ', trim(python_dir)
      write(unit_out,'(a,a)') ' +++++ MACE model      ', trim(model_path)
      write(unit_out,'(a)') ' +++++ MACE lattice (Angstrom)'
      write(unit_out,'(3F16.8)') lattice(1,:)
      write(unit_out,'(3F16.8)') lattice(2,:)
      write(unit_out,'(3F16.8)') lattice(3,:)
      write(unit_out,*)
    close(unit_out)
  end if

  return
end subroutine Set_MACE

subroutine Force_MACE
  use iso_c_binding, only: c_double, c_long
  use Parameters, &
    only: Natom, Nbead, Ista, Iend, r, fr, pot_bead, alabel, &
          AU2Ang, eV2AU, eVAng2AU, dp_inv, Lperiodic, lattice
  use mace_force_config, only: get_mace_paths, initialize_mace_interface, &
                               symbol_to_atomic_number
  use python_mace_interface, only: mace_calculate_energy_and_forces
  implicit none

  integer :: iatom, idir, imode
  integer :: bead_start, bead_end
  integer(c_long) :: atomic_numbers(Natom)
  real(c_double) :: positions(Natom, 3)
  real(c_double) :: cell(3, 3)
  real(c_double) :: energy
  real(c_double) :: forces(Natom, 3)
  character(len=256) :: python_dir, model_path

  call initialize_mace_interface()
  call get_mace_paths(python_dir, model_path)

  do iatom = 1, Natom
    atomic_numbers(iatom) = symbol_to_atomic_number(alabel(iatom))
  end do

  cell(:,:) = lattice(:,:)

  fr(:,:,:) = 0.0d0
  pot_bead(:) = 0.0d0

#ifdef _mpi_
  bead_start = Ista
  bead_end = Iend
#else
  bead_start = 1
  bead_end = Nbead
#endif

  do imode = bead_start, bead_end
    do iatom = 1, Natom
      do idir = 1, 3
        positions(iatom, idir) = r(idir, iatom, imode) * AU2Ang
      end do
    end do

    call mace_calculate_energy_and_forces( &
         trim(model_path), Natom, atomic_numbers, positions, cell, &
         energy, forces, pbc=Lperiodic)

    pot_bead(imode) = energy * eV2AU
    do iatom = 1, Natom
      do idir = 1, 3
        fr(idir, iatom, imode) = forces(iatom, idir) * eVAng2AU * dp_inv
      end do
    end do
  end do

  return
end subroutine Force_MACE

subroutine Finalize_MACE
  use mace_force_config, only: finalize_mace_interface
  implicit none

  call finalize_mace_interface()

  return
end subroutine Finalize_MACE
