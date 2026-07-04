module mace_force_config
  use iso_c_binding, only: c_long
  use python_mace_interface
  use utility, only: program_abort
  use element_table, only: symbol_to_atomic_number
  implicit none
  private
  public :: get_mace_model_path
  public :: initialize_mace_interface
  public :: finalize_mace_interface
  public :: symbol_to_atomic_number
  public :: mace_python_dir
  public :: mace_helper_path

  character(len=1024), save :: mace_python_dir = './MACE'
  character(len=1024), save :: mace_helper_path = './MACE/mace_function.py'
  character(len=1024), save :: mace_model_path = './MACE/small_mace_model.model'
  logical, save :: mace_ready = .false.
  logical, save :: mace_paths_ready = .false.

contains

  subroutine resolve_mace_paths()
    integer :: env_len, env_status

    if (mace_paths_ready) return

    call get_environment_variable('MACE_PYTHON_DIR', mace_python_dir, &
         length=env_len, status=env_status)
    if (env_status == -1) then
      call program_abort('ERROR!!! MACE_PYTHON_DIR is too long')
    end if
    if (env_status /= 0 .or. env_len <= 0 .or. len_trim(mace_python_dir) == 0) then
      call program_abort('ERROR!!! MACE_PYTHON_DIR is not set'//char(10)// &
                         'Please set MACE_PYTHON_DIR, e.g. '// &
                         'export MACE_PYTHON_DIR="/path/to/MACE"')
    end if

    mace_helper_path = trim(mace_python_dir)//'/mace_function.py'
    mace_paths_ready = .true.
  end subroutine resolve_mace_paths

  subroutine get_mace_model_path(model_path)
    character(len=*), intent(out) :: model_path
    integer :: env_len, env_status

    call resolve_mace_paths()

    call get_environment_variable('MACE_MODEL', mace_model_path, length=env_len, status=env_status)
    if (env_status == -1) then
      call program_abort('ERROR!!! MACE_MODEL is too long')
    end if
    if (env_status /= 0 .or. env_len <= 0 .or. len_trim(mace_model_path) == 0) then
      call program_abort('ERROR!!! MACE_MODEL is not set'//char(10)// &
                         'Please set MACE_MODEL, e.g. '// &
                         'export MACE_MODEL="/path/to/your/model.model"')
    end if

    model_path = trim(mace_model_path)
  end subroutine get_mace_model_path

  subroutine initialize_mace_interface()
    if (mace_ready) return

    call resolve_mace_paths()
    call mace_initialize(mace_python_dir)
    mace_ready = .true.
  end subroutine initialize_mace_interface

  subroutine finalize_mace_interface()
    if (.not. mace_ready) return

    call mace_finalize()
    mace_ready = .false.
  end subroutine finalize_mace_interface

end module mace_force_config

subroutine set_MACE
  use Parameters, only: Lperiodic, lattice, Fout, MyRank
  use mace_force_config, only: mace_python_dir, mace_helper_path, &
                               get_mace_model_path, initialize_mace_interface
  use utility, only: program_abort
  implicit none

  character(len=256) :: model_path
  integer :: access
  real(8) :: lattice_norm

  call get_mace_model_path(model_path)

  if (access(mace_python_dir, ' ') /= 0) then
    call program_abort('ERROR!!! There is no MACE directory: '//trim(mace_python_dir)// &
                       char(10)//'Please set MACE_PYTHON_DIR, e.g. '// &
                       'export MACE_PYTHON_DIR="/path/to/MACE"')
  end if
  if (access(mace_helper_path, ' ') /= 0) then
    call program_abort('ERROR!!! There is no MACE Python helper: '//trim(mace_helper_path)// &
                       char(10)//'MACE_PYTHON_DIR must contain mace_function.py')
  end if
  if (access(trim(model_path), ' ') /= 0) then
    call program_abort('ERROR!!! There is no MACE model: '//trim(model_path)// &
                       char(10)//'Please set MACE_MODEL, e.g. '// &
                       'export MACE_MODEL="/path/to/your/model.model"')
  end if

  lattice_norm = sum(abs(lattice(:,:)))
  if (Lperiodic .eqv. .true. .and. lattice_norm < 1.0d-12) then
    call program_abort('ERROR!!! Invalid $lattice for MACE periodic calculation')
  end if

  call initialize_mace_interface()

  return
end subroutine Set_MACE

subroutine Force_MACE
  use iso_c_binding, only: c_double, c_long
  use Parameters, &
    only: Natom, Nbead, Ista, Iend, r, fr, pot_bead, alabel, &
          AU2Ang, eV2AU, eVAng2AU, dp_inv, Lperiodic, lattice, &
          device
  use mace_force_config, only: get_mace_model_path, &
                               initialize_mace_interface, symbol_to_atomic_number
  use python_mace_interface, only: mace_calculate_energy_and_forces
  implicit none

  integer :: iatom, idir, imode
  integer(c_long) :: atomic_numbers(Natom)
  real(c_double) :: positions(Natom, 3)
  real(c_double) :: cell(3, 3)
  real(c_double) :: energy
  real(c_double) :: forces(Natom, 3)
  character(len=256) :: model_path

  call initialize_mace_interface()
  call get_mace_model_path(model_path)

  do iatom = 1, Natom
    atomic_numbers(iatom) = symbol_to_atomic_number(alabel(iatom))
  end do

  cell(:,:) = lattice(:,:)

  fr(:,:,:) = 0.0d0
  pot_bead(:) = 0.0d0

  do imode = Ista, Iend
    do iatom = 1, Natom
      do idir = 1, 3
        positions(iatom, idir) = r(idir, iatom, imode) * AU2Ang
      end do
    end do

    call mace_calculate_energy_and_forces( &
         trim(model_path), Natom, atomic_numbers, positions, cell, &
         energy, forces, pbc=Lperiodic, device=trim(device))

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
