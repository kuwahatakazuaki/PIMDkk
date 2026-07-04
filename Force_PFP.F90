module pfp_force_config
  use iso_c_binding, only: c_long
  use python_pfp_interface
  use utility, only: program_abort
  use element_table, only: symbol_to_atomic_number
  implicit none
  private
  public :: resolve_pfp_paths
  public :: get_pfp_calc_mode
  public :: initialize_pfp_interface
  public :: finalize_pfp_interface
  public :: symbol_to_atomic_number
  public :: pfp_python_dir
  public :: pfp_helper_path

  character(len=1024), save :: pfp_python_dir = './PFP'
  character(len=1024), save :: pfp_helper_path = './PFP/pfp_function.py'
  logical, save :: pfp_ready = .false.
  logical, save :: pfp_paths_ready = .false.

contains

  subroutine resolve_pfp_paths()
    integer :: env_len, env_status

    if (pfp_paths_ready) return

    call get_environment_variable('PFP_PYTHON_DIR', pfp_python_dir, &
         length=env_len, status=env_status)
    if (env_status == -1) then
      call program_abort('ERROR!!! PFP_PYTHON_DIR is too long')
    end if
    if (env_status /= 0 .or. env_len <= 0 .or. len_trim(pfp_python_dir) == 0) then
      call program_abort('ERROR!!! PFP_PYTHON_DIR is not set'//char(10)// &
                         'Please set PFP_PYTHON_DIR, e.g. '// &
                         'export PFP_PYTHON_DIR="/path/to/PFP"')
    end if

    pfp_helper_path = trim(pfp_python_dir)//'/pfp_function.py'
    pfp_paths_ready = .true.
  end subroutine resolve_pfp_paths

  function get_pfp_calc_mode() result(calc_mode)
    use Parameters, only: Lperiodic
    character(len=16) :: calc_mode

    ! NOTE: $PFP_calc_mode による明示上書きは Stage 3b (Parameter.F90/
    ! read_input.F90 へのフィールド追加) で配線する (D5)。現段階では Lperiodic
    ! からの自動選択のみ。
    if (Lperiodic) then
      calc_mode = 'CRYSTAL'
    else
      calc_mode = 'MOLECULE'
    end if
  end function get_pfp_calc_mode

  subroutine initialize_pfp_interface()
    if (pfp_ready) return

    call resolve_pfp_paths()
    call pfp_initialize(pfp_python_dir)
    pfp_ready = .true.
  end subroutine initialize_pfp_interface

  subroutine finalize_pfp_interface()
    if (.not. pfp_ready) return

    call pfp_finalize()
    pfp_ready = .false.
  end subroutine finalize_pfp_interface

end module pfp_force_config

subroutine set_PFP
  use Parameters, only: Lperiodic, lattice
  use pfp_force_config, only: pfp_python_dir, pfp_helper_path, &
                              resolve_pfp_paths, initialize_pfp_interface
  use utility, only: program_abort
  implicit none

  integer :: access
  real(8) :: lattice_norm

  call resolve_pfp_paths()

  if (access(pfp_python_dir, ' ') /= 0) then
    call program_abort('ERROR!!! There is no PFP directory: '//trim(pfp_python_dir)// &
                       char(10)//'Please set PFP_PYTHON_DIR, e.g. '// &
                       'export PFP_PYTHON_DIR="/path/to/PFP"')
  end if
  if (access(pfp_helper_path, ' ') /= 0) then
    call program_abort('ERROR!!! There is no PFP Python helper: '//trim(pfp_helper_path)// &
                       char(10)//'PFP_PYTHON_DIR must contain pfp_function.py')
  end if

  lattice_norm = sum(abs(lattice(:,:)))
  if (Lperiodic .eqv. .true. .and. lattice_norm < 1.0d-12) then
    call program_abort('ERROR!!! Invalid $lattice for PFP periodic calculation')
  end if

  call initialize_pfp_interface()

  return
end subroutine set_PFP

subroutine Force_PFP
  use iso_c_binding, only: c_double, c_long
  use Parameters, &
    only: Natom, Ista, Iend, r, fr, pot_bead, alabel, &
          AU2Ang, eV2AU, eVAng2AU, dp_inv, Lperiodic, lattice
  use pfp_force_config, only: get_pfp_calc_mode, initialize_pfp_interface, &
                              symbol_to_atomic_number
  use python_pfp_interface, only: pfp_calculate_energy_and_forces
  implicit none

  integer :: iatom, idir, imode
  integer(c_long) :: atomic_numbers(Natom)
  real(c_double) :: positions(Natom, 3)
  real(c_double) :: cell(3, 3)
  real(c_double) :: energy
  real(c_double) :: forces(Natom, 3)
  character(len=16) :: calc_mode

  call initialize_pfp_interface()
  calc_mode = get_pfp_calc_mode()

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

    call pfp_calculate_energy_and_forces( &
         Natom, atomic_numbers, positions, cell, energy, forces, &
         pbc=Lperiodic, calc_mode=trim(calc_mode))

    pot_bead(imode) = energy * eV2AU
    do iatom = 1, Natom
      do idir = 1, 3
        fr(idir, iatom, imode) = forces(iatom, idir) * eVAng2AU * dp_inv
      end do
    end do
  end do

  return
end subroutine Force_PFP

subroutine Finalize_PFP
  use pfp_force_config, only: finalize_pfp_interface
  implicit none
  call finalize_pfp_interface()
  return
end subroutine Finalize_PFP
