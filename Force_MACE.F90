module mace_force_config
  use iso_c_binding, only: c_long
  use Parameters, only: lattice, Lperiodic
  use python_mace_interface
  use utility, only: lowerchr, program_abort
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
    character(len=1024) :: exe_path
    character(len=1024) :: exe_dir
    character(len=8192) :: path_env
    character(len=1024) :: path_dir
    character(len=2048) :: candidate
    integer :: exe_len, arg_status, slash_pos
    integer :: env_len, env_status, env_last
    integer :: start_pos, colon_pos
    integer :: access

    if (mace_paths_ready) return

    exe_path = ''
    exe_dir = '.'
    call get_command_argument(0, exe_path, length=exe_len, status=arg_status)
    if (arg_status == 0 .and. exe_len > 0) then
      slash_pos = scan(trim(exe_path), '/', back=.true.)
      if (slash_pos > 0) then
        if (slash_pos == 1) then
          exe_dir = '/'
        else
          exe_dir = exe_path(:slash_pos-1)
        end if
      else
        ! argv[0] can be just "pimd.exe" when launched through PATH.
        call get_environment_variable('PATH', path_env, &
             length=env_len, status=env_status)
        if ((env_status == 0 .or. env_status == -1) .and. env_len > 0) then
          env_last = min(env_len, len(path_env))
          start_pos = 1
          do
            colon_pos = index(path_env(start_pos:env_last), ':')
            if (colon_pos == 0) then
              path_dir = path_env(start_pos:env_last)
            else if (colon_pos == 1) then
              path_dir = '.'
            else
              path_dir = path_env(start_pos:start_pos+colon_pos-2)
            end if

            candidate = trim(path_dir)//'/'//trim(exe_path)
            if (access(trim(candidate), ' ') == 0) then
              exe_dir = trim(path_dir)
              exit
            end if

            if (colon_pos == 0) exit
            start_pos = start_pos + colon_pos
            if (start_pos > env_last) exit
          end do
        end if
      end if
    end if

    mace_python_dir = trim(exe_dir)//'/MACE'
    mace_helper_path = trim(mace_python_dir)//'/mace_function.py'
    mace_model_path = trim(mace_python_dir)//'/small_mace_model.model'
    mace_paths_ready = .true.
  end subroutine resolve_mace_paths

  subroutine get_mace_model_path(model_path)
    character(len=*), intent(out) :: model_path
    integer :: env_len, env_status

    call resolve_mace_paths()

    call get_environment_variable('MACE_MODEL', mace_model_path, &
         length=env_len, status=env_status)
    if (env_status /= 0) mace_model_path = trim(mace_python_dir)//'/small_mace_model.model'

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
    call program_abort('ERROR!!! There is no MACE directory: '//mace_python_dir)
  end if
  if (access(mace_helper_path, ' ') /= 0) then
    call program_abort('ERROR!!! There is no MACE Python helper: '//mace_helper_path)
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
          AU2Ang, eV2AU, eVAng2AU, dp_inv, Lperiodic, lattice
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
