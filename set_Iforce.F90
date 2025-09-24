subroutine set_Iforce
  use Parameters
  use utility, only: program_abort
  implicit none

  select case(Iforce)
    case(1)
      call Set_mopac
    case(6)
      call Set_Gaussian_MPI_tk  ! Set for chk, rwf, etc.
    case(8)
      call Set_VASP
    case(9)
      call Set_siesta
    !case(10)
    !  call Set_Molcas_MPI
    case(21)
      call set_nnp_araidai
    case(22)
      call set_nnp_matlantis
#ifdef _LAMMPS_
    case(24)
      call set_LAMMPS
#endif
  end select

  return

contains

#ifdef _LAMMPS_
  subroutine set_LAMMPS
    implicit none
    character(len=*), parameter :: infile='lammps.data'
    character(len=1024) :: line
    integer :: ios, iu, i
    logical :: found_x, found_y, found_z, found_tilt
    real(8) :: xlo, xhi, ylo, yhi, zlo, zhi
    real(8) :: lx, ly, lz, xy, xz, yz
    real(8) :: lattice(3,3)

    found_x = .false.; found_y = .false.; found_z = .false.; found_tilt = .false.
    xlo = 0d0; xhi = 0d0; ylo = 0d0; yhi = 0d0; zlo = 0d0; zhi = 0d0
    xy  = 0d0; xz  = 0d0; yz  = 0d0

    ! --- open and scan file ---
    open(newunit=iu, file=trim(infile), status='old', action='read', iostat=ios)
    if (ios /= 0) stop 'Error: cannot open file: "lammps.dat"'

    do
      read(iu,'(A)', iostat=ios) line
      if (ios /= 0) exit

      ! find xlo xhi / ylo yhi / zlo zhi lines
      if (.not. found_x) then
        if (index(line, 'xlo') > 0 .and. index(line, 'xhi') > 0) then
          read(line,*, iostat=ios) xlo, xhi
          if (ios == 0) found_x = .true.
          cycle
        end if
      end if

      if (.not. found_y) then
        if (index(line, 'ylo') > 0 .and. index(line, 'yhi') > 0) then
          read(line,*, iostat=ios) ylo, yhi
          if (ios == 0) found_y = .true.
          cycle
        end if
      end if

      if (.not. found_z) then
        if (index(line, 'zlo') > 0 .and. index(line, 'zhi') > 0) then
           read(line,*, iostat=ios) zlo, zhi
           if (ios == 0) found_z = .true.
           cycle
        end if
      end if

      ! find tilt factors (xy xz yz)
      if (.not. found_tilt) then
        if (index(line, 'xy') > 0 .and. index(line, 'xz') > 0 .and. index(line, 'yz') > 0) then
          read(line,*, iostat=ios) xy, xz, yz
          if (ios == 0) found_tilt = .true.
          cycle
        end if
      end if
    end do
    close(iu)

    if (.not.(found_x .and. found_y .and. found_z)) then
       stop 'Error: failed to find all x/y/z bounds in file.'
    end if

    ! --- compute box lengths ---
    lx = xhi - xlo
    ly = yhi - ylo
    lz = zhi - zlo

    !! --- construct VASP cell vectors (POSCAR) ---
    lattice(1,:) = [ lx, 0d0, 0d0 ]
    lattice(2,:) = [ xy, ly,  0d0 ]
    lattice(3,:) = [ xz, yz,  lz ]

  end subroutine set_LAMMPS
#endif

  subroutine Set_siesta
    !use Parameters
    implicit none
    integer   :: i,j,k,imode

    do imode=ista,iend
      write(addresstmp(laddress+1:laddress+6),'(i5.5,A1)') imode,'/'
      call system('mkdir -p '//trim(addresstmp))
      call system('cp InputFile/* '//trim(addresstmp))
    enddo

  return
  end subroutine Set_siesta

end subroutine set_Iforce
