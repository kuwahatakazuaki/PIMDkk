module utility
  implicit none

contains
  subroutine program_abort(message)
    use mpi
    character(*) :: message
    integer :: ierr

    print *, message
    call mpi_barrier(MPI_COMM_WORLD, ierr)
    call mpi_abort(MPI_COMM_WORLD, -1, ierr)
    stop
  end subroutine program_abort

  subroutine gasdev(gasd)
    real(8), intent(inout) :: gasd
    integer, save :: iset = 0
    real(8) :: v1, v2, rsq, fac
    real(8), save :: gset
    real(8) :: ran
    integer :: i

    if ( iset == 0 ) then
      do
        call random_generator(1,ran)
        v1 = 2.0 * ran - 1.0
        call random_generator(1,ran)
        v2 = 2.0 * ran - 1.0
        rsq = v1 * v1 + v2 * v2
!        if ( rsq > 0.0d0 .and. rsq < 1.0d0 ) then
        if ( rsq < 1.0d0 ) then
          fac = dsqrt( -2.0 * dlog(rsq) / rsq )
          gset = v1 * fac
          gasd = v2 * fac
          exit
        end if
      end do
      iset = 1
    else
      gasd = gset
      iset = 0
    end if
  end subroutine gasdev


  subroutine set_input_file(Ifile,Def_file)
    integer :: leng
    character(:), allocatable, intent(out):: Ifile
    character(*), optional :: Def_file
    if ( command_argument_count() == 0) then
      Ifile = Def_file
!      print *, "Reading from default file as ", Ifile
!      print *, 'ERROR!! There is no input file'
    else
      call get_command_argument(1, length=leng)
        allocate(character(leng) :: Ifile)
        call get_command_argument(1, Ifile)
!      print *, "Reading from ", Ifile
    end if
    return
  end subroutine set_input_file

  real(8) function kinetic_energy()
    use global_variable
    real(8) :: kine
    integer :: i, j
    kine = 0.0d0
    do i = 1, Natom
      do j = 1, Nbead
        kine = kine + fictmass(i,j) * dot_product(vu(:,i,j),vu(:,i,j))
      end do
    end do
  ! --- Be care!! Original program doesn't have the factor of 0.5 ---
    kinetic_energy = 0.5*kine
  end function kinetic_energy

!  subroutine random_seed_clock()
!    integer :: nseed, clock
!    integer, allocatable :: seed(:)
!
!    call system_clock(clock)
!
!    call random_seed(size=nseed)
!    allocate(seed(nseed))
!    seed = clock
!    call random_seed(put=seed)
!  end subroutine random_seed_clock

end module utility


