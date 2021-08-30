program pimd
  use global_variable !,  only: Nproc, myrank, physmass
  use communication
  use mpi
  implicit none
  integer :: ierr
  integer :: i, j

  call mpi_init(ierr)
  call mpi_comm_size(MPI_COMM_WORLD, Nproc, ierr)
  call mpi_comm_rank(MPI_COMM_WORLD, myrank, ierr)

  if (myrank == 0) then
    open(newunit=Uout,file=Oname,status='replace')
      write(Uout,'(a)') ' ***********************'
      write(Uout,'(a)') '    Simulation Start!   '
      write(Uout,'(a)') ' ***********************'
      write(Uout,'(a)') '  Simulation Started at '
      call print_time()
      call read_input1
    close(Uout)
  end if
  call broad_parameters1
  call set_allocate

  if (myrank == 0) call read_input2
  call broad_parameters2
  call setup_whole
  call setup_mpi
  call setup_indivi

  if (myrank == 0) call check_input

  call normal_mode_matrix
  call init_mass

  if (Lrestart .eqv. .True.) then
!    if (myrank == 0) call restart_read
!    call broad_parameters3
  else
    if ( myrank == 0) call print_ini
    call init_position
    if ( myrank == 0 ) then
      call init_velocity
      call init_bath
    end if
    call broad_parameters3
    call normal_mode_trans ! x(i) = x(i) + sum_j tnm(i,j)*u(j)

    Istep = 0
    call get_force
    call normal_mode_force ! fu(i) = fu(i) + sum_j fx(j)*tnm(j,i)
    if ( myrank == 0 ) then
      call calc_hamil
    end if
  end if
  if ( myrank == 0 ) then
  call getforce_ref
  end if

  do Istep = Nrestart, Nstep
    if ( myrank == 0 ) then
      call integ_nhc_cent
    else
    end if
  end do

  if (myrank == 0) then
    open(newunit=Uout,file=Oname,status='old',position='append')
      write(Uout,'(a)') ' ***********************'
      write(Uout,'(a)') '     Simulation End!    '
      write(Uout,'(a)') ' ***********************'
      write(Uout,'(a)') '   Simulation Ended at  '
      call print_time()
    close(Uout)
  end if

  call mpi_barrier(MPI_COMM_WORLD,ierr)
  call mpi_finalize(ierr)

contains

!  do i = 0, Nproc-1
!    if (myrank == i) then
!      print *, "MyRank is ", myrank
!      print *, physmass(:)
!      print *, alabel(:)
!    end if
!    call mpi_barrier(MPI_COMM_WORLD, ierr)
!  end do

  subroutine print_time()
    integer :: newtime(8)
    character(len=10) :: date
    character(len=10) :: time
    character(len=10) :: zone

    call date_and_time(date, time, zone, newtime)
    write(Uout,'("  ",i4,a1,i2.2,a1,i2.2,a1,i2.2,a1,i2.2,a1,i2.2,a1,i3.3)') &
             newtime(1),'/',newtime(2),'/',newtime(3),' ', &
             newtime(5),':',newtime(6),':',newtime(7),':',newtime(8)
    write(Uout,'()')
  end subroutine print_time


end program pimd

