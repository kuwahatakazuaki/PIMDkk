program Path_Integral
  use Parameters
#ifdef _mpi_
  use MPI
#endif
  use utility, only: program_abort
  implicit none
  integer :: ierr

#ifdef _mpi_
  call MPI_INIT(IERR)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,Nproc,IERR)
  call MPI_COMM_RANK(MPI_COMM_WORLD,MyRank,IERR)
#else
  MyRank = 0
  NProc  = 1
#endif

  if ( MyRank == 0 ) then
    call read_parameter
    call print_start
  end if
  call Broad1
  call Set_Allocate

  if ( MyRank == 0 ) then
    call read_structure
  end if
  call Broad2
  call Calc_Constant
  call Check_Inp

! Isimulation = 0:PIMD, 1:RPMD, 2CMD, 10:CL
  select case(Isimulation)
    case(0:2)
      call PI_NEW_MPI
    case(10)
      call Classical
    case default
      print *, '"Simulation" is ', Isimulation
      call program_abort('ERROR!!! Wrong "Simulation" option')
  end select

  call Set_Deallocate
  call print_end
#ifdef _mpi_
  call MPI_FINALIZE(IERR)
#endif

stop

contains

  subroutine print_start
    use utility, only: get_time
    integer :: Uout
    if ( Lrestart .eqv. .False. ) then
      open(newunit=Uout,file=Fout,status='replace')
    else
      open(newunit=Uout,file=Fout,status='old',position='append')
    end if
    !open(newunit=Uout,file=Fout,status='replace')
      write(Uout,*) '***********************'
      write(Uout,*) '   Simulation Start!   '
      write(Uout,*) '***********************'
      write(Uout,*) ' Simulation Started at ', get_time()
      write(Uout,*)
    close(Uout)
  end subroutine print_start

  subroutine print_end
    use utility, only: get_time
    integer :: Uout
    if ( MyRank == 0 ) then
      open(newunit=Uout,file=Fout,status='old',position='append')
        write(Uout,*) '***********************'
        write(Uout,*) '    Simulation End!    '
        write(Uout,*) '***********************'
        write(Uout,*) '   Simulation Ended at ', get_time()
        write(Uout,*)
      close(Uout)
    end if
  end subroutine print_end

end program Path_Integral
