Program Path_Integral_MPI
  use Parameters
  use MPI
  use utility, only: program_abort
  implicit none

  call MPI_INIT(IERR)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,NProcs,IERR)
  call MPI_COMM_RANK(MPI_COMM_WORLD,MyRank,IERR)

  call print_start
  if ( MyRank == 0 ) then
    call read_parameter
  end if
  call Broad1
  call Set_Allocate

  if ( MyRank == 0 ) then
    call read_structure
  end if
  call Broad2
  call Calc_Constant
  call Check_Inp

  select case(Isimulation)
    case(0)
      call PI_NEW_MPI
    case(1)
    !  call RPMD_NEW_MPI
    !case(2)
    !  call CMD_NEW_MPI
    case(10)
      call Classical
    case default
      print *, '"Simulation" is ', Isimulation
      call program_abort('ERROR!!! Wrong "Simulation" option')
  end select

  call Set_Deallocate
  call print_end
  call MPI_FINALIZE(IERR)

stop

contains

  subroutine print_start
    use utility, only: get_time
    integer :: Uout
    if ( MyRank == 0 ) then
      !if ( Lrestart .eqv. .False. ) then
      !  open(newunit=Uout,file=Fout,status='replace')
      !else
      !  open(newunit=Uout,file=Fout,status='old',position='append')
      !end if
      open(newunit=Uout,file=Fout,status='replace')
        write(Uout,*) '***********************'
        write(Uout,*) '   Simulation Start!   '
        write(Uout,*) '***********************'
        write(Uout,*) ' Simulation Started at ', get_time()
        write(Uout,*)
      close(Uout)
    end if
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

End Program
