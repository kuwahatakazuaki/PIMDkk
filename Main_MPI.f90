Program Path_Integral_MPI

  Use Parameters
  Use MPI
  Implicit None
  Integer,dimension(8) :: newtime
  Character (Len=10) :: date
  Character (Len=10) :: time
  Character (Len=10) :: zone

  Call MPI_INIT(IERR)
  Call MPI_COMM_SIZE(MPI_COMM_WORLD,NProcs,IERR)
  Call MPI_COMM_RANK(MPI_COMM_WORLD,MyRank,IERR)

  If(MyRank==0) Then
     Write(*,*) '***********************'
     Write(*,*) '   Simulation Start!   '
     Write(*,*) '***********************'
     Write(*,*) ' Simulation Started at '
     Call Date_and_Time(date,time,zone,newtime)
     write(*,'(i4,a1,i2.2,a1,i2.2,a1,i2.2,a1,i2.2,a1,i2.2,a1,i3.3)') &
              newtime(1),'/',newtime(2),'/',newtime(3),' ', &
              newtime(5),':',newtime(6),':',newtime(7),':',newtime(8)
     Call Read_Inp1
  EndIf
  Call Broad1
  Call Set_Allocate

  If(MyRank==0) Then
     Call Read_Inp2
  EndIf
  Call Broad2

  select case(Simulation)
    case(0)
      call PI_NEW_MPI
!    case(1)
!      call PI_Harmonic_Oscillator
!    case(2)
!      call PI_Spline
    case(3)
      call CMD_NEW_MPI
    case(4)
      call RPMD_NEW_MPI
    case(10)
      call Classical
!    case(99)
!      call SHOOT_NEW_MPI
    case default
      print *, '"Simulation" is ', Simulation
      stop 'ERROR!!! Wrong "Simulation" option'
  end select

  If(Simulation==10) Then
    Call Set_Deallocate_Classical
  Else
    Call Set_Deallocate
  EndIf
!   Call Set_Deallocate

  If(MyRank==0) Then
     Write(*,*) '***********************'
     Write(*,*) '    Simulation End!    '
     Write(*,*) '***********************'
     Write(*,*) '  Simulation Ended at  '

     Call Date_and_Time(date,time,zone,newtime)
     write(*,'(i4,a1,i2.2,a1,i2.2,a1,i2.2,a1,i2.2,a1,i2.2,a1,i3.3)') &
           newtime(1),'/',newtime(2),'/',newtime(3),' ', &
           newtime(5),':',newtime(6),':',newtime(7),':',newtime(8)
  EndIf
  Call MPI_FINALIZE(IERR)

stop
End Program
