subroutine read_input1
  use global_variable
  use utility
  use mpi
  implicit none
  integer :: Uinp, ios, ierr
  integer :: i
  character(len=90) :: line
  character(len=:), allocatable :: Ifile


  call set_input_file(Ifile, 'input')

  open(newunit=Uinp, file=Ifile, iostat=ios,status='old')
    if (ios /= 0) then
      print *, 'Failed to open ', Ifile
      call mpi_abort(MPI_COMM_WORLD,1,ierr)
      stop "Filed to open the input file"
    end if

! +++ Reading required parameters +++
    do
      read(Uinp,'(a)',iostat=ios) line
      if (ios < 0) exit

      if     (index(line,'$natom') > 0 ) then
        read(Uinp,*) Natom
      elseif (index(line,'$Temperature') > 0) then
        read(Uinp,*) temperature
      elseif (index(line,'$nbead') > 0) then
        read(Uinp,*) Nbead
      elseif (index(line,'$nstep') > 0) then
        read(Uinp,*) Nstep
      elseif (index(line,'$dt') > 0) then
        read(Uinp,*) dt
      elseif (index(line,'$nforce') > 0) then
        read(Uinp,*) Nforce
      elseif (index(line,'$Simulation') > 0) then
        read(Uinp,*) simulation
      elseif (index(line,'$nrestart') > 0) then
        read(Uinp,*) Lrestart
      elseif (index(line,'$path_result') > 0) then
!        read(Uinp,'(a)') path_result_temp
        read(Uinp,'(a)') line
        path_result = trim(line)
!        print *, path_result
!        call program_abort("HERE")
      elseif (index(line,'$seed') > 0) then
        do i = 1, 4
          read(Uinp,*) iseed(i)
        end do
      end if
    end do
! +++ Reading required parameters +++

! +++ Check the input parameters
!    if ( .not. allocated(path_result)) then
!    if ( trim(path_result_temp) == "0" ) then
!      print *, 'ERROR!! path_result is not exits'
!      call mpi_abort(MPI_COMM_WORLD,1,ierr)
    if ( Natom == 0) then
      print *, 'ERROR!! "natom" is not defiled!'
      call mpi_abort(MPI_COMM_WORLD,1,ierr)
    else if ( Nbead == 0) then
      print *, 'ERROR!! "nbead" is not defiled!'
      call mpi_abort(MPI_COMM_WORLD,1,ierr)
    else if ( nforce == 0) then
      print *, 'ERROR!! "nforce" is not defiled!'
      call mpi_abort(MPI_COMM_WORLD,1,ierr)
    else if ( simulation == -1) then
      print *, 'ERROR!! "simulation" is not defiled!'
      call mpi_abort(MPI_COMM_WORLD,1,ierr)
    end if
! +++ End Check the input parameters


  rewind(Uinp)
! +++ Reading optional parameters +++
    do
      read(Uinp,*,iostat=ios) line
      if (ios < 0) exit

      if     (index(line,'$nref') > 0 ) then
        read(Uinp,*) Nref
      elseif (index(line,'$Nys') > 0) then
        read(Uinp,*) Nys
      elseif (index(line,'$nnhc') > 0) then
        read(Uinp,*) Nnhc
      elseif (index(line,'$Order') > 0) then
        read(Uinp,*) order
      elseif (index(line,'$gamma') > 0) then
        read(Uinp,*) gamma
      elseif (index(line,'$nensemble') > 0) then
        read(Uinp,*) Nensemble
      elseif (index(line,'$theory') > 0) then
        read(Uinp,*) Ntheory
      elseif (index(line,'$ncent') > 0) then
        read(Uinp,*) Ncent
      elseif (index(line,'$ncolor') > 0) then
        read(Uinp,*) Ncolor
      elseif (index(line,'$angstrom') > 0) then
        read(Uinp,*) Langstrom
      elseif (index(line,'$randomc') > 0) then
        read(Uinp,*) Lrandomc
      elseif (index(line,'$ngengau') > 0) then
        read(Uinp,*) Lgengau
      elseif (index(line,'$path_scr') > 0) then
!        read(Uinp,'(a)') path_scr_temp
        read(Uinp,'(a)') line
        path_scr = trim(line)
      elseif (index(line,'$write.charge') > 0) then
        read(Uinp,*) Lcharge
      elseif (index(line,'$write.hfcc') > 0) then
        read(Uinp,*) Lhfcc
      elseif (index(line,'$write.homolumo') > 0) then
        read(Uinp,*) Lhomolumo
      elseif (index(line,'$write.dipole') > 0) then
        read(Uinp,*) Ldipole
      elseif (index(line,'$write.pop') > 0) then
        read(Uinp,*) Lpop
      elseif (index(line,'$write.force') > 0) then
        read(Uinp,*) Lforce
      end if
    end do
  close(Uinp)
! +++ End Reading optional parameters +++

  dt = dt * facttime

end subroutine read_input1


subroutine read_input2
  use global_variable
  use utility
  use mpi
  implicit none
  integer :: Uinp=20, ios, ierr
  integer :: i
  character(len=90) :: line
  character(len=:), allocatable :: Ifile
  call set_input_file(Ifile, 'input')

! +++ Reading coordinate +++
  open(newunit=Uinp, file=Ifile, iostat=ios,status='old')
    do
!      read(Uinp,'(a)',iostat=ios) line
      read(Uinp,'(a)',end=100) line
!      if (ios < 0) exit

      if     (index(line,'$Coords') > 0 ) then
        do i = 1, Natom
          read(Uinp,*,err=400) alabel(i), physmass(i), u(1,i,1), u(2,i,1), u(3,i,1)
        end do
        exit
      end if
    end do
  close(Uinp)
! +++ End Reading coordinate +++
!print *, "In Reading"
!do i = 1, Natom
!  print *, "alabel ", alabel(i)
!end do

  physmass(:) = physmass(:) * factmass
  if ( Langstrom .eqv. .True.) u(:,:,1) = u(:,:,1) * AngtoAU

return
100 continue
  print *, 'ERROR!! "$Coords" is not exist'
  call mpi_abort(MPI_COMM_WORLD,1,ierr)
400 continue
  print *, 'ERROR!! during the reading of the coordinate'
  call mpi_abort(MPI_COMM_WORLD,1,ierr)
end subroutine read_input2

