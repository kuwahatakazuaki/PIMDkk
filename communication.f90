module communication
  use global_variable
  use mpi
  use utility
  implicit none
  integer :: Nrecv
  integer, allocatable :: ireq(:,:)

contains

subroutine setup_mpi
  integer :: i, j, n, n0
  integer :: ierr
  integer :: listeachtmp(Nproc)
  integer :: listeach(Nproc)
  integer :: Nhmod, Neach, Nsendrecv
  integer, allocatable :: recvlist(:), recvIndex(:), recvSize(:)
  integer :: sendlist, sendIndex, sendSize

  Nhmod = mod(Nbead, Nproc)
  Neach = (Nbead-Nhmod) / Nproc

  do i = 1, Nhmod
    listeach(i) = Neach + 1
  end do
  do i = Nhmod+1, Nproc
    listeach(i) = Neach
  end do
  Neach = listeach(myrank+1)

  Ista = 1
  do i = 1, myrank
    Ista = Ista + listeach(i)
  end do
  Iend = Ista + Neach - 1
! print *, "myrank, Isend, Iend",myrank, Ista, Iend

  i = 1
  j = Nproc
  do
    j = (j+1) / 2
    if ( j == 1) then
      Nsendrecv = i
      exit
    else
      i = i + 1
    end if
  end do

  allocate( recvlist(Nsendrecv) )
  allocate( recvIndex(Nsendrecv))
  allocate( recvSize(Nsendrecv))

  Nrecv = 0
  n = 1
  sendlist  = 0
  sendIndex = Ista
  sendSize = 0
  do i = myrank+1, Nproc
    listeachtmp(i) = listeach(i)
  end do

  do i = 1, Nsendrecv
    n0 = n
    n  = n0 * 2
    if ( (mod(myrank,n) == 0) .and. (myrank+n0 < Nproc) ) then
      Nrecv = Nrecv + 1
      recvlist(Nrecv)  = myrank + n0
      recvIndex(Nrecv) = Ista + listeachtmp(myrank+1)
      recvSize(Nrecv) = listeachtmp(myrank+n0+1)
! print *, myrank, " get from", recvlist(Nrecv), recvIndex(Nrecv), recvSize(Nrecv)
    else if ( mod(myrank,n) == n0 ) then
      sendlist  = myrank - n0
      sendIndex = Ista
      sendSize = listeachtmp(myrank+1)
! print *, myrank, " send to ", sendlist, sendIndex, sendSize
    end if
    do j = myrank+1, Nproc-n0,n
      listeachtmp(j) = listeachtmp(j) + listeachtmp(j+n0)
    end do
  end do

! print *, "myrank Nrecv", myrank, Nrecv
  allocate(ireq(6,Nrecv+1))



!!! Debug
!!r(:,:,:) = 0.0d0
!!r(:,:,myrank+1) = myrank+1
!!do i = 1, Nproc
!!  if (myrank+1 == i) then
!!    print *, myrank, ":",r(:,:,i)
!!  end if
!!  call mpi_barrier(MPI_COMM_WORLD, ierr)
!!end do
!!  call mpi_barrier(MPI_COMM_WORLD, ierr)
!!
!!if ( myrank == 0 ) then
!!  print *, "myrank in 0"
!!  do i = 1, Nbead
!!    print *, r(:,:,i)
!!  end do
!!end if
!!  call mpi_barrier(MPI_COMM_WORLD, ierr)
!!! Debug

block
  integer :: Irecv, Isend
  integer :: sender, receiver
  integer :: Ssend, Srecv
  integer :: ierr
  do i = 1, Nrecv
    sender = recvlist(i)
    Isend  = recvIndex(i)
    Srecv  = recvSize(i)
    call mpi_recv_init(f(1,1,Isend), Srecv*3*Natom, MPI_DOUBLE_PRECISION, sender,sender+1000, MPI_COMM_WORLD, ireq(1,i), ierr)
    call mpi_recv_init(energy(Isend),Srecv,         MPI_DOUBLE_PRECISION, sender,sender+2000, MPI_COMM_WORLD, ireq(2,i), ierr)
    call mpi_recv_init(r(1,1,Isend), Srecv*3*Natom, MPI_DOUBLE_PRECISION, sender,sender+3000, MPI_COMM_WORLD, ireq(3,i), ierr)
    if ( Ldipole .eqv. .True. ) then
      call mpi_recv_init(dipole(1,Isend), Srecv*4, MPI_DOUBLE_PRECISION, sender,sender+4000, MPI_COMM_WORLD, ireq(4,i), ierr)
    end if
    if ( Lcharge .eqv. .True. ) then
      call mpi_recv_init(charge(1,Isend), Srecv*Natom, MPI_DOUBLE_PRECISION, sender,sender+5000, MPI_COMM_WORLD, ireq(5,i), ierr)
    end if
  end do
  if ( myrank /= 0 ) then
    receiver = sendlist
    call mpi_send_init(f(1,1,Ista),sendSize*3*Natom,MPI_DOUBLE_PRECISION,receiver,myrank+1000,MPI_COMM_WORLD,ireq(1,Nrecv+1),ierr)
    call mpi_send_init(energy(Ista),sendSize,       MPI_DOUBLE_PRECISION,receiver,myrank+2000,MPI_COMM_WORLD,ireq(2,Nrecv+1),ierr)
    call mpi_send_init(r(1,1,Ista),sendSize*3*Natom,MPI_DOUBLE_PRECISION,receiver,myrank+3000,MPI_COMM_WORLD,ireq(3,Nrecv+1),ierr)
    if ( Ldipole .eqv. .True. ) then
      call mpi_send_init(dipole(1,Ista),sendSize*4,MPI_DOUBLE_PRECISION,receiver,myrank+4000,MPI_COMM_WORLD,ireq(4,Nrecv+1),ierr)
    end if
    if ( Lcharge .eqv. .True. ) then
      call mpi_send_init(charge(1,Ista),sendSize*Natom,MPI_DOUBLE_PRECISION,receiver,myrank+5000,MPI_COMM_WORLD,ireq(5,Nrecv+1),ierr)
    end if
  end if

!call mpi_barrier(MPI_COMM_WORLD, ierr)
!do i = 1, Nrecv
!  call mpi_start(ireq(1,i), ierr)
!  call mpi_start(ireq(2,i), ierr)
!  call mpi_start(ireq(3,i), ierr)
!end do
!do i = 1, Nrecv
!  call mpi_wait(ireq(1,i), MPI_STATUS_IGNORE, ierr)
!  call mpi_wait(ireq(2,i), MPI_STATUS_IGNORE, ierr)
!  call mpi_wait(ireq(3,i), MPI_STATUS_IGNORE, ierr)
!end do
!if ( myrank /= 0 ) then
!  call mpi_start(ireq(1,Nrecv+1), ierr)
!  call mpi_start(ireq(2,Nrecv+1), ierr)
!  call mpi_start(ireq(3,Nrecv+1), ierr)
!  call mpi_wait( ireq(1,Nrecv+1), MPI_STATUS_IGNORE, ierr)
!  call mpi_wait( ireq(2,Nrecv+1), MPI_STATUS_IGNORE, ierr)
!  call mpi_wait( ireq(3,Nrecv+1), MPI_STATUS_IGNORE, ierr)
!end if

!if ( myrank == 0 ) then
!  print *, "after myrank in 0"
!  do i = 1, Nbead
!    print *, r(:,:,i)
!  end do
!end if
end block

end subroutine setup_mpi

subroutine comm_coord
  integer :: i
  integer :: ierr
  do i = 1, Nrecv
    call mpi_start(ireq(3,i), ierr)
  end do
  do i = 1, Nrecv
    call mpi_wait(ireq(3,i), MPI_STATUS_IGNORE, ierr)
  end do
  if ( myrank /= 0 ) then
    call mpi_start(ireq(3,Nrecv+1), ierr)
    call mpi_wait( ireq(3,Nrecv+1), MPI_STATUS_IGNORE, ierr)
  end if
end subroutine comm_coord

subroutine comm_output
  integer :: i, ierr
  do i = 1, Nrecv
    call mpi_start(ireq(1,i), ierr)
    call mpi_start(ireq(2,i), ierr)
    if ( Ldipole .eqv. .True. ) then
      call mpi_start(ireq(4,i), ierr)
    end if
    if ( Lcharge .eqv. .True. ) then
      call mpi_start(ireq(5,i), ierr)
    end if
  end do
  do i = 1, Nrecv
    call mpi_wait(ireq(1,i), MPI_STATUS_IGNORE, ierr)
    call mpi_wait(ireq(2,i), MPI_STATUS_IGNORE, ierr)
    if ( Ldipole .eqv. .True. ) then
      call mpi_wait(ireq(4,i), MPI_STATUS_IGNORE, ierr)
    end if
    if ( Lcharge .eqv. .True. ) then
      call mpi_wait(ireq(5,i), MPI_STATUS_IGNORE, ierr)
    end if
  end do
  if ( myrank /= 0 ) then
    call mpi_start(ireq(1,Nrecv+1), ierr)
    call mpi_start(ireq(2,Nrecv+1), ierr)
    call mpi_wait( ireq(1,Nrecv+1), MPI_STATUS_IGNORE, ierr)
    call mpi_wait( ireq(2,Nrecv+1), MPI_STATUS_IGNORE, ierr)
    if ( Ldipole .eqv. .True. ) then
      call mpi_start(ireq(4,Nrecv+1), ierr)
      call mpi_wait( ireq(4,Nrecv+1), MPI_STATUS_IGNORE, ierr)
    end if
    if ( Lcharge .eqv. .True. ) then
      call mpi_start(ireq(5,Nrecv+1), ierr)
      call mpi_wait( ireq(5,Nrecv+1), MPI_STATUS_IGNORE, ierr)
    end if
  end if
end subroutine comm_output

end module


