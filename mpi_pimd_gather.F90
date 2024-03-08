subroutine mpi_pimd_gather
  use Parameters
#ifdef _mpi_
  use MPI
  use utility, only: program_abort
  implicit none
  integer :: i, j, k, Irank
  integer :: mstatus(MPI_STATUS_SIZE)
  integer :: ierr
  integer :: MyNbead, Ncomm

  MyNbead = Iend - Ista + 1

!if (MyRank == 0) print *, "Before"
!do Irank = 0, Nprocs
!  if (MyRank == Irank) then
!    print *, "MyRank Ista, Iend, MyNbead", MyRank, Ista, Iend, MyNbead
!    do j = 1, Nbead
!      do i = 1, Natom
!        print *, j, fr(:,i,j)
!      end do
!    end do
!  end if
!  call MPI_Barrier(MPI_COMM_WORLD,Ierr)
!end do

  if ( MyRank == 0 ) then
    call MPI_Gather(MPI_IN_PLACE,3*Natom*MyNbead,MPI_DOUBLE_PRECISION, &
                    fr(1,1,Ista),3*Natom*MyNbead,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  else
    call MPI_Gather(fr(1,1,Ista),3*Natom*MyNbead,MPI_DOUBLE_PRECISION, &
                    fr(1,1,Ista),3*Natom*MyNbead,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  end if

!if (MyRank == 0) then
!print *, "After"
!do j = 1, Nbead
!  do i = 1, Natom
!    print *, j, fr(:,i,j)
!  end do
!end do
!end if
!call program_abort("HERE_MPI")

  if ( Iforce == 8 ) then
    if ( MyRank == 0 ) then
    call MPI_Gather(MPI_IN_PLACE,listeach(MyRank+1),MPI_DOUBLE_PRECISION, & 
                    pressure(ista),listeach(MyRank+1),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    else
    call MPI_Gather(pressure(ista),listeach(MyRank+1),MPI_DOUBLE_PRECISION, & 
                    pressure(ista),listeach(MyRank+1),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    end if
  end if

  !Do ii=1,nrecv
  !   Call MPI_Start(ireqa(1,ii),IERR)
  !   Call MPI_Start(ireqa(2,ii),IERR)
  !   Call MPI_Start(ireqa(3,ii),IERR)
  !Enddo
  !Do i=1,nrecv
  !   Call MPI_Wait(ireqa(1,i),mstatus,IERR)
  !   Call MPI_Wait(ireqa(2,i),mstatus,IERR)
  !   Call MPI_Wait(ireqa(3,i),mstatus,IERR)
  !Enddo
  !if(myrank/=0) then
  !   ii=nrecv+1
  !   Call MPI_Start(ireqa(1,ii),IERR)
  !   Call MPI_Start(ireqa(2,ii),IERR)
  !   Call MPI_Start(ireqa(3,ii),IERR)
  !   Call MPI_Wait(ireqa(1,ii),mstatus,IERR)
  !   Call MPI_Wait(ireqa(2,ii),mstatus,IERR)
  !   Call MPI_Wait(ireqa(3,ii),mstatus,IERR)
  !endif

  !Do ii=1,nrecv
  !   Call MPI_Start(ireqa(4,ii),IERR)
  ! !if(nocharge==0) then
  ! if( Lsave_charge .eqv. .True.) then
  !   Call MPI_Start(ireqa(8,ii),IERR)
  ! endif
  ! !if(nohfcc==0) then
  ! if( Lsave_hfcc .eqv. .True. ) then
  !   Call MPI_Start(ireqa(9,ii),IERR)
  ! endif
  !Enddo

  !Do i=1,nrecv
  !     Call MPI_Wait(ireqa(4,i),mstatus,IERR)
  !   if( Lsave_charge .eqv. .True.) then
  !     Call MPI_Wait(ireqa(8,i),mstatus,IERR)
  !   endif
  !   if( Lsave_hfcc .eqv. .True. ) then
  !     Call MPI_Wait(ireqa(9,i),mstatus,IERR)
  !   endif
  !Enddo

  !if(myrank/=0) then
  !   ii=nrecv+1
  !   Call MPI_Start(ireqa(4,ii),IERR)
  !  if( Lsave_charge .eqv. .True.) then
  !   Call MPI_Start(ireqa(8,ii),IERR)
  !  endif
  !  if( Lsave_hfcc .eqv. .True. ) then
  !   Call MPI_Start(ireqa(9,ii),IERR)
  !  endif
  !   Call MPI_Wait(ireqa(4,ii),mstatus,IERR)
  !  if( Lsave_charge .eqv. .True.) then
  !   Call MPI_Wait(ireqa(8,ii),mstatus,IERR)
  !  endif
  !  if( Lsave_hfcc .eqv. .True. ) then
  !   Call MPI_Wait(ireqa(9,ii),mstatus,IERR)
  !  endif
  !endif


  !if( Lsave_dipole .eqv. .True. ) then
  !  Do ii=1,nrecv
  !     Call MPI_Start(ireqa(5,ii),IERR)
  !     Call MPI_Start(ireqa(6,ii),IERR)
  !     Call MPI_Start(ireqa(7,ii),IERR)
  !     Call MPI_Start(ireqa(10,ii),IERR)
  !  Enddo
  !  Do i=1,nrecv
  !     Call MPI_Wait(ireqa(5,i),mstatus,IERR)
  !     Call MPI_Wait(ireqa(6,i),mstatus,IERR)
  !     Call MPI_Wait(ireqa(7,i),mstatus,IERR)
  !     Call MPI_Wait(ireqa(10,i),mstatus,IERR)
  !  Enddo
  !  if(myrank/=0) then
  !     ii=nrecv+1
  !     Call MPI_Start(ireqa(5,ii),IERR)
  !     Call MPI_Start(ireqa(6,ii),IERR)
  !     Call MPI_Start(ireqa(7,ii),IERR)
  !     Call MPI_Start(ireqa(10,ii),IERR)
  !     Call MPI_Wait(ireqa(5,ii),mstatus,IERR)
  !     Call MPI_Wait(ireqa(6,ii),mstatus,IERR)
  !     Call MPI_Wait(ireqa(7,ii),mstatus,IERR)
  !     Call MPI_Wait(ireqa(10,ii),mstatus,IERR)
  !  endif
  !endif
#endif

return
end subroutine mpi_pimd_gather
