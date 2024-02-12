  Subroutine Start_Recv_Send_MPI_tk
    Use Parameters
    Use MPI

    Implicit None
    Integer :: i,j,k,ii,jj,isend,irecv,n
    Integer :: mstatus(MPI_STATUS_SIZE)

    if(myrank/=0) then
       ii=nrecv+1
       Call MPI_Start(ireqb(1,ii),IERR)
       Call MPI_Start(ireqb(2,ii),IERR)
       Call MPI_Start(ireqb(3,ii),IERR)
       Call MPI_Wait(ireqb(1,ii),mstatus,IERR)
       Call MPI_Wait(ireqb(2,ii),mstatus,IERR)
       Call MPI_Wait(ireqb(3,ii),mstatus,IERR)
    endif
    Do ii=nrecv,1,-1
       Call MPI_Start(ireqb(1,ii),IERR)
       Call MPI_Start(ireqb(2,ii),IERR)
       Call MPI_Start(ireqb(3,ii),IERR)
    Enddo
    Do ii=nrecv,1,-1
       Call MPI_Wait(ireqb(1,ii),mstatus,IERR)
       Call MPI_Wait(ireqb(2,ii),mstatus,IERR)
       Call MPI_Wait(ireqb(3,ii),mstatus,IERR)
    Enddo

    Return
  End Subroutine
