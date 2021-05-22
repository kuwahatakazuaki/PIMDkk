  Subroutine Recv_Send_MPI_tk4

    Use Parameters
    Use Parameter_tk
    Use MPI
    Implicit None

    Integer :: i,j,k,ii,jj,isend,irecv,n

    if(myrank/=0) then
       irecv=sendlist
       j=sendlist2
       n=j*natom
!       print*,myrank," recv coordinates from ",irecv,j
       Call MPI_Irecv(x(1,ista),n,MPI_DOUBLE_PRECISION,irecv,900000+myrank,MPI_COMM_WORLD,ireq(1),IERR)
       Call MPI_Irecv(y(1,ista),n,MPI_DOUBLE_PRECISION,irecv,1000000+myrank,MPI_COMM_WORLD,ireq(2),IERR)
       Call MPI_Irecv(z(1,ista),n,MPI_DOUBLE_PRECISION,irecv,1100000+myrank,MPI_COMM_WORLD,ireq(3),IERR)
       Call MPI_Wait(ireq(1),mstatus,IERR)
       Call MPI_Wait(ireq(2),mstatus,IERR)
       Call MPI_Wait(ireq(3),mstatus,IERR)
    endif
    Do ii=nrecv,1,-1
       isend=recvlist(ii)
       i=recvlist1(ii)
       j=recvlist2(ii)
       n=j*natom
!       print*,myrank," send coordinates to ",isend,i,j
       Call MPI_Isend(x(1,i),n,MPI_DOUBLE_PRECISION,isend,900000+isend,MPI_COMM_WORLD,ireq(1),IERR)
       Call MPI_Isend(y(1,i),n,MPI_DOUBLE_PRECISION,isend,1000000+isend,MPI_COMM_WORLD,ireq(2),IERR)
       Call MPI_Isend(z(1,i),n,MPI_DOUBLE_PRECISION,isend,1100000+isend,MPI_COMM_WORLD,ireq(3),IERR)
       Call MPI_Wait(ireq(1),mstatus,IERR)
       Call MPI_Wait(ireq(2),mstatus,IERR)
       Call MPI_Wait(ireq(3),mstatus,IERR)
    Enddo

    Return
  End Subroutine
