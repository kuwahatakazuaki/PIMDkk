  Subroutine Init_Recv_Send_MPI_tk

    Use Parameters
    Use Parameter_tk
    Use MPI
    Implicit None

    Integer :: i,j,k,ii,jj,isend,irecv,n

    if(myrank/=0) then
       ii=nrecv+1
       irecv=sendlist
       j=sendlist2
       n=j*natom
!       print*,myrank," recv coordinates from ",irecv,j
       Call MPI_Recv_init(x(1,ista),n,MPI_DOUBLE_PRECISION,irecv,900000+myrank,MPI_COMM_WORLD,ireqb(1,ii),IERR)
       Call MPI_Recv_init(y(1,ista),n,MPI_DOUBLE_PRECISION,irecv,1000000+myrank,MPI_COMM_WORLD,ireqb(2,ii),IERR)
       Call MPI_Recv_init(z(1,ista),n,MPI_DOUBLE_PRECISION,irecv,1100000+myrank,MPI_COMM_WORLD,ireqb(3,ii),IERR)
    endif
    Do ii=nrecv,1,-1
       isend=recvlist(ii)
       i=recvlist1(ii)
       j=recvlist2(ii)
       n=j*natom
!       print*,myrank," send coordinates to ",isend,i,j
       Call MPI_Send_init(x(1,i),n,MPI_DOUBLE_PRECISION,isend,900000+isend,MPI_COMM_WORLD,ireqb(1,ii),IERR)
       Call MPI_Send_init(y(1,i),n,MPI_DOUBLE_PRECISION,isend,1000000+isend,MPI_COMM_WORLD,ireqb(2,ii),IERR)
       Call MPI_Send_init(z(1,i),n,MPI_DOUBLE_PRECISION,isend,1100000+isend,MPI_COMM_WORLD,ireqb(3,ii),IERR)
    Enddo

    Return
  End Subroutine
