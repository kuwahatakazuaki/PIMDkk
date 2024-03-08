Subroutine Init_Recv_Send_MPI_tk
  Use Parameters
#ifdef _mpi_
  Use MPI
  Implicit None
  Integer :: i,j,k,ii,jj,isend,irecv,n
  integer :: ierr

  if(myrank/=0) then
    ii=nrecv+1
    irecv=sendlist
    j=sendlist2
    n=j*natom
    call MPI_Recv_init(r(1,1,ista),3*n,MPI_DOUBLE_PRECISION,irecv,90000+myrank,MPI_COMM_WORLD,ireqb(1,ii),IERR)
  endif
  Do ii=nrecv,1,-1
    isend=recvlist(ii)
    i=recvlist1(ii)
    j=recvlist2(ii)
    n=j*natom
    Call MPI_Send_init(r(1,1,i),3*n,MPI_DOUBLE_PRECISION,isend,90000+isend,MPI_COMM_WORLD,ireqb(1,ii),IERR)
  Enddo
#endif

  Return
End Subroutine
