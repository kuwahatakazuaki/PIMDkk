Subroutine Init_Send_Recv_MPI_tk

  Use Parameters
  Use Parameter_tk
  Use MPI
  Implicit None

  Integer :: i,j,k,ii,jj,isend,irecv,n

!    print*,myrank," test ",nrecv
!    wait(i)
  Do ii=1,nrecv
     isend=recvlist(ii)
     i=recvlist1(ii)
     j=recvlist2(ii)
     n=j*natom
!       print*,myrank," recv properties from ",isend,i,j
     Call MPI_Recv_init(fx(1,i),   n,MPI_DOUBLE_PRECISION,isend,200000+isend,MPI_COMM_WORLD,ireqa(1,ii),IERR)
     Call MPI_Recv_init(fy(1,i),   n,MPI_DOUBLE_PRECISION,isend,300000+isend,MPI_COMM_WORLD,ireqa(2,ii),IERR)
     Call MPI_Recv_init(fz(1,i),   n,MPI_DOUBLE_PRECISION,isend,400000+isend,MPI_COMM_WORLD,ireqa(3,ii),IERR)
     Call MPI_Recv_init(Eenergy(i),j,MPI_DOUBLE_PRECISION,isend,100000+isend,MPI_COMM_WORLD,ireqa(4,ii),IERR)
     if(nocharge==0) then
       Call MPI_Recv_init(charge(1,i),n,MPI_DOUBLE_PRECISION,isend,800000+isend,MPI_COMM_WORLD,ireqa(8,ii),IERR)
     endif
     if(nohfcc==0) then
       Call MPI_Recv_init(hfcc(1,i),n,MPI_DOUBLE_PRECISION,isend,900000+isend,MPI_COMM_WORLD,ireqa(9,ii),IERR)
     endif

    if(nodipole==0) then
        Call MPI_Recv_init(dipolex(i),j,MPI_DOUBLE_PRECISION,isend,500000+isend,MPI_COMM_WORLD,ireqa(5,ii),IERR)
        Call MPI_Recv_init(dipoley(i),j,MPI_DOUBLE_PRECISION,isend,600000+isend,MPI_COMM_WORLD,ireqa(6,ii),IERR)
        Call MPI_Recv_init(dipolez(i),j,MPI_DOUBLE_PRECISION,isend,700000+isend,MPI_COMM_WORLD,ireqa(7,ii),IERR)
     endif
  Enddo

  if(myrank/=0) then
     ii=nrecv+1
     irecv=sendlist
     j=sendlist2
     n=j*natom
!       print*,myrank," send properties to ",irecv,j
     Call MPI_Send_init(fx(1,ista),   n,MPI_DOUBLE_PRECISION,irecv,200000+myrank,MPI_COMM_WORLD,ireqa(1,ii),IERR)
     Call MPI_Send_init(fy(1,ista),   n,MPI_DOUBLE_PRECISION,irecv,300000+myrank,MPI_COMM_WORLD,ireqa(2,ii),IERR)
     Call MPI_Send_init(fz(1,ista),   n,MPI_DOUBLE_PRECISION,irecv,400000+myrank,MPI_COMM_WORLD,ireqa(3,ii),IERR)
     Call MPI_Send_init(Eenergy(ista),j,MPI_DOUBLE_PRECISION,irecv,100000+myrank,MPI_COMM_WORLD,ireqa(4,ii),IERR)
     if(nocharge==0) then
      Call MPI_Send_init(charge(1,ista),n,MPI_DOUBLE_PRECISION,irecv,800000+myrank,MPI_COMM_WORLD,ireqa(8,ii),IERR)
     endif
     if(nohfcc==0) then
      Call MPI_Send_init(hfcc(1,ista),n,MPI_DOUBLE_PRECISION,irecv,900000+myrank,MPI_COMM_WORLD,ireqa(9,ii),IERR)
     endif

     if(nodipole==0) then
       Call MPI_Send_init(dipolex(ista),j,MPI_DOUBLE_PRECISION,irecv,500000+myrank,MPI_COMM_WORLD,ireqa(5,ii),IERR)
       Call MPI_Send_init(dipoley(ista),j,MPI_DOUBLE_PRECISION,irecv,600000+myrank,MPI_COMM_WORLD,ireqa(6,ii),IERR)
       Call MPI_Send_init(dipolez(ista),j,MPI_DOUBLE_PRECISION,irecv,700000+myrank,MPI_COMM_WORLD,ireqa(7,ii),IERR)
     endif

  endif

Return
End Subroutine
