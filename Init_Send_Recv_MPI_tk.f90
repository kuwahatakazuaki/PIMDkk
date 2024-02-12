Subroutine Init_Send_Recv_MPI_tk
  Use Parameters
  !Use Parameter_tk
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
!      print*,myrank," recv properties from ",isend,i,j
    !Call MPI_Recv_init(fx(1,i),n,MPI_DOUBLE_PRECISION,isend,20000+isend,MPI_COMM_WORLD,ireqa(1,ii),IERR)
    !Call MPI_Recv_init(fy(1,i),n,MPI_DOUBLE_PRECISION,isend,30000+isend,MPI_COMM_WORLD,ireqa(2,ii),IERR)
    !Call MPI_Recv_init(fz(1,i),n,MPI_DOUBLE_PRECISION,isend,40000+isend,MPI_COMM_WORLD,ireqa(3,ii),IERR)
    Call MPI_Recv_init(fr(1,1,i),3*n,MPI_DOUBLE_PRECISION,isend,20000+isend,MPI_COMM_WORLD,ireqa(1,ii),IERR)
    Call MPI_Recv_init(Eenergy(i),j,MPI_DOUBLE_PRECISION,isend,10000+isend,MPI_COMM_WORLD,ireqa(4,ii),IERR)
    if( Lsave_charge .eqv. .True.) then
      Call MPI_Recv_init(charge(1,i),n,MPI_DOUBLE_PRECISION,isend,80000+isend,MPI_COMM_WORLD,ireqa(8,ii),IERR)
    endif
    if( Lsave_hfcc .eqv. .True. ) then
      Call MPI_Recv_init(hfcc(1,i),n,MPI_DOUBLE_PRECISION,isend,90000+isend,MPI_COMM_WORLD,ireqa(9,ii),IERR)
    endif

    if( Lsave_dipole .eqv. .True. ) then
        !Call MPI_Recv_init(dipolex(i),j,MPI_DOUBLE_PRECISION,isend,50000+isend,MPI_COMM_WORLD,ireqa(5,ii),IERR)
        !Call MPI_Recv_init(dipoley(i),j,MPI_DOUBLE_PRECISION,isend,60000+isend,MPI_COMM_WORLD,ireqa(6,ii),IERR)
        !Call MPI_Recv_init(dipolez(i),j,MPI_DOUBLE_PRECISION,isend,70000+isend,MPI_COMM_WORLD,ireqa(7,ii),IERR)
        Call MPI_Recv_init(dipoler(1,i),3*j,MPI_DOUBLE_PRECISION,isend,100000+isend,MPI_COMM_WORLD,ireqa(10,ii),IERR)
     endif

!     if ( Nforce == 8 ) then
!        Call MPI_Recv_init(pressure(i),j,MPI_DOUBLE_PRECISION,isend,110000+isend,MPI_COMM_WORLD,ireqa(11,ii),IERR)
!     end if
  Enddo

  if(myrank/=0) then
     ii=nrecv+1
     irecv=sendlist
     j=sendlist2
     n=j*natom
!       print*,myrank," send properties to ",irecv,j
     !Call MPI_Send_init(fx(1,ista),n,MPI_DOUBLE_PRECISION,irecv,20000+myrank,MPI_COMM_WORLD,ireqa(1,ii),IERR)
     !Call MPI_Send_init(fy(1,ista),n,MPI_DOUBLE_PRECISION,irecv,30000+myrank,MPI_COMM_WORLD,ireqa(2,ii),IERR)
     !Call MPI_Send_init(fz(1,ista),n,MPI_DOUBLE_PRECISION,irecv,40000+myrank,MPI_COMM_WORLD,ireqa(3,ii),IERR)
     Call MPI_Send_init(fr(1,1,ista),3*n,MPI_DOUBLE_PRECISION,irecv,20000+myrank,MPI_COMM_WORLD,ireqa(1,ii),IERR)
     Call MPI_Send_init(Eenergy(ista),j,MPI_DOUBLE_PRECISION,irecv,10000+myrank,MPI_COMM_WORLD,ireqa(4,ii),IERR)
     if( Lsave_charge .eqv. .True.) then
      Call MPI_Send_init(charge(1,ista),n,MPI_DOUBLE_PRECISION,irecv,80000+myrank,MPI_COMM_WORLD,ireqa(8,ii),IERR)
     endif
     if( Lsave_hfcc .eqv. .True. ) then
      Call MPI_Send_init(hfcc(1,ista),n,MPI_DOUBLE_PRECISION,irecv,90000+myrank,MPI_COMM_WORLD,ireqa(9,ii),IERR)
     endif

     if( Lsave_dipole .eqv. .True. ) then
       !Call MPI_Send_init(dipolex(ista),j,MPI_DOUBLE_PRECISION,irecv,50000+myrank,MPI_COMM_WORLD,ireqa(5,ii),IERR)
       !Call MPI_Send_init(dipoley(ista),j,MPI_DOUBLE_PRECISION,irecv,60000+myrank,MPI_COMM_WORLD,ireqa(6,ii),IERR)
       !Call MPI_Send_init(dipolez(ista),j,MPI_DOUBLE_PRECISION,irecv,70000+myrank,MPI_COMM_WORLD,ireqa(7,ii),IERR)
       Call MPI_Send_init(dipoler(1,ista),3*j,MPI_DOUBLE_PRECISION,irecv,100000+myrank,MPI_COMM_WORLD,ireqa(10,ii),IERR)
     endif

!     if ( Nforce == 8 ) then
!        Call MPI_Send_init(pressure(ista),j,MPI_DOUBLE_PRECISION,irecv,110000+myrank,MPI_COMM_WORLD,ireqa(11,ii),IERR)
!     end if

  endif

Return
End Subroutine
