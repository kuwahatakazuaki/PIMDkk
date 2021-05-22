  Subroutine Send_Recv_MPI_tk2

    Use Parameters
    Use Parameter_tk
    Use MPI
    Implicit None

    Integer :: i,j,k,ii,jj,isend,irecv,n,n0

    n=1
    Do ii=myrank+1,nprocs
      listeachtmp(ii)=listeach(ii)
    Enddo
    Do ii=1,nsendrecv
      n0=n
      n=n0*2
      if((mod(myrank,n)==0) .AND. (myrank+n0<nprocs)) then
       isend=myrank+n0
       j=listeachtmp(isend+1)
       i=ista+listeachtmp(myrank+1)
!       print*,myrank," recv properties from ",isend,i,j
       Call MPI_Irecv(fx(1,i),j*natom,MPI_DOUBLE_PRECISION,isend,200000+isend,MPI_COMM_WORLD,ireq(1),IERR)
       Call MPI_Irecv(fy(1,i),j*natom,MPI_DOUBLE_PRECISION,isend,300000+isend,MPI_COMM_WORLD,ireq(2),IERR)
       Call MPI_Irecv(fz(1,i),j*natom,MPI_DOUBLE_PRECISION,isend,400000+isend,MPI_COMM_WORLD,ireq(3),IERR)
       Call MPI_Irecv(Eenergy(i),j,MPI_DOUBLE_PRECISION,isend,100000+isend,MPI_COMM_WORLD,ireq(4),IERR)
!      if(nodipole==0) then
         Call MPI_Irecv(dipolex(i),j,MPI_DOUBLE_PRECISION,isend,500000+isend,MPI_COMM_WORLD,ireq(5),IERR)
         Call MPI_Irecv(dipoley(i),j,MPI_DOUBLE_PRECISION,isend,600000+isend,MPI_COMM_WORLD,ireq(6),IERR)
         Call MPI_Irecv(dipolez(i),j,MPI_DOUBLE_PRECISION,isend,700000+isend,MPI_COMM_WORLD,ireq(7),IERR)
!       endif
!       if(nocharge==0) then
         Call MPI_Irecv(charge(1,i),j*natom,MPI_DOUBLE_PRECISION,isend,800000+isend,MPI_COMM_WORLD,ireq(8),IERR)
!      endif
       Do jj=myrank+1,nprocs-n0,n
         listeachtmp(jj)=listeachtmp(jj)+listeachtmp(jj+n0)
       Enddo
       Call MPI_Wait(ireq(1),mstatus,IERR)
       Call MPI_Wait(ireq(2),mstatus,IERR)
       Call MPI_Wait(ireq(3),mstatus,IERR)
       Call MPI_Wait(ireq(4),mstatus,IERR)
       Call MPI_Wait(ireq(5),mstatus,IERR)
       Call MPI_Wait(ireq(6),mstatus,IERR)
       Call MPI_Wait(ireq(7),mstatus,IERR)
       Call MPI_Wait(ireq(8),mstatus,IERR)
      else if(mod(myrank,n)==n0) then
       irecv=myrank-n0
       j=listeachtmp(myrank+1)
!       print*,myrank," send properties to ",irecv,j
       Call MPI_Isend(fx(1,ista),j*natom,MPI_DOUBLE_PRECISION,irecv,200000+myrank,MPI_COMM_WORLD,ireq(1),IERR)
       Call MPI_Isend(fy(1,ista),j*natom,MPI_DOUBLE_PRECISION,irecv,300000+myrank,MPI_COMM_WORLD,ireq(2),IERR)
       Call MPI_Isend(fz(1,ista),j*natom,MPI_DOUBLE_PRECISION,irecv,400000+myrank,MPI_COMM_WORLD,ireq(3),IERR)
       Call MPI_Isend(Eenergy(ista),j,MPI_DOUBLE_PRECISION,irecv,100000+myrank,MPI_COMM_WORLD,ireq(4),IERR)
!       if(nodipole==0) then
        Call MPI_Isend(dipolex(ista),j,MPI_DOUBLE_PRECISION,irecv,500000+myrank,MPI_COMM_WORLD,ireq(5),IERR)
        Call MPI_Isend(dipoley(ista),j,MPI_DOUBLE_PRECISION,irecv,600000+myrank,MPI_COMM_WORLD,ireq(6),IERR)
        Call MPI_Isend(dipolez(ista),j,MPI_DOUBLE_PRECISION,irecv,700000+myrank,MPI_COMM_WORLD,ireq(7),IERR)
!       endif
!       if(nocharge==0) then
        Call MPI_Isend(charge(1,ista),j*natom,MPI_DOUBLE_PRECISION,irecv,800000+myrank,MPI_COMM_WORLD,ireq(8),IERR)
!       endif
       Call MPI_Wait(ireq(1),mstatus,IERR)
       Call MPI_Wait(ireq(2),mstatus,IERR)
       Call MPI_Wait(ireq(3),mstatus,IERR)
       Call MPI_Wait(ireq(4),mstatus,IERR)
       Call MPI_Wait(ireq(5),mstatus,IERR)
       Call MPI_Wait(ireq(6),mstatus,IERR)
       Call MPI_Wait(ireq(7),mstatus,IERR)
       Call MPI_Wait(ireq(8),mstatus,IERR)
      endif
    Enddo


    Return
  End Subroutine
