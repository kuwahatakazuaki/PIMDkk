  Subroutine Recv_Send_MPI_tk3
    Use Parameters
    Use MPI

    Implicit None
    Integer :: i,j,k,ii,jj,isend,irecv,n,n0,i0
    Integer :: mstatus(MPI_STATUS_SIZE)

   if(myrank==0) then
    i0=numeach*na3+1
    Do ii=numeach+1,nbead
      Do jj=1,natom
       work(i0)=x(jj,ii)
       work(i0+1)=y(jj,ii)
       work(i0+2)=z(jj,ii)
       i0=i0+3
      Enddo
    Enddo
   endif

    n0=2**nsendrecv
    k=nprocs
    Do ii=1,nsendrecv
      n=n0
      n0=n/2
      if(mod(myrank,n)==0) then
       isend=myrank+n0
       j=0
       i=1
       Do jj=myrank+1,isend
         i=i+listeach(jj)*na3
       Enddo
       Do jj=isend+1,k
         j=j+listeach(jj)*na3
       Enddo
!       print*,myrank," send to ",isend,i,j
       Call MPI_Isend(work(i),j,MPI_DOUBLE_PRECISION,isend,900000+isend,MPI_COMM_WORLD,ireq0,IERR)
       k=isend
       Call MPI_Wait(ireq0,mstatus,IERR)
      else if(mod(myrank,n)==n0) then
       irecv=myrank-n0
       j=0
       Do jj=myrank+1,myrank+n0
         j=j+listeach(jj)*na3
       Enddo
!       print*,myrank," recv from ",irecv,j
       Call MPI_Irecv(work,j,MPI_DOUBLE_PRECISION,irecv,900000+myrank,MPI_COMM_WORLD,ireq0,IERR)
       k=myrank+n0
       Call MPI_Wait(ireq0,mstatus,IERR)
      endif
!      Call MPI_BARRIER(MPI_COMM_WORLD,IERR)
    Enddo

   if(myrank/=0) then
    i0=1
    Do ii=ista,iend
      Do jj=1,natom
       x(jj,ii)=work(i0)
       y(jj,ii)=work(i0+1)
       z(jj,ii)=work(i0+2)
       i0=i0+3
      Enddo
    Enddo
   endif

    Return
  End Subroutine
