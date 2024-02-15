  Subroutine Send_Recv_MPI_tk3

    Use Parameters
    Use Parameter_tk
    Use MPI
    Implicit None

    Integer :: i,j,k,ii,jj,isend,irecv,n,i0

    n=1
    Do ii=myrank+1,nprocs
      listeachtmp(ii)=listeach(ii)
    Enddo
    i0=0
    Do ii=ista,iend
      work(i0+1)=Eenergy(ii)
      work(i0+2)=dipolex(ii)
      work(i0+3)=dipoley(ii)
      work(i0+4)=dipolez(ii)
      i0=i0+4
      Do jj=1,natom
       work(jj+i0)=fx(jj,ii)
      Enddo
      i0=i0+natom
      Do jj=1,natom
       work(jj+i0)=fy(jj,ii)
      Enddo
      i0=i0+natom
      Do jj=1,natom
       work(jj+i0)=fz(jj,ii)
      Enddo
      i0=i0+natom
      Do jj=1,natom
       work(jj+i0)=charge(jj,ii)
      Enddo
      i0=i0+natom
    Enddo
    Do ii=1,nsendrecv
      n=n*2
      if(mod(myrank,n)==0) then
       isend=myrank+n-1
       j=listeachtmp(isend+1)*npacksize
       i=listeachtmp(myrank+1)*npacksize+1
       Call MPI_Irecv(work(i),j,MPI_DOUBLE_PRECISION,isend,100000+isend,MPI_COMM_WORLD,ireq0,IERR)
       Do jj=myrank+1,nsendrecv,n
         listeachtmp(jj)=listeachtmp(jj)+listeachtmp(jj+n-1)
       Enddo
       Call MPI_Wait(ireq0,mstatus,IERR)
      else if(mod(myrank,n)==n-1) then
       irecv=myrank-n+1
       j=listeachtmp(myrank+1)*npacksize
       Call MPI_Isend(work,j,MPI_DOUBLE_PRECISION,irecv,100000+myrank,MPI_COMM_WORLD,ireq0,IERR)
       Call MPI_Wait(ireq0,mstatus,IERR)
      endif
    Enddo

   If(myrank==0) then
    i0=0
    Do ii=1,nbead
      Eenergy(ii)=work(i0+1)
      dipolex(ii)=work(i0+2)
      dipoley(ii)=work(i0+3)
      dipolez(ii)=work(i0+4)
      i0=i0+4
      Do jj=1,natom
       fx(jj,ii)=work(jj+i0)
      Enddo
      i0=i0+natom
      Do jj=1,natom
       fy(jj,ii)=work(jj+i0)
      Enddo
      i0=i0+natom
      Do jj=1,natom
       fz(jj,ii)=work(jj+i0)
      Enddo
      i0=i0+natom
      Do jj=1,natom
       charge(jj,ii)=work(jj+i0)
      Enddo
      i0=i0+natom
    Enddo
   Endif

    Return
  End Subroutine
