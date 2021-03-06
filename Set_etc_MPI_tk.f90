Subroutine Set_etc_MPI_tk

  Use MPI
  Use Parameters
  Use Parameter_tk
  Implicit None
  Integer   :: i,j,k,n,n0


  Allocate (listeach(nprocs))
  Allocate (listeachtmp(nprocs))

    nhmod=mod(nbead,nprocs)
    numeach=(nbead-nhmod)/nprocs

    Do i=1,nhmod
      listeach(i)=numeach+1
    Enddo
    Do i=nhmod+1,nprocs
      listeach(i)=numeach
    Enddo
    numeach=listeach(myrank+1)
    ista=1
    Do i=1,myrank
      ista=ista+listeach(i)
    Enddo
    iend=ista+listeach(myrank+1)-1

    na3=3*natom
    na31=na3+1
!    If (myrank==0) Allocate(work2(nbead*na31)))
!    Allocate(work(numeach*na31))
    npacksize=(natom*4)+4
    Allocate(work((nbead-ista+1)*npacksize))


    j=nprocs
    Do i=1,1000
     j=(j+1)/2
     if(j==1) then
       nsendrecv=i
       Exit
     endif
    Enddo

    Allocate (recvlist(nsendrecv))
    Allocate (recvlist1(nsendrecv))
    Allocate (recvlist2(nsendrecv))

    nrecv=0
    n=1
    sendlist=0
    sendlist1=ista
    sendlist2=0
    Do i=myrank+1,nprocs
      listeachtmp(i)=listeach(i)
    Enddo
    Do i=1,nsendrecv
      n0=n
      n=n0*2
      if((mod(myrank,n)==0) .AND. (myrank+n0<nprocs)) then
        nrecv=nrecv+1
        recvlist(nrecv)=myrank+n0
        recvlist1(nrecv)=ista+listeachtmp(myrank+1)
        recvlist2(nrecv)=listeachtmp(myrank+n0+1)
!        print*,myrank," to ",recvlist(nrecv),recvlist1(nrecv),recvlist2(nrecv)
      else if(mod(myrank,n)==n0) then
        sendlist=myrank-n0
        sendlist1=ista
        sendlist2=listeachtmp(myrank+1)
!        print*,myrank," from ",sendlist,sendlist1,sendlist2
      endif
      Do j=myrank+1,nprocs-n0,n
         listeachtmp(j)=listeachtmp(j)+listeachtmp(j+n0)
      Enddo
    Enddo
!    Do i=1,nrecv
!        recvlistx(i)=recvlist(nrecv-i+1)
!        recvlistx1(i)=recvlist1(nrecv-i+1)
!        recvlistx2(i)=recvlist2(nrecv-i+1)
!    Enddo


!for MPI_Send_init ver
   Allocate(ireqa(9,nrecv+1))
   Allocate(ireqb(3,nrecv+1))
   Call Init_Send_Recv_MPI_tk
   Call Init_Recv_Send_MPI_tk

  Return
End Subroutine

