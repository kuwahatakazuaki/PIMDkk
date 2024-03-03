subroutine set_pallarel
!subroutine Set_etc_MPI_tk
  use Parameters

#ifdef _mpi_
  use MPI
  implicit none
  integer   :: i,j,k,n,n0
  integer :: Nhmod, numeach

  Allocate (listeach(nprocs))
  Allocate (listeachtmp(nprocs))

  Nhmod=mod(Nbead,nprocs)
  numeach=(Nbead-Nhmod)/nprocs

  Do i=1,Nhmod
    listeach(i)=numeach+1
  Enddo
  Do i=Nhmod+1,nprocs
    listeach(i)=numeach
  Enddo
  numeach=listeach(myrank+1)
  ista=1
  Do i=1,myrank
    ista=ista+listeach(i)
  Enddo
  iend=ista+listeach(myrank+1)-1
#else
  implicit none
  ista = 1
  iend = Nbead
#endif

! === till here ===

#ifdef _mpi_
  na3=3*natom
  na31=na3+1
  npacksize=(natom*4)+4
  Allocate(work((Nbead-ista+1)*npacksize))

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
!      print*,myrank," to ",recvlist(nrecv),recvlist1(nrecv),recvlist2(nrecv)
    else if(mod(myrank,n)==n0) then
      sendlist=myrank-n0
      sendlist1=ista
      sendlist2=listeachtmp(myrank+1)
!      print*,myrank," from ",sendlist,sendlist1,sendlist2
    endif
    Do j=myrank+1,nprocs-n0,n
       listeachtmp(j)=listeachtmp(j)+listeachtmp(j+n0)
    Enddo
  Enddo

!for MPI_Send_init ver
  Allocate(ireqa(11,nrecv+1))
  Allocate(ireqb(3,nrecv+1))
  call Init_Send_Recv_MPI_tk
  call Init_Recv_Send_MPI_tk
#endif

  return
end subroutine set_pallarel

