  Subroutine Start_Send_Recv_MPI_tk

    Use Parameters
    Use Parameter_tk
    Use MPI
    Implicit None
    Integer :: i,j,k,ii,jj,isend,irecv,n

!    print*,myrank," test ",nrecv
!    wait(i)
    Do ii=1,nrecv
       Call MPI_Start(ireqa(1,ii),IERR)
       Call MPI_Start(ireqa(2,ii),IERR)
       Call MPI_Start(ireqa(3,ii),IERR)
    Enddo
    Do i=1,nrecv
       Call MPI_Wait(ireqa(1,i),mstatus,IERR)
       Call MPI_Wait(ireqa(2,i),mstatus,IERR)
       Call MPI_Wait(ireqa(3,i),mstatus,IERR)
    Enddo
    if(myrank/=0) then
       ii=nrecv+1
       Call MPI_Start(ireqa(1,ii),IERR)
       Call MPI_Start(ireqa(2,ii),IERR)
       Call MPI_Start(ireqa(3,ii),IERR)
       Call MPI_Wait(ireqa(1,ii),mstatus,IERR)
       Call MPI_Wait(ireqa(2,ii),mstatus,IERR)
       Call MPI_Wait(ireqa(3,ii),mstatus,IERR)
    endif

    Do ii=1,nrecv
       Call MPI_Start(ireqa(4,ii),IERR)
     if(nocharge==0) then
       Call MPI_Start(ireqa(8,ii),IERR)
     endif
     if(nohfcc==0) then
       Call MPI_Start(ireqa(9,ii),IERR)
     endif
    Enddo

    Do i=1,nrecv
         Call MPI_Wait(ireqa(4,i),mstatus,IERR)
       if(nocharge==0) then
         Call MPI_Wait(ireqa(8,i),mstatus,IERR)
       endif
       if(nohfcc==0) then
         Call MPI_Wait(ireqa(9,i),mstatus,IERR)
       endif
    Enddo

    if(myrank/=0) then
       ii=nrecv+1
       Call MPI_Start(ireqa(4,ii),IERR)
      if(nocharge==0) then
       Call MPI_Start(ireqa(8,ii),IERR)
      endif
      if(nohfcc==0) then
       Call MPI_Start(ireqa(9,ii),IERR)
      endif
       Call MPI_Wait(ireqa(4,ii),mstatus,IERR)
      if(nocharge==0) then
       Call MPI_Wait(ireqa(8,ii),mstatus,IERR)
      endif
      if(nohfcc==0) then
       Call MPI_Wait(ireqa(9,ii),mstatus,IERR)
      endif
    endif

    if(nodipole==0) then
      Do ii=1,nrecv
         Call MPI_Start(ireqa(5,ii),IERR)
         Call MPI_Start(ireqa(6,ii),IERR)
         Call MPI_Start(ireqa(7,ii),IERR)
      Enddo
      Do i=1,nrecv
         Call MPI_Wait(ireqa(5,i),mstatus,IERR)
         Call MPI_Wait(ireqa(6,i),mstatus,IERR)
         Call MPI_Wait(ireqa(7,i),mstatus,IERR)
      Enddo
      if(myrank/=0) then
         ii=nrecv+1
         Call MPI_Start(ireqa(5,ii),IERR)
         Call MPI_Start(ireqa(6,ii),IERR)
         Call MPI_Start(ireqa(7,ii),IERR)
         Call MPI_Wait(ireqa(5,ii),mstatus,IERR)
         Call MPI_Wait(ireqa(6,ii),mstatus,IERR)
         Call MPI_Wait(ireqa(7,ii),mstatus,IERR)
      endif
   endif

Return
End Subroutine
