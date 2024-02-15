Subroutine Start_Send_Recv_MPI_tk
  Use Parameters
  Use MPI
  Implicit None
  Integer :: i,j,k,ii,jj,isend,irecv,n
  Integer :: mstatus(MPI_STATUS_SIZE)
  integer :: ierr
  real(8) :: tmp

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
   !if(nocharge==0) then
   if( Lsave_charge .eqv. .True.) then
     Call MPI_Start(ireqa(8,ii),IERR)
   endif
   !if(nohfcc==0) then
   if( Lsave_hfcc .eqv. .True. ) then
     Call MPI_Start(ireqa(9,ii),IERR)
   endif
  Enddo

  Do i=1,nrecv
       Call MPI_Wait(ireqa(4,i),mstatus,IERR)
     !if(nocharge==0) then
     if( Lsave_charge .eqv. .True.) then
       Call MPI_Wait(ireqa(8,i),mstatus,IERR)
     endif
     !if(nohfcc==0) then
     if( Lsave_hfcc .eqv. .True. ) then
       Call MPI_Wait(ireqa(9,i),mstatus,IERR)
     endif
  Enddo

  if(myrank/=0) then
     ii=nrecv+1
     Call MPI_Start(ireqa(4,ii),IERR)
    !if(nocharge==0) then
    if( Lsave_charge .eqv. .True.) then
     Call MPI_Start(ireqa(8,ii),IERR)
    endif
    !if(nohfcc==0) then
    if( Lsave_hfcc .eqv. .True. ) then
     Call MPI_Start(ireqa(9,ii),IERR)
    endif
     Call MPI_Wait(ireqa(4,ii),mstatus,IERR)
    !if(nocharge==0) then
    if( Lsave_charge .eqv. .True.) then
     Call MPI_Wait(ireqa(8,ii),mstatus,IERR)
    endif
    !if(nohfcc==0) then
    if( Lsave_hfcc .eqv. .True. ) then
     Call MPI_Wait(ireqa(9,ii),mstatus,IERR)
    endif
  endif

  if ( Iforce == 8 ) then
    if ( MyRank == 0 ) then
    call MPI_Gather(MPI_IN_PLACE,listeach(MyRank+1),MPI_DOUBLE_PRECISION, & 
                    pressure(ista),listeach(MyRank+1),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    else
    call MPI_Gather(pressure(ista),listeach(MyRank+1),MPI_DOUBLE_PRECISION, & 
                    pressure(ista),listeach(MyRank+1),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    end if
    !tmp = pressure(ista)
    !call MPI_Gather(tmp,listeach(MyRank+1),MPI_DOUBLE_PRECISION, & 
    !                pressure(ista),listeach(MyRank+1),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  end if

  !if(nodipole==0) then
  if( Lsave_dipole .eqv. .True. ) then
    Do ii=1,nrecv
       Call MPI_Start(ireqa(5,ii),IERR)
       Call MPI_Start(ireqa(6,ii),IERR)
       Call MPI_Start(ireqa(7,ii),IERR)
       Call MPI_Start(ireqa(10,ii),IERR)
    Enddo
    Do i=1,nrecv
       Call MPI_Wait(ireqa(5,i),mstatus,IERR)
       Call MPI_Wait(ireqa(6,i),mstatus,IERR)
       Call MPI_Wait(ireqa(7,i),mstatus,IERR)
       Call MPI_Wait(ireqa(10,i),mstatus,IERR)
    Enddo
    if(myrank/=0) then
       ii=nrecv+1
       Call MPI_Start(ireqa(5,ii),IERR)
       Call MPI_Start(ireqa(6,ii),IERR)
       Call MPI_Start(ireqa(7,ii),IERR)
       Call MPI_Start(ireqa(10,ii),IERR)
       Call MPI_Wait(ireqa(5,ii),mstatus,IERR)
       Call MPI_Wait(ireqa(6,ii),mstatus,IERR)
       Call MPI_Wait(ireqa(7,ii),mstatus,IERR)
       Call MPI_Wait(ireqa(10,ii),mstatus,IERR)
    endif
  endif

Return
End Subroutine
