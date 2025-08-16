module mpi_module
  use Parameters
  use utility, only: program_abort
#ifdef _mpi_
  use mpi
  implicit none
  integer :: mstatus(MPI_STATUS_SIZE)
  integer :: ierr

contains

  subroutine MyMPI_gather_fr
    integer :: i, j, k, MyNbead
    integer :: Irank
    MyNbead = Iend - Ista + 1
    if ( MyRank == 0 ) then
      call MPI_Gather(MPI_IN_PLACE,3*Natom*MyNbead,MPI_DOUBLE_PRECISION, &
                      fr(1,1,Ista),3*Natom*MyNbead,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    else
      call MPI_Gather(fr(1,1,Ista),3*Natom*MyNbead,MPI_DOUBLE_PRECISION, &
                      fr(1,1,Ista),3*Natom*MyNbead,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    end if
    if ( MyRank == 0 ) then
      call MPI_Gather(MPI_IN_PLACE, MyNbead,MPI_DOUBLE_PRECISION, &
                      pot_bead(Ista),MyNbead,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    else
      call MPI_Gather(pot_bead(Ista),MyNbead,MPI_DOUBLE_PRECISION, &
                      pot_bead(Ista),MyNbead,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    end if
  end subroutine MyMPI_gather_fr

  subroutine MyMPI_gather_others
    integer :: i, j, k, MyNbead
    integer :: Irank
    MyNbead = Iend - Ista + 1

    if ( Lsave_charge .eqv. .True. ) then
      if ( MyRank == 0 ) then
        call MPI_Gather(MPI_IN_PLACE,  Natom*MyNbead,MPI_DOUBLE_PRECISION,&
                        charge(1,Ista),Natom*MyNbead,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      else
        call MPI_Gather(charge(1,Ista),Natom*MyNbead,MPI_DOUBLE_PRECISION,&
                        charge(1,Ista),Natom*MyNbead,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      end if
    end if

    if ( Lsave_dipole .eqv. .True. ) then
      if ( MyRank == 0 ) then
        call MPI_Gather(MPI_IN_PLACE,   3*MyNbead,MPI_DOUBLE_PRECISION,&
                        dipoler(1,Ista),3*MyNbead,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      else
        call MPI_Gather(dipoler(1,Ista),3*MyNbead,MPI_DOUBLE_PRECISION,&
                        dipoler(1,Ista),3*MyNbead,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      end if
    end if
!print *, MyRank, "dipoler"

    if ( Lsave_hfcc .eqv. .True. ) then
      if ( MyRank == 0 ) then
        call MPI_Gather(MPI_IN_PLACE,Natom*MyNbead,MPI_DOUBLE_PRECISION,&
                        hfcc(1,Ista),Natom*MyNbead,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      else
        call MPI_Gather(hfcc(1,Ista),Natom*MyNbead,MPI_DOUBLE_PRECISION,&
                        hfcc(1,Ista),Natom*MyNbead,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      end if
    end if

  end subroutine MyMPI_gather_others

  subroutine MyMPI_scatter_r
    integer :: i, j, k, MyNbead
    MyNbead = Iend - Ista + 1
    if ( MyRank == 0 ) then
      call MPI_Scatter(r(1,1,Ista),3*Natom*MyNbead,MPI_DOUBLE_PRECISION, &
                       MPI_IN_PLACE,3*Natom*MyNbead,MPI_DOUBLE_PRECISION, &
                       0,MPI_COMM_WORLD,ierr)
    else
      call MPI_Scatter(r(1,1,Ista),3*Natom*MyNbead,MPI_DOUBLE_PRECISION, &
                       r(1,1,Ista),3*Natom*MyNbead,MPI_DOUBLE_PRECISION, &
                       0,MPI_COMM_WORLD,ierr)
    end if
  end subroutine MyMPI_scatter_r

  subroutine MyMPI_bcast_r
    call MPI_Bcast(r(1,1,1),3*Natom*Nbead,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,Ierr)
  end subroutine MyMPI_bcast_r

  !subroutine MyMPI_gather_r
  !  integer :: i, j, k, MyNbead
  !  MyNbead = Iend - Ista + 1
  !  if ( MyRank == 0 ) then
  !    call MPI_Gather(MPI_IN_PLACE,3*Natom*MyNbead,MPI_DOUBLE_PRECISION, &
  !                    r(1,1,Ista), 3*Natom*MyNbead,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  !  else
  !    call MPI_Gather(r(1,1,Ista),3*Natom*MyNbead,MPI_DOUBLE_PRECISION, &
  !                    r(1,1,Ista),3*Natom*MyNbead,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  !  end if
  !end subroutine MyMPI_gather_r

  !subroutine MyMPI_allgather_r
  !  integer :: i, j, k, MyNbead
  !  integer :: Irank
  !  MyNbead = Iend - Ista + 1
  !  call MPI_Allgather(r(1,1,Ista),3*Natom*MyNbead,MPI_DOUBLE_PRECISION, &
  !                     r(1,1,Ista),3*Natom*MyNbead,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
  !end subroutine MyMPI_allgather_r

  subroutine para_range(N1,N2,Nproc,Irank,Ista,Iend)
    integer, intent(in) :: N1, N2, Nproc, Irank
    integer, intent(out) :: Ista, Iend
    integer :: Itemp1, Itemp2
    Itemp1 = (N2-N1+1) / Nproc
    Itemp2 = mod(N2-N1+1,Nproc)
    Ista  = Irank * Itemp1 + N1 + min(Irank,Itemp2)
    Iend = Ista + Itemp1 - 1
    if (Itemp2 > Irank) Iend = Iend + 1
  end subroutine para_range

#endif

end module mpi_module

