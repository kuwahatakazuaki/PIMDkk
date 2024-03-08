module mpi_module
  implicit none
#ifdef _mpi_
  use mpi
#endif

contains
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

!  subroutine program_abort(message)
!#ifdef _mpi_
!    use mpi
!#endif
!    character(*) :: message
!    integer :: ierr
!    print *, message
!#ifdef _mpi_
!    call mpi_abort(MPI_COMM_WORLD, -1, ierr)
!#else
!    stop
!#endif
!  end subroutine program_abort

end module mpi_module

