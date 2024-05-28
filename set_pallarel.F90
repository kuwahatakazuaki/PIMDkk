subroutine set_pallarel
  use Parameters

#ifdef _mpi_
  use MPI
  implicit none
  integer   :: i,j,k,n,n0
  integer :: Nhmod, numeach

  Allocate (listeach(Nproc))
  Allocate (listeachtmp(Nproc))

  Nhmod=mod(Nbead,Nproc)
  numeach=(Nbead-Nhmod)/Nproc

  Do i=1,Nhmod
    listeach(i)=numeach+1
  Enddo
  Do i=Nhmod+1,Nproc
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

  return
end subroutine set_pallarel

