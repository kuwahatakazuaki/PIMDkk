subroutine set_pallarel
  use Parameters

#ifdef _mpi_
  use MPI
  implicit none
  integer   :: i,j,k
  integer :: Nhmod, numeach

  Allocate (listeach(Nproc))
  Allocate (listeachtmp(Nproc))

  Nhmod=mod(Nbead,Nproc)
  numeach=(Nbead-Nhmod)/Nproc

  do i=1,Nhmod
    listeach(i)=numeach+1
  end do

  do i=Nhmod+1,Nproc
    listeach(i)=numeach
  end do
  numeach=listeach(myrank+1)

  ista=1
  do i=1,myrank
    ista=ista+listeach(i)
  end do
  iend=ista+listeach(myrank+1)-1

#else

  implicit none
  ista = 1
  iend = Nbead

#endif

  return
end subroutine set_pallarel

