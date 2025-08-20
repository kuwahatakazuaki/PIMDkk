subroutine Broad1
  use Parameters
#ifdef _mpi_
  use MPI
#endif
  implicit none
  integer :: ierr

#ifdef _mpi_
  call MPI_BCAST(Natom,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(temperature,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(Nbead,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(Nstep,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(dt,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(Isimulation,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(Nref,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(out_step,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(Nys,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(Nnhc,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(gamma1,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(Iforce,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  call MPI_Bcast(Iseeds,4,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(NCent,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(freq1,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(dir_result,80,MPI_CHARACTER,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(dir_scr,80,MPI_CHARACTER,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(version,80,MPI_CHARACTER,0,MPI_COMM_WORLD,IERR)

  call MPI_Bcast(Lsave_force,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
  call MPI_Bcast(Lsave_npa,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
  call MPI_Bcast(Lsave_charge,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
  call MPI_Bcast(Lsave_dipole,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
  call MPI_Bcast(Lsave_hfcc,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
  call MPI_Bcast(Lsave_energy,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
  call MPI_Bcast(Langstrom,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
  call MPI_Bcast(Lperiodic,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
  call MPI_Bcast(Lrestart,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
  call MPI_Bcast(Lrandom_coor,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)

  call MPI_BCAST(Iumb,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(umb_cons,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  !call MPI_BCAST(umb_pot,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(umb_atom1,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(umb_atom2,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(umb_atom3,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
#endif

  return
end subroutine Broad1

subroutine Broad2
  use Parameters
#ifdef _mpi_
  use MPI
  implicit none
  integer :: ierr
  real(8) :: tempr(3,Natom)

  Call MPI_BCAST(physmass,Natom,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  !Call MPI_BCAST(atom_num,Natom,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  Call MPI_BCAST(alabel,2*Natom,MPI_CHARACTER,0,MPI_COMM_WORLD,IERR)

  if ( MyRank == 0 ) then
    tempr(:,:) = ur(:,:,1)
  EndIf

  call MPI_Bcast(tempr,3*Natom,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)

  ur(:,:,1) = tempr(:,:)
#endif

  return
end subroutine Broad2


subroutine Broad3
  use Parameters
#ifdef _mpi_
  use MPI
  implicit none
  integer :: ierr

  call MPI_BCAST(Irestep,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  call MPI_Bcast(ur(1,1,1),3*Natom*Nbead,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  call MPI_Bcast(vur(1,1,1),3*Natom*Nbead,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  call MPI_Bcast(fur(1,1,1),3*Natom*Nbead,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)

  call MPI_Bcast(vrbath(1,1,1,1),3*Natom*Nnhc*Nbead,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  call MPI_Bcast(rbath(1,1,1,1),3*Natom*Nnhc*Nbead,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  call MPI_Bcast(frbath(1,1,1,1),3*Natom*Nnhc*Nbead,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)

  select case(Ncent)
    case(1)
      call MPI_Bcast(rbc11,Nnhc,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
      call MPI_Bcast(vbc11,Nnhc,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
      call MPI_Bcast(fbc11,Nnhc,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    case(3)
      call MPI_Bcast(frbc31(1,1,1),3*Natom*Nnhc,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
      call MPI_Bcast(vrbc31(1,1,1),3*Natom*Nnhc,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
      call MPI_Bcast(rbc31(1,1,1),3*Natom*Nnhc,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  end select
#endif

  return
end subroutine Broad3

