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
  call MPI_BCAST(Nys,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(Nnhc,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(gamma,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(Iforce,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  call MPI_Bcast(Iseeds,4,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(NCent,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(freq1,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(address,80,MPI_CHARACTER,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(address2,80,MPI_CHARACTER,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(version,80,MPI_CHARACTER,0,MPI_COMM_WORLD,IERR)
  call MPI_Bcast(Lsave_force,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
  call MPI_Bcast(Lsave_npa,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
  call MPI_Bcast(Lsave_charge,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
  call MPI_Bcast(Lsave_dipole,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
  call MPI_Bcast(Lsave_hfcc,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
  call MPI_Bcast(Langstrom,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
  call MPI_Bcast(Lrestart,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
  call MPI_Bcast(Lrandom_coor,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(umbrella_sampling,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(umbrella_width,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(umbrella_height,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(umbrella_constant,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(umbrella_atom1,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(umbrella_atom2,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(umbrella_atom3,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
#endif

  return
End Subroutine

subroutine Broad2
  use Parameters
#ifdef _mpi_
  use MPI
  implicit none
  integer :: ierr
  real(8) :: tempr(3,Natom)

  Call MPI_BCAST(physmass,natom,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  Call MPI_BCAST(atom_num,natom,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  Call MPI_BCAST(alabel,2*natom,MPI_CHARACTER,0,MPI_COMM_WORLD,IERR)

  if ( MyRank == 0 ) then
    !x1(:)=ux(:,1)
    !y1(:)=uy(:,1)
    !z1(:)=uz(:,1)
    tempr(:,:) = ur(:,:,1)
  EndIf

  !Call MPI_BCAST(x1,natom,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  !Call MPI_BCAST(y1,natom,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  !Call MPI_BCAST(z1,natom,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  call MPI_Bcast(tempr,3*natom,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)

  !ux(:,1)=x1(:)
  !uy(:,1)=y1(:)
  !uz(:,1)=z1(:)
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

  Call MPI_BCAST(Nrstep,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  !Call MPI_BCAST(ux(1,1),natom*nbead,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  !Call MPI_BCAST(uy(1,1),natom*nbead,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  !Call MPI_BCAST(uz(1,1),natom*nbead,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  call MPI_Bcast(ur(1,1,1),3*natom*nbead,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  !Call MPI_BCAST(vux(1,1),natom*nbead,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  !Call MPI_BCAST(vuy(1,1),natom*nbead,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  !Call MPI_BCAST(vuz(1,1),natom*nbead,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  call MPI_Bcast(vur(1,1,1),3*natom*nbead,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  !Call MPI_BCAST(fux(1,1),natom*nbead,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  !Call MPI_BCAST(fuy(1,1),natom*nbead,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  !Call MPI_BCAST(fuz(1,1),natom*nbead,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  call MPI_Bcast(fur(1,1,1),3*natom*nbead,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)

  !Call MPI_BCAST(vxbath(1,1,1),natom*nnhc*nbead,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  !Call MPI_BCAST(vybath(1,1,1),natom*nnhc*nbead,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  !Call MPI_BCAST(vzbath(1,1,1),natom*nnhc*nbead,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  call MPI_Bcast(vrbath(1,1,1,1),3*natom*nnhc*nbead,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  !Call MPI_BCAST(xbath(1,1,1),natom*nnhc*nbead,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  !Call MPI_BCAST(ybath(1,1,1),natom*nnhc*nbead,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  !Call MPI_BCAST(zbath(1,1,1),natom*nnhc*nbead,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  call MPI_Bcast(rbath(1,1,1,1),3*natom*nnhc*nbead,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  !Call MPI_BCAST(fxbath(1,1,1),natom*nnhc*nbead,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  !Call MPI_BCAST(fybath(1,1,1),natom*nnhc*nbead,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  !Call MPI_BCAST(fzbath(1,1,1),natom*nnhc*nbead,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  call MPI_Bcast(frbath(1,1,1,1),3*natom*nnhc*nbead,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)

  select case(Ncent)
    case(1)
  !If(Ncent==1) Then
      call MPI_Bcast(rbc11,nnhc,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
      call MPI_Bcast(vbc11,nnhc,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
      call MPI_Bcast(fbc11,nnhc,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    case(3)
  !EndIf
  !If(Ncent==3) Then
    !Call MPI_BCAST(fxbc31(1,1),natom*nnhc,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    !Call MPI_BCAST(fybc31(1,1),natom*nnhc,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    !Call MPI_BCAST(fzbc31(1,1),natom*nnhc,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_Bcast(frbc31(1,1,1),3*natom*nnhc,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    !Call MPI_BCAST(vxbc31(1,1),natom*nnhc,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    !Call MPI_BCAST(vybc31(1,1),natom*nnhc,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    !Call MPI_BCAST(vzbc31(1,1),natom*nnhc,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_Bcast(vrbc31(1,1,1),3*natom*nnhc,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    !Call MPI_BCAST(xbc31(1,1),natom*nnhc,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    !Call MPI_BCAST(ybc31(1,1),natom*nnhc,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    !Call MPI_BCAST(zbc31(1,1),natom*nnhc,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_Bcast(rbc31(1,1,1),3*natom*nnhc,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  !EndIf
  end select
#endif

  return
end subroutine

