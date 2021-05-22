subroutine broad_parameters1
  use global_variable
  use mpi
  implicit none
  integer ierr

  call mpi_bcast(Natom,       1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(Nbead,       1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(Nstep,       1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(Nref,        1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(Nnhc,        1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(Nys,         1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(Ncolor,      1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(order,       1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(Nensemble,   1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(simulation,  1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(Ntheory,     1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(Nforce,      1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(Ncent,       1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(iseed,       4, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

  call mpi_bcast(temperature, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(dt,          1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(gamma,       1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(freq1,       1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

  call mpi_bcast(Lqmmm,       1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(Lgengau,     1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(Lsavevel,    1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(Lsavelog,    1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(Lsavechk,    1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(Langstrom,   1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(Lrandomc,    1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(Lcharge,     1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(Ldipole,     1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(Lhomolumo,   1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(Lpop,        1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(Lhfcc,       1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(Lsaveforce,  1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(Lumbrella,   1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(Lrestart,    1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(version_gaussian, 10, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)

!  call mpi_bcast(path_result_temp, 90, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
!  call mpi_bcast(path_scr_temp,    90, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
block
  character(len=120)  :: path_result_temp
  character(len=120)  :: path_scr_temp
  if (myrank == 0) then
    path_result_temp = path_result
    path_scr_temp    = path_scr
  end if
  call mpi_bcast(path_result_temp, 120, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(path_scr_temp,    120, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
  if (myrank /= 0) then
    path_result = trim(path_result_temp)
    path_scr    = trim(path_scr_temp)
  end if
end block

  if (Lumbrella .eqv. .True.) then
    call mpi_bcast(umbrella_atom1, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(umbrella_atom2, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(umbrella_atom3, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  end if


!! Kuwahata 2019/08/04
!  Call MPI_BCAST(umbrella_sampling,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
!  Call MPI_BCAST(umbrella_width,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
!  Call MPI_BCAST(umbrella_height,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
!! End Kuwahata 2020/0929

return
!  Call MPI_BCAST(nshoot, 1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
!  Call MPI_BCAST(ntshoot,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
!  Call MPI_BCAST(zshoot,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
!  Call MPI_BCAST(xshootmin,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
!  Call MPI_BCAST(xshootmax,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
!  Call MPI_BCAST(yshootmin,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
!  Call MPI_BCAST(yshootmax,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
!  Call MPI_BCAST(eshoot,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
!  Call MPI_BCAST(hdel,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
end subroutine broad_parameters1


subroutine broad_parameters2
  use global_variable
  use mpi
  implicit none
  integer :: i
  integer :: ierr
  real(8) :: r_temp(3,Natom)

  call mpi_bcast(physmass,  Natom, MPI_DOUBLE_PRECISION,0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(alabel,  2*Natom, MPI_CHARACTER,       0, MPI_COMM_WORLD, ierr)


  if (myrank == 0) then
    do i = 1, 3
      r_temp(i,:) = u(i,:,1)
    end do
  end if

  call mpi_bcast(r_temp,   3*Natom, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

  do i = 1, 3
    u(i,:,1) = r_temp(i,:)
  end do

!  if (myrank == 1) then
!    print *, "myrank and",myrank, "u is "
!    do i = 1, Natom
!      print *, u(:,i,1)
!    end do
!  end if

end subroutine broad_parameters2

subroutine broad_parameters3
  use global_variable
  use mpi
  implicit none
  integer :: i, j, ierr

!  do i = 0, Nproc-1
!    if (myrank == i) then
!      print *, "myrank ", myrank
!      do j = 1, Nbead
!        print *, vbath(:,:,2,j)
!      end do
!    end if
!    call mpi_barrier(MPI_COMM_WORLD,ierr)
!  end do


  call mpi_bcast(    u, 3*Natom*Nbead,      MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
  call mpi_bcast(   vu, 3*Natom*Nbead,      MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
  call mpi_bcast(   fu, 3*Natom*Nbead,      MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
  call mpi_bcast( bath, 3*Natom*Nnhc*Nbead, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
  call mpi_bcast(vbath, 3*Natom*Nnhc*Nbead, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
  call mpi_bcast(fbath, 3*Natom*Nnhc*Nbead, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )

  if ( Ncent == 1 ) then

    if ( Ncolor == 1 ) then
      call mpi_bcast( rbc11, Nnhc, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
      call mpi_bcast( vbc11, Nnhc, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
      call mpi_bcast( fbc11, Nnhc, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
    else if ( Ncolor >= 2 ) then
      call mpi_bcast( rbc1, Nnhc*Ncolor, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
      call mpi_bcast( vbc1, Nnhc*Ncolor, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
      call mpi_bcast( fbc1, Nnhc*Ncolor, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
    end if

  else if ( Ncent == 3) then

    if ( Ncolor == 1 ) then
      call mpi_bcast( rbc31, 3*Natom*Nnhc, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
      call mpi_bcast( vbc31, 3*Natom*Nnhc, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
      call mpi_bcast( fbc31, 3*Natom*Nnhc, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
    else if ( Ncolor >= 2 ) then
      call mpi_bcast( rbc3,  3*Natom*Nnhc*Ncolor, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
      call mpi_bcast( vbc3,  3*Natom*Nnhc*Ncolor, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
      call mpi_bcast( fbc3,  3*Natom*Nnhc*Ncolor, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
    end if

  end if

!  do i = 0, Nproc-1
!    if (myrank == i) then
!      print *, "myrank ", myrank
!      do j = 1, Nbead
!        print *, vbath(:,:,2,j)
!      end do
!    end if
!    call mpi_barrier(MPI_COMM_WORLD,ierr)
!  end do
end subroutine broad_parameters3


!Subroutine Broad2
!  Use Parameters
!  Use MPI
!  Implicit None
!  Double Precision,dimension(natom) :: x1,y1,z1
!
!  If(MyRank==0) Then
!    Do iatom=1,natom
!       x1(iatom)=ux(iatom,1)
!       y1(iatom)=uy(iatom,1)
!       z1(iatom)=uz(iatom,1)
!    EndDo
!  EndIf
!
!  Call MPI_BCAST(x1,natom,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
!  Call MPI_BCAST(y1,natom,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
!  Call MPI_BCAST(z1,natom,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
!
!  Do iatom=1,natom
!     ux(iatom,1)=x1(iatom)
!     uy(iatom,1)=y1(iatom)
!     uz(iatom,1)=z1(iatom)
!  EndDo
!
!  Return
!End Subroutine
!
!
!  Subroutine Broad3
!    Use Parameters
!    Use MPI
!    Implicit None
!!   Double Precision,dimension(natom,nnhc) :: x1,y1,z1,vx1,vy1,vz1
!!   Double Precision,dimension(nnhc) :: r1,v1
!    
!
!    Call MPI_BCAST(ux(1,1),natom*nbead,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
!    Call MPI_BCAST(uy(1,1),natom*nbead,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
!    Call MPI_BCAST(uz(1,1),natom*nbead,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
!    Call MPI_BCAST(vux(1,1),natom*nbead,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
!    Call MPI_BCAST(vuy(1,1),natom*nbead,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
!    Call MPI_BCAST(vuz(1,1),natom*nbead,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
!    Call MPI_BCAST(fux(1,1),natom*nbead,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
!    Call MPI_BCAST(fuy(1,1),natom*nbead,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
!    Call MPI_BCAST(fuz(1,1),natom*nbead,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
!
!    Call MPI_BCAST(vxbath(1,1,1),natom*nnhc*nbead,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
!    Call MPI_BCAST(vybath(1,1,1),natom*nnhc*nbead,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
!    Call MPI_BCAST(vzbath(1,1,1),natom*nnhc*nbead,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
!    Call MPI_BCAST(xbath(1,1,1),natom*nnhc*nbead,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
!    Call MPI_BCAST(ybath(1,1,1),natom*nnhc*nbead,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
!    Call MPI_BCAST(zbath(1,1,1),natom*nnhc*nbead,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
!    Call MPI_BCAST(fxbath(1,1,1),natom*nnhc*nbead,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
!    Call MPI_BCAST(fybath(1,1,1),natom*nnhc*nbead,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
!    Call MPI_BCAST(fzbath(1,1,1),natom*nnhc*nbead,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
!
!    If(NCent==3) Then
!     
!       If(NColor==1) Then
!
!          Call MPI_BCAST(fxbc31(1,1),natom*nnhc,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
!          Call MPI_BCAST(fybc31(1,1),natom*nnhc,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
!          Call MPI_BCAST(fzbc31(1,1),natom*nnhc,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
!          Call MPI_BCAST(vxbc31(1,1),natom*nnhc,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
!          Call MPI_BCAST(vybc31(1,1),natom*nnhc,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
!          Call MPI_BCAST(vzbc31(1,1),natom*nnhc,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
!          Call MPI_BCAST(xbc31(1,1),natom*nnhc,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
!          Call MPI_BCAST(ybc31(1,1),natom*nnhc,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
!          Call MPI_BCAST(zbc31(1,1),natom*nnhc,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
!          
!       Else
!
!          Call MPI_BCAST(fxbc3(1,1,1),natom*nnhc*ncolor,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
!          Call MPI_BCAST(fybc3(1,1,1),natom*nnhc*ncolor,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
!          Call MPI_BCAST(fzbc3(1,1,1),natom*nnhc*ncolor,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
!          Call MPI_BCAST(vxbc3(1,1,1),natom*nnhc*ncolor,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
!          Call MPI_BCAST(vybc3(1,1,1),natom*nnhc*ncolor,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
!          Call MPI_BCAST(vzbc3(1,1,1),natom*nnhc*ncolor,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
!          Call MPI_BCAST(xbc3(1,1,1),natom*nnhc*ncolor,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
!          Call MPI_BCAST(ybc3(1,1,1),natom*nnhc*ncolor,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
!          Call MPI_BCAST(zbc3(1,1,1),natom*nnhc*ncolor,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
!
!       EndIf
!
!    EndIf
!
!    If(NCent==1) Then
!
!       If(NColor==1) Then
!
!          Call MPI_BCAST(rbc11,nnhc,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
!          Call MPI_BCAST(vbc11,nnhc,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
!          Call MPI_BCAST(fbc11,nnhc,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
!
!       Else
!
!          Call MPI_BCAST(rbc1(1,1),nnhc*ncolor,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
!          Call MPI_BCAST(vbc1(1,1),nnhc*ncolor,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
!          Call MPI_BCAST(fbc1(1,1),nnhc*ncolor,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
!
!       EndIf
!    
!    EndIf
!
!    Return
!  End Subroutine
!!
