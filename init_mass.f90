subroutine init_mass
  use global_variable
  use mpi
  implicit none
  integer :: i, imode
  real(8) :: dp, di
  integer :: ierr
  dp = dble(Nbead)

  dnmmass(:,1) = 0.0d0

  do i = 1, Natom
    dnmmass(i,Nbead) = 4.0*dp*physmass(i)
    do imode = 1, (Nbead-2)/2
      di = dble(imode)
      dnmmass(i,2*imode) = 2.0 * (1.0 - dcos(2*pi*di/dp)) * dp * physmass(i)
      dnmmass(i,2*imode+1) = dnmmass(i,2*imode)
    end do
  end do

! simulation =>  0: PIMD, 1: RPMD, 3: CMD
  if ( simulation == 1 ) then ! RPMD
    do i = 1, Natom
      fictmass(i,:) = physmass(i)
    end do
  else
    do i = 1, Natom
      fictmass(i,1) = physmass(i)
      do imode = 2, Nbead
!        fictmass(i,imode) = dnmmass(i,imode)
        fictmass(i,imode) = gamma*gamma * dnmmass(i,imode)
      end do
    end do
  end if

!  print *, "myrank ", myrank
!  if ( myrank == 0 ) then
!    do i = 1, Natom
!      print *, dnmmass(i,:)
!    end do
!    print *, ""
!    do i = 1, Natom 
!      print *, fictmass(i,:)
!    end do
!  end if
!  call mpi_barrier(MPI_COMM_WORLD,ierr)
!  stop

end subroutine init_mass


