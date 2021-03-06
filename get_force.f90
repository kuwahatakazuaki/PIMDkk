subroutine get_force
  use global_variable
  use print_out
  use utility

  select case(Nforce)
!    case(1)
!      call force_mopac
!    case(2)
!      call force_siesta
    case(6)
      call force_gaussian
    case(11)
      call force_morse
    case default
      call program_abort('ERROR!!! Wrong "NForce" option')
  end select
contains

  subroutine force_morse
    use mpi
    use communication
    implicit none

  ! Parameters in Atomic unit
  ! Hydrogen molecule
    real(8), parameter :: r0     = 1.41014d0
    real(8), parameter :: De     = 0.1745d0
    real(8), parameter :: curvat = 1.0213d0

    integer :: Imode, i, j
    real(8) :: rij(3), fij(3), power, dis

    call comm_coord

    f(:,:,:) = 0.0d0
    do Imode = Ista, Iend
      do i = 1, Natom
        do j = i+1, Natom
          rij(:) = r(:,i,Imode) - r(:,j,Imode)
          dis = dsqrt( dot_product(rij(:), rij(:)) )
          power = (-1) * curvat * (dis - r0)
          fij(:) = 2 * curvat * De * exp(power) * (exp(power)-1) * rij(:) / dis
          f(:,i,Imode) = f(:,i,Imode) + fij(:)
          f(:,j,Imode) = f(:,j,Imode) - fij(:)
        end do
      end do
      f(:,:,Imode) = f(:,:,Imode) / dble(Nbead)
    end do

!if ( Myrank == 0 ) then
!  do j = 1, Nbead
!    do i = 1, Natom
!      print *, f(:,i,j)
!    end do
!  end do
!end if

    call comm_output
    if ( Myrank == 0 ) call print_result
  end subroutine force_morse




! ++++++++++++++++++++++
! +++ force_gaussian +++
! ++++++++++++++++++++++
  subroutine force_gaussian
    use mpi
    use communication
    implicit none
    integer :: i,Imode, Iposi
    integer :: Utemp, Ierr
    integer :: dummyI
    real(8) :: dummyD
    character :: dummyC
    real(8) :: dp_inv
    character(len=120) :: line
    character(len=:), allocatable  :: key1, key2, key3, key4, key5
    character(len=:), allocatable :: path_scr_sub
    character(len=5)   :: Cnum

    select case (Ntheory)
      case(0)
        key1  = 'SCF Done'
      case(1)
        key1  = 'EUMP2'
      case default
        stop 'ERROR!!! Wrong "theory" option'
    end select

    key2  = 'Number     Number              X              Y              Z'
    key3  = 'Dipole moment'
    key4  = 'Mulliken charges'
    key5  = 'Isotropic Fermi Contact Couplings'

    call comm_coord

    do Imode = Ista, Iend
      write(Cnum,'(I5.5)') Imode
      path_scr_sub = path_scr//'/'//Cnum//'/'

      call system('cat '//path_scr_sub//'gauss.tmp1 > '//path_scr_sub//'gauss.com')

      open(newunit=Utemp,file=path_scr_sub//'gauss.com',status='old',position='append')
        do i = 1, Natom
          write(Utemp,*) alabel(i), r(:,i,Imode) * AUtoAng
        end do
        write(Utemp,*)
      close(Utemp)

  !    call system('cat '//path_scr_sub//'gauss.xyz >> '//path_scr_sub//'gauss.com')

      call system(path_scr//'/'//'g0xrun_p '//path_scr_sub//'gauss.com '//path_scr_sub//'gauss.log '//path_scr_sub)

      open(newunit=Utemp,file=path_scr_sub//'gauss.log',status='old',iostat=Ierr)
        if ( Ierr /= 0 ) then
          call program_abort("ERROR!! opne the gaussian output")
        end if

  !  +++ Reading "SCF Done" +++
        do
          read(Utemp,'(a)',end=401) line
          if ( index(line,key1) > 0 ) exit
        end do
        if ( Ntheory == 0 ) then
          Iposi = index(line,'=')
          read(line(Iposi+2:Iposi+17),*) energy(Imode)
        else if ( Ntheory == 1 ) then
          read(line(38:60),*) energy(Imode)
        end if
  !  +++ End Reading "SCF Done" +++

  !  +++ Reading "Mulliken charge" +++
        if ( Lcharge .eqv. .True. ) then
          do
            read(Utemp,'(a)', end=404) line
            if ( index(line,key4) > 0) exit
          end do
          read(Utemp,'()')
          do i = 1, Natom
            read(Utemp,*) dummyI, dummyC, charge(i,Imode)
          end do
        end if
  !  +++ End Reading "Mulliken charge" +++

  !  +++ Reading "Dipole moment" +++
        if ( Ldipole .eqv. .True. ) then
          do
            read(Utemp,'(a)', end=403) line
            if ( index(line,key3) > 0) exit
          end do
          read(Utemp,*) dummyC, dipole(1,Imode), dummyC, dipole(2,Imode), dummyC, dipole(3,Imode), dummyC, dipole(4,Imode)
        end if
  !  +++ End Reading "Dipole moment" +++

  !  +++ Reading "Atomic Force" +++
        do
          read(Utemp,'(a)') line
          if ( index(line,key2) > 0 ) exit
        end do
        read(Utemp,'()')
        do i = 1, Natom
          read(Utemp,*) dummyI, dummyI, f(:,i,Imode)
        end do
  !  +++ End Reading "Atomic Force" +++

      close(Utemp)
      call system('rm -rf '//path_scr_sub//'Gau*')
      dp_inv = 1 / dble(Nbead)
      f(:,:,Imode) = f(:,:,Imode) * dp_inv
    end do

    call comm_output
  
  !print *, "After"
  !do j = 1, Nproc
  !  if ( Myrank == j -1 ) then
  !    print *, Myrank
  !    do i = 1, Natom
  !      print *, "alabel", alabel(i)
  !    end do
  !  end if
  !call mpi_barrier(MPI_COMM_WORLD,Ierr)
  !end do
  
  !if ( Myrank == 0 ) then
  !  print *, Myrank
  !  do j = 1, Nbead
  !!    do i = 1, Natom
  !!      print *, f(:,i,j)
  !!    end do
  !    print *, energy(j)
  !    print *, charge(:,j)
  !    print *, dipole(:,j)
  !  end do
  !end if
  !call mpi_barrier(MPI_COMM_WORLD, Ierr)
  
   if  ( Myrank == 0 ) call print_result
  
    return
    401 call program_abort('ERROR!!: We can not find "SCF Done" or "EUMP2" in Gaussian output')
    402 call program_abort('ERROR!!: We can not find "Force" in Gaussian output')
    403 call program_abort('ERROR!!: We can not find "Dipole moment" in Gaussian output')
    404 call program_abort('ERROR!!: We can not find "Mulliken charges" in Gaussian output')
    405 call program_abort('ERROR!!: We can not find "Fermi Contact" in Gaussian output')
  end subroutine force_gaussian
! ++++++++++++++++++++++
! +++ force_gaussian +++
! ++++++++++++++++++++++


end subroutine get_force

