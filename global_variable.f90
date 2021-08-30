module global_variable
  integer :: Istep
! +++ Constant parameters +++
!    real(8), parameter :: facttime  = 1.d0 / 0.024188843d0
    real(8), parameter :: facttime  = 1.d0/0.024d0
    real(8), parameter :: factmass  = 1.6605402d-27/9.1095d-31

!    real(8), parameter :: AngtoAU   = 1.d0/0.529177249d0
    real(8), parameter :: AngtoAU   = 1.d0/0.52918d0
    real(8), parameter :: bohr      = AngtoAU

!    real(8), parameter :: AUtoAng   = 0.529177249d0
    real(8), parameter :: AUtoAng   = 0.52918d0
    real(8), parameter :: bohr_inv  = AUtoAng

!    real(8), parameter :: eVtoAU    = 1.0d0/27.21162
    real(8), parameter :: KtoAU    = 1.98624d-3/627.51d0 ! Boltzmann constant K to AU (kcal/K hartree/kcal)
    real(8), parameter :: boltz     = KtoAU
    real(8), parameter :: pi = 3.14159265358979d0
    integer :: Uout
    character(len=7), parameter :: Oname='out.dat'
! +++ Constant parameters +++


! +++ Required parameters +++
  integer :: Nproc, Myrank
  integer :: Natom = 0, Nbead = 0, Nstep = 0, Nforce = 0
  integer :: simulation = -1  ! 0: PIMD, 1: RPMD, 3: CMD
  real(8) :: temperature = 0.0d0, dt = 0.0d0
  character(len=:), allocatable :: path_result
  integer :: iseed(4) = 0
! +++ End Required parameters +++

! +++ Optional parameters +++
  character(len=:), allocatable :: path_scr
  character(len=10)  :: version_gaussian
  integer :: Nref      = 5
!  integer :: iys       = 5
  integer :: Nys       = 5
  integer :: Nnhc      = 4  ! Length of Nose-Hoover chain
  integer :: order     = 4
  integer :: Nrestart  = 1
  real(8) :: gamma     = 1.0
  integer :: Nensemble = 1  ! Nensemble = 0: NVE, 1: NVT
  integer :: Ntheory   = 0  ! 0: DFT, 1:MP2
  integer :: Ncent     = 3  ! Type of thermostat = 1: atom, 3: component
  integer :: Ncolor    = 1
  real(8) :: freq1     = 10.0d0
  integer :: umbrella_atom1 = 0
  integer :: umbrella_atom2 = 0
  integer :: umbrella_atom3 = 0
  real(8) :: umbrella_height= 0.0d0

  logical :: Lrestart   = .False.
  logical :: Langstrom  = .True.
  logical :: Lrandomc   = .False.
  logical :: Lcharge    = .False.
  logical :: Lhfcc      = .False.
  logical :: Lhomolumo  = .False.
  logical :: Ldipole    = .True.
  logical :: Lpop       = .False.
  logical :: Lforce     = .False.
  logical :: Lumbrella  = .False.
  logical :: Lsavevel   = .False.
  logical :: Lsavelog   = .False.
  logical :: Lsavechk   = .False.
  logical :: Lgengau    = .False.
  logical :: Lqmmm      = .False.
  logical :: Lsaveforce = .False.

! +++ Internal parameter +++ 
  real(8), allocatable :: r(:,:,:)   ! r(3,Natom,Nbead)
  real(8), allocatable :: f(:,:,:)   ! f(3,Natom,Nbead)
  real(8), allocatable :: u(:,:,:)   ! u(3,Natom,Nbead)
  real(8), allocatable :: vu(:,:,:)   ! u(3,Natom,Nbead)
  real(8), allocatable :: fu(:,:,:), fu_ref(:,:,:)
  real(8), allocatable :: tnm(:,:), tnminv(:,:)
  real(8), allocatable :: physmass(:)
  real(8), allocatable :: energy(:)
  real(8), allocatable :: bath(:,:,:,:)   !  bath(3,Natom, Nnhc, Nbead)
  real(8), allocatable :: vbath(:,:,:,:)  ! vbath(3,Natom, Nnhc, Nbead)
  real(8), allocatable :: fbath(:,:,:,:)  ! fbath(3,Natom, Nnhc, Nbead)
  real(8), allocatable :: dipole(:,:)     ! dipole(4,Nbead)
  real(8), allocatable :: charge(:,:)     ! charge(Natom,Nbead)
  real(8) :: beta
  real(8) :: omega_system
! +++ End Internal parameter +++ 

  real(8), allocatable :: rbc11(:),      vbc11(:),      fbc11(:)
  real(8), allocatable :: rbc1(:,:),     vbc1(:,:),     fbc1(:,:)
  real(8), allocatable :: rbc31(:,:,:),  vbc31(:,:,:),  fbc31(:,:,:)
  real(8), allocatable :: rbc3(:,:,:,:), vbc3(:,:,:,:), fbc3(:,:,:,:)
  character(len=2), allocatable :: alabel(:)

  real(8) :: dt_ref !, dp, dp_inv
  real(8) :: omega_p2, omega2
  real(8) :: gnkt, gkt
  real(8), allocatable :: dnmmass(:,:), fictmass(:,:), qmass(:), ysweight(:)
  real(8), allocatable :: qmcent11(:),  qmcent1(:,:), qmcent31(:), qmcent3(:,:)
! +++ End Optional parameters +++

! +++ MPI parameters +++
  integer :: Ista, Iend
! +++ End MPI parameters +++
end module

