module Parameters
  Double Precision, Parameter  :: pi        = 3.14159265358979d0
  ! . . > fs -- > a.u.
!  Double Precision, Parameter  :: facttime  = 1.d0/0.024d0
  Double Precision, Parameter  :: facttime  = 1.d0 / 0.024188843d0   ! fstoAU
  !! Mass amu -- > A.U.
  Double Precision, Parameter  :: factmass  = 1.6605402d-27/9.1093897d-31
  !Double Precision, Parameter  :: factmass  = 1.6605402d-27/9.1095d-31  ! 1822.866458
  ! . . > Boltzmann constant
!  Double Precision, Parameter  :: boltz     = 1.98624d-3/627.51d0 ! Boltzmann constant K to AU (kcal/K hartree/kcal)
!  Double Precision, Parameter  :: boltz     = 8.617333262145d-5/27.211396132 ! Boltzmann constant K to AU (eV/K hartree/eV) from NIST
  ! . . > Boltzmann constant
  !Double Precision, Parameter  :: bohr      = 1.d0/0.529177249d0  ! AngtoAU
  !Double Precision, Parameter  :: bohr_inv      = 0.529177249d0   ! AUtoAng

! +++ Constants for conversion +++
  real(8), parameter :: eVtoAU  = 1.0d0/27.21162d0
  real(8), parameter :: AngtoAU = 1.0d0/0.529177249d0
  real(8), parameter :: AUtoAng = 0.529177249d0   ! bohr_inv
  real(8), parameter :: AUtoJ   = 4.35974434d-18
  real(8), parameter :: KtoAU   = 8.617333262145d-5/27.211396132d0 ! Boltzmann constant K to AU (eV/K hartree/eV) from NIST
  real(8), parameter :: eVAng2AU= eVtoAU * AUtoAng
  integer, parameter :: LINELEN = 120


  integer              :: Nref, Nys, Nnhc
  real(8), allocatable :: r(:,:,:), fr(:,:,:), ur(:,:,:), vur(:,:,:)
  real(8), allocatable :: fur(:,:,:), fur_ref(:,:,:)
  real(8), allocatable :: rbath(:,:,:,:), vrbath(:,:,:,:), frbath(:,:,:,:)
  real(8), Allocatable :: rbc11(:),vbc11(:),fbc11(:)
  real(8), allocatable :: rbc31(:,:,:), vrbc31(:,:,:), frbc31(:,:,:)
  real(8), allocatable :: tnm(:,:), tnminv(:,:)
  real(8), allocatable :: u(:,:), uinv(:,:)
  real(8), allocatable :: pot(:), physmass(:)
  real(8), allocatable :: dnmmass(:,:),fictmass(:,:),qmass(:),ysweight(:)
  real(8), allocatable :: qmcent11(:)
  real(8), allocatable :: qmcent31(:)
  real(8), allocatable :: dipoler(:,:), atom_num(:)
  real(8), Allocatable :: charge(:,:),nbo(:,:),Eenergy(:),homo(:),lumo(:), hfcc(:,:)
  real(8)              :: gamma, gamma2, omega_system
  real(8)              :: omega_p2, omega2
  real(8)              :: gkt, gnkt, dp_inv
  real(8)              :: ebath, ebath_cent, dkinetic, qkinetic
  character(Len=2), allocatable :: alabel(:)
  character(len=8)     :: name_simulation
! !YK Type of Ensemble NVE=0, NVT with Nose-Hoover Chain=1
!YK Method and Color of NHC on centroid
  Integer                       :: Ncent
!YK Switch for Restarting from previous MD, Number of the previous MD step
  integer                       :: Nrstep = 0
  Integer                       :: Nproc, MyRank
!YK Added to Specify How to Calculate Force
  Integer                       :: Iforce
!YK Added to Set Random Number Generator Seed
  integer :: Iseeds(4)
  Double Precision              :: beta, temperature, dt, dt_ref
  Double Precision              :: potential, hamiltonian, temp

  Integer                       :: Natom, Nbead, Nstep, Isimulation
  character(len=9)   :: Finp = "input.inp"
  character(len=7)   :: Fout = "std.out"
  character(Len=80)  :: address
  character(Len=80)  :: address2
  Character(len=81)  :: address0
  Character(len=87)  :: addresstmp
  Double Precision :: E_Virial
  real(8) :: freq1 = 10.0d0

  integer :: istepsv !, nrandomc
  logical :: Lsave_force  = .False.
  logical :: Lsave_npa    = .False.
  logical :: Lsave_charge = .False.
  logical :: Lsave_dipole = .False.
  logical :: Lsave_hfcc   = .False.
  logical :: Lsave_energy = .False.
  logical :: Langstrom    = .True.
  logical :: Lperiodic    = .False.
  logical :: Lrestart     = .False.
  logical :: Lrandom_coor = .False.

! Kuwahata 2019/08/04 add some key parameters
  character (Len=80) :: version = "g16"
  real(8), allocatable :: pressure(:)
  real(8) :: virial, PV
  integer :: Iumb = 0
  integer :: umb_atom1 = 0, umb_atom2 = 0, umb_atom3 = 0
  real(8) :: umb_cons = 1d-5, umb_pot
! End Kuwahata 2020/10/06

  integer   :: Ista, Iend, laddress
  integer, Allocatable :: listeach(:),listeachtmp(:),ireqa(:,:),ireqb(:,:)
end module Parameters

!  integer  ::  laddress,na31,na3,nsendrecv,ireq(9),ireq0,npacksize
!  integer  ::  nrecv, sendlist,sendlist1,sendlist2
!  integer, Allocatable :: recvlist(:),recvlist1(:),recvlist2(:)
!!YK Set the method for electronic structure calculation
!    Integer                       :: theory
!!YK If additional basis sets for g03 or g09 calculations are necessary
!    Integer                       :: NGenGau
!!YK Switch for QM/MM calculation
!    Integer                       :: NQMMM

