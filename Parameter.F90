Module Parameters
   Double Precision, Parameter  :: pi        = 3.14159265358979d0
!   !! Mass amu -- > A.U.
!   Double Precision, Parameter  :: factmass  = 1.6605402d-27 / 9.1093897d-31
!   ! . . > Boltzmann constant
!   Double Precision, Parameter  :: boltz     = 0.316682968d-5

    ! . . > fs -- > a.u.
!    Double Precision, Parameter  :: facttime  = 1.d0/0.024d0
   Double Precision, Parameter  :: facttime  = 1.d0 / 0.024188843d0   ! fstoAU
    !! Mass amu -- > A.U.
    Double Precision, Parameter  :: factmass  = 1.6605402d-27/9.1095d-31  ! 1822.866458
    ! . . > Boltzmann constant
!    Double Precision, Parameter  :: boltz     = 1.98624d-3/627.51d0 ! Boltzmann constant K to AU (kcal/K hartree/kcal)
    Double Precision, Parameter  :: boltz     = 8.617333262145d-5/27.211396132 ! Boltzmann constant K to AU (eV/K hartree/eV) from NIST
    ! . . > Boltzmann constant
!    Double Precision, Parameter  :: bohr      = 1.d0/0.52918d0
    Double Precision, Parameter  :: bohr      = 1.d0/0.529177249d0
  ! +++ Constants for conversion +++
    real(8), parameter :: eVtoAU    = 1.0d0/27.21162
    real(8), parameter :: AngtoAU = 1/0.529177249d0
    real(8), parameter :: AUtoAng = 0.529177249d0
    real(8), parameter :: AUtoJ   = 4.35974434d-18

    Double Precision, Parameter  :: bohr_inv      = 0.529177249d0
!    Double Precision, Parameter  :: bohr_inv      = 0.52918d0  ! AUtoAng

    real(8), allocatable :: r(:,:,:), fr(:,:,:), ur(:,:,:), vur(:,:,:)
    real(8), allocatable :: fur(:,:,:), fur_ref(:,:,:)
    real(8), allocatable :: rbath(:,:,:,:), vrbath(:,:,:,:), frbath(:,:,:,:)
    Double Precision, Allocatable :: rbc11(:),vbc11(:),fbc11(:)
    !Double Precision, Allocatable :: rbc1(:,:),vbc1(:,:),fbc1(:,:)
    real(8), allocatable :: rbc31(:,:,:), vrbc31(:,:,:), frbc31(:,:,:)
    !real(8), allocatable :: rbc3(:,:,:,:), vrbc3(:,:,:,:), frbc3(:,:,:,:)
    Double Precision, Allocatable :: tnm(:,:), tnminv(:,:)
    Double Precision, Allocatable :: u(:,:), uinv(:,:)
    Double Precision, Allocatable :: pot(:), physmass(:)
    Double Precision, Allocatable :: dnmmass(:,:),fictmass(:,:),qmass(:),ysweight(:)
    Double Precision, Allocatable :: qmcent11(:)
    !Double Precision, Allocatable :: qmcent1(:,:)
    Double Precision, Allocatable :: qmcent31(:)
    !Double Precision, Allocatable :: qmcent3(:,:)
    real(8), allocatable :: dipoler(:,:), atom_num(:)
    Double Precision, Allocatable :: charge(:,:),nbo(:,:),Eenergy(:),homo(:),lumo(:), hfcc(:,:)
    character (Len=2),Allocatable :: alabel(:)
    character(len=8) :: simulation

    Double Precision              :: gamma, gamma2, omega_system
    Double Precision              :: omega_p, omega_p2, omega2, vsigma, usigma
    Double Precision              :: gkt, gnkt, dp, dp_inv
    Double Precision              :: factor_color
!YK Made Some Variables Global for Ham_Temp
    Double Precision              :: ebath,ebath_cent,dkinetic,qkinetic
    Integer                       :: Nref, Nys, Nnhc
! !YK Type of Ensemble NVE=0, NVT with Nose-Hoover Chain=1
!     Integer                       :: Iensemble
!YK Method and Color of NHC on centroid
    Integer                       :: NCent
!YK Switch for Restarting from previous MD, Number of the previous MD step
    integer                       :: Nrstep = 0
!YK Added for MPI
    Integer                       :: NProcs,MyRank!,IERR
!YK Added to Specify How to Calculate Force
    Integer                       :: Iforce
!YK Added to Set Random Number Generator Seed
    integer :: Iseeds(4)
    Double Precision              :: beta, temperature, dt, dt_ref
    Double Precision              :: potential, hamiltonian, temp!, potential_old

    Integer                       :: Natom, Nbead, Nstep, Isimulation
    character(len=9) :: Finp = "input.inp"
    character(len=7) :: Fout = "std.out"

!YK Read Labels from the Input
    character (Len=80)            :: address
    character (Len=80)            :: address2
    Character(len=81)  :: address0
    Character(len=87)  :: addresstmp
    Double Precision :: E_Virial! ,freq1
    real(8) :: freq1 = 10.0d0

    integer :: istepsv !, nrandomc
    logical :: Lsave_force  = .False.
    logical :: Lsave_npa    = .False.
    logical :: Lsave_charge = .False.
    logical :: Lsave_dipole = .False.
    logical :: Lsave_hfcc   = .False.
    logical :: Langstrom    = .True.
    logical :: Lrestart     = .False.
    logical :: Lrandom_coor = .False.

! Kuwahata 2019/08/04 add some key parameters
    character (Len=80) :: version = "g16"
    real(8), allocatable :: pressure(:)
    real(8) :: virial, PV
    integer :: umbrella_sampling = 0
    real(8) :: umbrella_width, umbrella_height, umbrella_constant
    integer :: umbrella_atom1 = 0, umbrella_atom2 = 0, umbrella_atom3 = 0
! End Kuwahata 2020/10/06

  Integer   :: numeach, ista, iend, nhmod
  Integer, Allocatable :: listeach(:),listeachtmp(:),ireqa(:,:),ireqb(:,:)
  Integer, Allocatable :: recvlist(:),recvlist1(:),recvlist2(:)
  Integer  ::  laddress,na31,na3,nsendrecv,ireq(9),ireq0,npacksize
  Integer  ::  nrecv, sendlist,sendlist1,sendlist2
  Double Precision, Allocatable ::  work(:), work2(:)
End Module

!!YK Set the method for electronic structure calculation
!    Integer                       :: theory
!!YK If additional basis sets for g03 or g09 calculations are necessary
!    Integer                       :: NGenGau
!!YK
!    Integer                       :: Order
!!YK Switch for QM/MM calculation
!    Integer                       :: NQMMM
