module Parameters
  real(8), parameter :: pi        = 3.14159265358979d0

! +++ Constants for conversion +++
  ! fs -- > a.u.
  real(8), parameter :: fs2AU  = 1.d0 / 0.024188843d0   ! fs2AU
  ! Mass amu -- > A.U.
  real(8), parameter :: factmass  = 1.6605402d-27/9.1093897d-31 ! mass2AU
  real(8), parameter :: eVtoAU    = 1.0d0/27.21162d0
  real(8), parameter :: AngtoAU   = 1.0d0/0.529177249d0
  real(8), parameter :: AUtoAng   = 0.529177249d0   ! bohr_inv
  real(8), parameter :: AUtoJ     = 4.35974434d-18
  real(8), parameter :: KtoAU     = 8.617333262145d-5/27.211396132d0 ! Boltzmann constant K to AU (eV/K hartree/eV) from NIST
  real(8), parameter :: eVAng2AU  = eVtoAU * AUtoAng
  real(8), parameter :: AU2kcalPmol = 627.509469d0 ! 1 Hartree = 627.509469 kcal/mol
  real(8), parameter :: kcalPmol2AU = 1.0d0/AU2kcalPmol
  integer, parameter :: LINELEN = 120

  integer              :: Natom, Nbead, Nstep, Isimulation
  integer              :: Nref, Nys, Nnhc, out_step = 1
  real(8), allocatable :: r(:,:,:), fr(:,:,:), ur(:,:,:), vur(:,:,:)
  real(8), allocatable :: fur(:,:,:), fur_ref(:,:,:)
  real(8), allocatable :: rbath(:,:,:,:), vrbath(:,:,:,:), frbath(:,:,:,:)
  real(8), Allocatable :: rbc11(:),vbc11(:),fbc11(:)
  real(8), allocatable :: rbc31(:,:,:), vrbc31(:,:,:), frbc31(:,:,:)
  real(8), allocatable :: tnm(:,:), tnminv(:,:)
  real(8), allocatable :: u(:,:), uinv(:,:)
  real(8), allocatable :: pot(:), physmass(:)
  real(8), allocatable :: dnmmass(:,:),fictmass(:,:),qmass(:),ysweight(:)
  real(8), allocatable :: qmcent11(:), qmcent31(:)
  real(8), allocatable :: dipoler(:,:), atom_num(:)
  real(8), allocatable :: charge(:,:),nbo(:,:),Eenergy(:),homo(:),lumo(:), hfcc(:,:)
  real(8)              :: gamma1 = 1.0d0, gamma2, omega_system
  real(8)              :: omega_p2, omega2
  real(8)              :: gkt, gnkt, dp_inv, E_Virial
  real(8)              :: ebath, ebath_cent, dkinetic, qkinetic
  real(8)              :: beta, temperature, dt, dt_ref
  real(8)              :: potential, hamiltonian, temp
  real(8) :: freq1 = 10.0d0

  integer              :: Ncent  ! Type of the thermostat
  integer              :: Irestep = 0
  integer              :: Nproc, MyRank
  integer              :: Iforce ! Type of force
  integer              :: Iseeds(4) ! Random Number Generator Seed
  integer              :: Ista, Iend, laddress
  integer, allocatable :: listeach(:),listeachtmp(:) !,ireqa(:,:),ireqb(:,:)

  character(Len=2), allocatable :: alabel(:)
  character(len=8)     :: name_simulation
  character(len=9)     :: Finp = "input.inp"
  character(len=7)     :: Fout = "std.out"
  character(len=7)     :: Ferr = "std.err"
  character(Len=80)    :: dir_result
  character(Len=80)    :: dir_scr
  Character(len=81)    :: address0
  Character(len=87)    :: addresstmp

  integer :: istepsv = 0
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
  character (Len=20) :: version = "g16"
  real(8), allocatable :: pressure(:)
  real(8) :: virial, PV
  integer :: Iumb = 0
  integer :: umb_atom1 = 0, umb_atom2 = 0, umb_atom3 = 0
  real(8) :: umb_cons = 1d-5, umb_pot
! End Kuwahata 2020/10/06

  integer :: Icons = 0, cons_atom1, cons_atom2, cons_atom3
  real(8) :: cons_strenght, cons_val
end module Parameters

module Struct_
  implicit none
  TYPE xyz_
     REAL(8) :: x, y, z
  END TYPE xyz_

  TYPE abc_
     REAL(8) :: a, b, c
  END TYPE abc_

  TYPE force_
     REAL(8) :: x, y, z
     REAL(8) :: norm
  END TYPE force_
end module Struct_

!!YK Set the method for electronic structure calculation
!    Integer                       :: theory
!!YK If additional basis sets for g03 or g09 calculations are necessary
!    Integer                       :: NGenGau
!!YK Switch for QM/MM calculation
!    Integer                       :: NQMMM

