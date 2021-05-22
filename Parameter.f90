Module Parameters

    !! General Parameters
!   Double Precision, Parameter  :: pi        = 3.14159265358979d0
!   ! . . > fs -- > a.u.
!   Double Precision, Parameter  :: facttime  = 1.d0 / 0.024188843d0
!   !! Mass amu -- > A.U.
!   Double Precision, Parameter  :: factmass  = 1.6605402d-27 / 9.1093897d-31
!   ! . . > Boltzmann constant
!   Double Precision, Parameter  :: boltz     = 0.316682968d-5
!   ! . . > Boltzmann constant
!   Double Precision, Parameter  :: bohr      = 0.529177249D+00

    ! . . > fs -- > a.u.
    Double Precision, Parameter  :: facttime  = 1.d0/0.024d0
!   Double Precision, Parameter  :: facttime  = 1.d0 / 0.024188843d0
    !! Mass amu -- > A.U.
    Double Precision, Parameter  :: factmass  = 1.6605402d-27/9.1095d-31
    ! . . > Boltzmann constant
    Double Precision, Parameter  :: boltz     = 1.98624d-3/627.51d0 ! Boltzmann constant K to AU (kcal/K hartree/kcal)
!    Double Precision, Parameter  :: boltz     = 8.617333262145d-5/27.211396132 ! Boltzmann constant K to AU (eV/K hartree/eV) from NIST
    ! . . > Boltzmann constant
    Double Precision, Parameter  :: bohr      = 1.d0/0.52918d0
!    Double Precision, Parameter  :: bohr      = 1.d0/0.529177249d0
  !! +++ Constants for conversion +++
  !real(8), parameter :: eVtoAU    = 1.0d0/27.21162
  !real(8), parameter :: AngtoAU = 1/0.529177249d0
  !real(8), parameter :: AUtoAng = 0.529177249d0

!  Kuwahata 2019/12/08
!    Double Precision, Parameter  :: bohr_inv      = 0.529177249d0
    Double Precision, Parameter  :: bohr_inv      = 0.52918d0  ! AUtoAng
!  End Kuwahata 2019/12/08

    Double Precision, Allocatable ::  x(:,:),  y(:,:),  z(:,:)
    Double Precision, Allocatable :: fx(:,:), fy(:,:), fz(:,:)
    Double Precision, Allocatable :: ux(:,:), uy(:,:), uz(:,:)  ! coordinate
    Double Precision, Allocatable :: vux(:,:),vuy(:,:),vuz(:,:)
    Double Precision, Allocatable :: fux(:,:),fuy(:,:),fuz(:,:)
    Double Precision, Allocatable :: fux_ref(:,:),fuy_ref(:,:),fuz_ref(:,:)
    Double Precision, Allocatable :: tnm(:,:), tnminv(:,:)
    Double Precision, Allocatable :: u(:,:), uinv(:,:)
    Double Precision, Allocatable :: pot(:), physmass(:)
    Double Precision, Allocatable :: xbath(:,:,:),ybath(:,:,:),zbath(:,:,:)
    Double Precision, Allocatable :: vxbath(:,:,:),vybath(:,:,:),vzbath(:,:,:)
    Double Precision, Allocatable :: fxbath(:,:,:),fybath(:,:,:),fzbath(:,:,:)
    Double Precision, Allocatable :: rbc11(:),vbc11(:),fbc11(:)
    Double Precision, Allocatable :: rbc1(:,:),vbc1(:,:),fbc1(:,:)
    Double Precision, Allocatable :: xbc31(:,:),vxbc31(:,:),fxbc31(:,:)
    Double Precision, Allocatable :: ybc31(:,:), vybc31(:,:), fybc31(:,:)
    Double Precision, Allocatable :: zbc31(:,:), vzbc31(:,:), fzbc31(:,:)
    Double Precision, Allocatable :: xbc3(:,:,:), vxbc3(:,:,:), fxbc3(:,:,:)
    Double Precision, Allocatable :: ybc3(:,:,:), vybc3(:,:,:), fybc3(:,:,:)
    Double Precision, Allocatable :: zbc3(:,:,:), vzbc3(:,:,:), fzbc3(:,:,:)
    Double Precision, Allocatable :: dnmmass(:,:),fictmass(:,:),qmass(:),ysweight(:)
    Double Precision, Allocatable :: qmcent11(:)
    Double Precision, Allocatable :: qmcent1(:,:)
    Double Precision, Allocatable :: qmcent31(:)
    Double Precision, Allocatable :: qmcent3(:,:)
    Double Precision, Allocatable :: dipole(:),dipolex(:),dipoley(:),dipolez(:)
    Double Precision, Allocatable :: charge(:,:),nbo(:,:),Eenergy(:),homo(:),lumo(:)
    Double Precision, Allocatable :: hfcc(:,:)
    Double Precision, Allocatable :: hvelin(:),hvelout(:),hinitx(:),hinity(:)
    character (Len=2),Allocatable :: alabel(:)
    integer         , Allocatable :: no_atom(:)
    integer         , Allocatable :: nhspin(:)
    Double Precision, Allocatable :: hess(:,:,:), pot_ti(:)
! Kuwahata Only used in Harmonic potential ??
    Double Precision, Allocatable :: fx_ti(:,:), fy_ti(:,:), fz_ti(:,:)      !
    Double Precision, Allocatable :: fx_org(:,:), fy_org(:,:), fz_org(:,:)   !
! Kuwahata End Only used in Harmonic potential ??

!YK
!   Double Precision, External         :: gasdev
    Double Precision              :: gamma, gamma2, omega_system
    Double Precision              :: omega_p, omega_p2, omega2, vsigma, usigma
    Double Precision              :: gkt, gnkt, dp, dp_inv
    Double Precision              :: factor_color,pi
!YK Made Some Variables Global for Ham_Temp
    Double Precision              :: ebath,ebath_cent,dkinetic,qkinetic
    Integer                       :: nref, nys, nnhc
    Integer                       :: imode, iatom
    Integer                       :: angstrom
!YK
    Integer                       :: icolor,inhc
!YK
    Integer                       :: Order
!YK Type of Ensemble NVE=0, NVT with Nose-Hoover Chain=1
    Integer                       :: Nensemble
!YK Set the method for electronic structure calculation
    Integer                       :: theory
!YK Switch for QM/MM calculation
    Integer                       :: NQMMM
!YK If additional basis sets for g03 or g09 calculations are necessary
    Integer                       :: NGenGau
!YK Frequency of saving gaussian inputs, outputs, and chk files
    Integer                       :: NSavevel, NSavelog, NSavechk
!YK Method and Color of NHC on centroid
    Integer                       :: NCent, NColor
!YK Switch for Restarting from previous MD, Number of the previous MD step
    Integer                       :: NRestart,Nrstep
!YK Added for MPI
    Integer                       :: NProcs,MyRank,IERR
!YK Added to Specify How to Calculate Force
    Integer                       :: NForce
!YK Added to Set Random Number Generator Seed
    Integer                       :: ISEED1,ISEED2,ISEED3,ISEED4
!YK
    Double Precision              :: beta, temperature, dt, dt_ref
    Double Precision              :: temp_step
    Double Precision              :: potential, hamiltonian, temp, potential_old

    Integer                       :: natom, nbead, nstep, Simulation
!YK IO parameters
    Integer, Parameter            :: irst=99  ! restart
    Integer, Parameter            :: igetf=98
    Integer, Parameter            :: igetx=97
    Integer, Parameter            :: igetd=96
    Integer, Parameter            :: inbo=95
    Integer, Parameter            :: iham=94
    Integer, Parameter            :: imopac=93
    Integer, Parameter            :: igauss=92
    Integer, Parameter            :: igete=91
    Integer, Parameter            :: igetxyz=90
    Integer, Parameter            :: igetc=89
    Integer, Parameter            :: igethl=88
    Integer, Parameter            :: igetv=85
    Integer, Parameter            :: igetvb=84

!YK Read Labels from the Input
    character (Len=80)            :: address
    character (Len=80)            :: address2
!YK
    Double Precision :: hdel, zshoot, xshootmin, xshootmax
    Double Precision :: yshootmin,yshootmax,eshoot,vhydinit
    Double Precision :: dkinetic_noncent,qkinetic_noncent,temp_noncent,total_noncent
    Double Precision :: dkinetic_cent,qkinetic_cent,temp_cent,total_cent
    Double Precision :: E_Virial,freq1
    Integer, Allocatable :: index_hyd(:)
    Integer, Allocatable :: ncontact(:)

    Integer          :: nshoot,nhydtot,nhyddel,nmshot,nhydads
    Integer          :: nhydrfl,nhydthr,ntshoot,natom0,nrandomc
    Integer          :: istepsv,nocharge,nodipole,nohomo,npop,nohfcc
! Kuwahata 2019/08/04 add some key parameters
    character (Len=80) :: version_gaussian
    logical :: Save_force = .False.
    logical :: umbrella_sampling = .False.
    real(8) :: umbrella_width, umbrella_height ! in K ?
    integer :: umbrella_atom1, umbrella_atom2
! End Kuwahata 2020/10/06

End Module
