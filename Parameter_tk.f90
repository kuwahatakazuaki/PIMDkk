Module Parameter_tk

!  Use iso_c_binding
   Use MPI

  Double Precision, Allocatable ::  r_tmp(:,:), fr_tmp(:,:), char_tmp(:)
  Integer, Allocatable ::  atm_type(:)
  Integer   :: numeach, ista, iend, nhmod
  Integer, Allocatable :: listeach(:),listeachtmp(:),ireqa(:,:),ireqb(:,:)
  Integer, Allocatable :: recvlist(:),recvlist1(:),recvlist2(:)
!  Integer, Allocatable :: recvlistx(:),recvlistx1(:),recvlistx2(:)
!  Double Precision :: dp_inv, bohr_inv, dip_tmp(3)
  Double Precision :: dip_tmp(3)
! after fortran2003
!  Character(len=:), Allocatable  :: addresstmp, address0
  Character(len=81)  :: address0
  Character(len=87)  :: addresstmp
  Integer  ::  laddress,kproc,na31,na3,nsendrecv,ireq(9),ireq0,npacksize
  Integer  ::  nrecv, sendlist,sendlist1,sendlist2

    Integer :: mstatus(MPI_STATUS_SIZE)
!    Integer, Allocatable :: ireqe(:), ireqfx(:), ireqfy(:), ireqfz(:)
!    Integer, Allocatable :: ireqc(:), ireqdx(:), ireqdy(:), ireqdz(:)
  Double Precision, Allocatable ::  work(:), work2(:)

  TYPE gms_t
    SEQUENCE
    Integer ::  nhedder
    Integer ::  ncpus
    Integer ::  rst(3)
    Integer ::  pid
    Character, Pointer :: gmspath
    Character, Pointer :: pregms
    Character, Pointer :: setgms
    Character, Pointer :: scr
!    Character, Pointer :: scr2
    Character, Pointer :: mkscr
    Character, Pointer :: killgms
    Character, Pointer :: scrdir
    Character, Pointer :: inp
    Character, Pointer :: f10
    Character, Pointer :: f10_ln
    Character, Pointer :: job
    Character, Pointer :: xopt
    Character, Pointer :: gmsrun
    Character, Pointer :: gmstype
    Character, Pointer :: gmsrun1
    Character, Pointer :: gmsrun2
!    Character, Pointer :: hedd0
!    Character, Pointer :: heddh
!    Character, Pointer :: heddr
!    Character, Pointer :: heddE
    Character, Pointer :: hedd(:)
    Character, Pointer :: tmp
  END TYPE gms_t
  TYPE atom_t
    SEQUENCE
    Integer ::  ntype
    Double Precision, Pointer ::  nvelec
    Double Precision, Pointer ::  mass
    Double Precision, Pointer ::  smass
    Character, Pointer :: atm_name(:)
    Character, Pointer :: atm_name2(:)
    Double Precision ::  bohr
  END TYPE atom_t

  TYPE (gms_t)  :: gmsenv
  TYPE (atom_t)  :: atmparam

  End Module

