! Last modified; 2025/01/13 by M.A.

  MODULE Precision_
    IMPLICIT NONE
    INTEGER, PARAMETER :: sp = SELECTED_REAL_KIND( 6)   ! Single precision
    INTEGER, PARAMETER :: dp = kind(1d0)                ! Double precision
    INTEGER, PARAMETER :: qp = SELECTED_REAL_KIND(33)   ! Quad   precision
    INTEGER, PARAMETER :: wp = dp
  END MODULE Precision_

  MODULE Constants_
    USE Precision_
    IMPLICIT NONE

  ! General parameters
    INTEGER,     PARAMETER :: MaxDiff    = 100                                    ! Parameter for differentiation
    REAL(wp),    PARAMETER :: pi         = 3.1415926535897932384626433832795_wp   ! Circular constant
    REAL(wp),    PARAMETER :: EulerGamma = 0.577215664901532860606512090082_wp    ! Euler-Mascheroni constant
    COMPLEX(wp), PARAMETER :: iunit      = (0.0_wp,1.0_wp)                        ! Imaginary unit

  ! Physical constants
    REAL(wp), PARAMETER :: kB = 3.166764d-6       ! Boltzmann constant (a.u./K)
    REAL(wp), PARAMETER :: mu = 1.66053886d-27    ! Atomic mass unit [kg]
    REAL(wp), PARAMETER :: qe = 1.602176462d-19   ! Elementary charge [C] => e [J]
    REAL(wp), PARAMETER :: x0 = 1.0d-10           ! Scaling length = 1 [A]

  ! Unit conversion
    REAL(wp), PARAMETER :: au2Ang      = 0.52917721092_wp    ! 1 Bohr = 0.52917721092  Angstrom (Length)
    REAL(wp), PARAMETER :: au2nm       = 0.052917721092_wp   ! 1 Bohr = 0.052917721092 Angstrom (Length)
    REAL(wp), PARAMETER :: Ry2eV       = 13.60569253_wp      ! 1 Ry      = 13.60569253 eV (Energy)
    REAL(wp), PARAMETER :: au2eV       = 27.21138386_wp      ! 1 Hartree = 27.21138386 eV (Energy)
    REAL(wp), PARAMETER :: au2kcalPmol = 627.509469_wp       ! 1 Hartree = 627.509469 kcal/mol (Energy)
    REAL(wp), PARAMETER :: au2kJPmol   = 2625.49962_wp       ! 1 Hartree = 2625.49962 kJ/mol (Energy)
    REAL(wp), PARAMETER :: au2cmi      = 219474.6313705_wp   ! 1 Hartree = 219474.6313705 cm^-1 (Energy)
    REAL(wp), PARAMETER :: au2fs       = 0.02418885_wp       ! 1 a.u. = 2.418885e-2 fs (Time)
    REAL(wp), PARAMETER :: au2u        = 0.00054857994_wp    ! 1 a.u. = 5.4857994e-4 u (Mass)
                                                             !   1 u    = 1.660538782e-27 kg (unified atomic mass unit)
                                                             !   1 a.u. = 9.1093826e-31 kg (electron mass)
    REAL(wp), PARAMETER :: au2N        = 8.2388570d-8        ! 1 a.u. = 1 Hartree/Bohr = 4.359814e-8 N
    REAL(wp), PARAMETER :: au2J        = 4.3598140d-18       ! 1 a.u. = 1 Hartree = 4.3598140e-18 J

  END MODULE Constants_

  MODULE LAMMPSCalculator_
    USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_int32_t, c_int64_t, c_double
    USE LIBLAMMPS
    USE Precision_, ONLY: wp
    USE Constants_, ONLY: au2eV, au2kcalPmol
    USE Struct_,    ONLY: xyz_, abc_, force_
    IMPLICIT NONE

    TYPE lammps_calculator_
       TYPE(lammps)              :: lammps             ! LAMMPS instance
       REAL(wp)                  :: potential_energy   ! Potential energy (eV)
       TYPE(force_), ALLOCATABLE :: forces(:)          ! Forces on atoms  (eV/Ang)

       ! --- Add Kuwa 2025/02/25 ---
       integer, allocatable :: atom_id(:)
     CONTAINS
       PRIVATE
       PROCEDURE, PUBLIC :: run
       PROCEDURE, PUBLIC :: clear
       PROCEDURE, PUBLIC :: close
    END TYPE lammps_calculator_

  ! Local variables
    INTEGER(c_int32_t), POINTER, DIMENSION(:)   :: atom_id => NULL()  ! Atom ID
    REAL   (c_double) , POINTER, DIMENSION(:,:) :: xyz     => NULL()  ! Atomic coordinates   (Ang)
    REAL   (c_double) , POINTER                 :: energy  => NULL()  ! Potential energy     (kcal/mol)
    REAL   (c_double) , POINTER, DIMENSION(:,:) :: forces  => NULL()  ! Force acting on atom (kcal/mol/Ang)

    INTERFACE LAMMPSCalculator
       MODULE PROCEDURE constructor
    END INTERFACE LAMMPSCalculator

    PRIVATE :: wp
    PRIVATE :: au2eV, au2kcalPmol
    PRIVATE :: xyz_, abc_, force_
    PRIVATE :: xyz, energy, forces

  CONTAINS  !**************************************************************************************!

    TYPE(lammps_calculator_) FUNCTION constructor( lammps_file  &
                                                   ) RESULT(calculator)
      IMPLICIT NONE

      CHARACTER(LEN=*), OPTIONAL :: lammps_file

      INTEGER :: ierr, i
      INTEGER(KIND=c_int64_t), POINTER :: total_atoms

    ! For LAMMPS
      CHARACTER(LEN=16), PARAMETER :: command_line_args(5) =  &   !
                                       [ CHARACTER(LEN=16) :: 'lmp_nano', '-log', 'none', '-screen', 'none' ]

    ! I/O check
      ierr = 0
      IF ( PRESENT(lammps_file) .EQV. .FALSE. ) THEN
         WRITE(6,'(2X,A)') "'lammps_file' is mandatory in 'LAMMPSCalculator'!!!"
         ierr = 1
      END IF
      IF ( ierr == 1 ) STOP

    ! Create a LAMMPS instance (and initialize MPI)
      calculator%lammps = lammps( command_line_args(1:5) )

    ! Call lammps_file() to have LAMMPS read and process commands from a file
      CALL calculator%lammps%file( 'lammps.in' )

    ! Total number of atoms from LAMMPS
      total_atoms = calculator%lammps%extract_global( 'natoms' )

    ! Link atom ID inside LAMMPS with a pointer variable `atom_id'
      atom_id = calculator%lammps%extract_atom( 'id' )

    ! Link atomic coordinates inside LAMMPS with a pointer variable `xyz'
      xyz = calculator%lammps%extract_atom( 'x' )

    ! Link potential energy (kcal/mol) inside LAMMPS with a pointer variable `energy'
      energy = calculator%lammps%extract_compute( 'thermo_pe', 0, 0 )

    ! Link forces (kcal/mol/Ang) inside LAMMPS with a pointer variable `forces'
      forces = calculator%lammps%extract_atom( 'f' )

    ! For check
      IF ( .FALSE. ) THEN
         WRITE(6,*) energy
         DO i = 1, total_atoms
            WRITE(6,'(2X,I6,2X,1P,3(E14.7,X),2X,3(E14.7,X),I6)')  &
                 i, xyz(1,i), xyz(2,i), xyz(3,i), forces(1,i), forces(2,i), forces(3,i), atom_id(i)
         END DO
!         STOP
      END IF
      
      RETURN
    END FUNCTION constructor

  !================================================================================================!

    SUBROUTINE run( self, cartesian_coordinates )
      IMPLICIT NONE

      CLASS(lammps_calculator_), INTENT(INOUT) :: self
      TYPE(xyz_),                OPTIONAL      :: cartesian_coordinates(:)

      INTEGER  :: ierr, i
      INTEGER  :: num_atoms
      REAL(wp) :: coefficient

    ! I/O check
      ierr = 0
      IF ( PRESENT(cartesian_coordinates) .EQV. .TRUE. ) THEN
         num_atoms = SIZE( cartesian_coordinates(1:) )
      ELSE
         WRITE(6,'(2X,A)') "'cartesian_coordinates' is mandatory in 'LAMMPSCalculator%run'!!!"
         ierr = 1
      END IF
      IF ( ierr == 1 ) STOP

    !-----------------------------------------------------------------------------------------!

    ! Update new atomic coordinates (Ang)
      DO i = 1, num_atoms
         xyz(1,i) = cartesian_coordinates(atom_id(i))%x
         xyz(2,i) = cartesian_coordinates(atom_id(i))%y
         xyz(3,i) = cartesian_coordinates(atom_id(i))%z
      END DO

    ! For check
      IF ( .FALSE. ) THEN
         DO i = 1, num_atoms
            WRITE(6,'(2X,I6,3(F15.6,X),2X,3(F15.6,X),I6)')  &
                 i, xyz(1,i), xyz(2,i), xyz(3,i),           &
                 cartesian_coordinates(atom_id(i))%x,       &
                 cartesian_coordinates(atom_id(i))%y,       &
                 cartesian_coordinates(atom_id(i))%z,       &
                 atom_id(i)
         END DO
!         STOP
      END IF

    ! Execute LAMMPS with the atomic configuration
      CALL self%lammps%command( 'run 0' )

    ! For check
      IF ( .FALSE. ) THEN
         WRITE(6,*) energy
         DO i = 1, num_atoms
            WRITE(6,'(2X,I6,2X,1P,3(E14.7,X),2X,3(E14.7,X),I6)')  &
                 i, xyz(1,i), xyz(2,i), xyz(3,i),                 &
                 forces(1,i), forces(2,i), forces(3,i),           &
                 atom_id(i)
         END DO
!         STOP
      END IF

    !-----------------------------------------------------------------------------------------!

    ! Unit conversion of potential energy: kcal/mol -> eV
      coefficient = 1.0_wp / au2kcalPmol * au2eV
      self%potential_energy = energy * coefficient

    ! Unit conversion of forces: kcal/mol/Ang -> eV/Ang
      coefficient = coefficient
      IF ( ALLOCATED( self%forces ) .EQV. .FALSE. ) ALLOCATE( self%forces(1:num_atoms) )
      DO i = 1, num_atoms
         self%forces(atom_id(i))%x    = forces(1,i) * coefficient
         self%forces(atom_id(i))%y    = forces(2,i) * coefficient
         self%forces(atom_id(i))%z    = forces(3,i) * coefficient
         self%forces(atom_id(i))%norm = SQRT( forces(1,i)**2  &
                                            + forces(2,i)**2  &
                                            + forces(3,i)**2 ) * coefficient
      END DO

      IF ( ALLOCATED( self%atom_id ) .EQV. .FALSE. ) ALLOCATE( self%atom_id(1:num_atoms) )
      do i = 1, num_atoms
      self%atom_id(i) = atom_id(i)
      end do

      RETURN
    END SUBROUTINE run
    
  !================================================================================================!

    SUBROUTINE clear( self )
      IMPLICIT NONE
      CLASS(lammps_calculator_), INTENT(INOUT) :: self

      IF ( ALLOCATED( self%forces ) .EQV. .TRUE. ) DEALLOCATE( self%forces )

      RETURN
    END SUBROUTINE clear

  !================================================================================================!

    SUBROUTINE close( self )
      IMPLICIT NONE
      CLASS(lammps_calculator_), INTENT(INOUT) :: self

      NULLIFY( atom_id, xyz, energy, forces )
      CALL clear( self )
      CALL self%lammps%close( .TRUE. )

      RETURN
    END SUBROUTINE close

  END MODULE LAMMPSCalculator_
