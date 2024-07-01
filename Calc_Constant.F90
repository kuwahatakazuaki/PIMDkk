Subroutine Calc_Constant
  use Parameters
  Implicit None
  integer :: iatom
  real(8) :: omega_p, dp

  !! fs  -- > A.U.
  dt = dt * facttime
  !! amu -- > A.U.
  physmass(1:natom) = physmass(1:natom)*factmass


  beta   = 1.d0/(KtoAU*temperature)
  !! Frequency of the system
  omega_system = 2.d0*pi/(freq1*facttime)

  gnkt = 3.d0*dble(Natom)/beta
  gkt = 1.d0/beta

  address0=trim(address2)//'/'
  laddress=len_trim(address2)+1
  addresstmp=trim(trim(address2)//'/')

!
!     /*   parameters for path integral simulation   */
!
    dp = dble(Nbead)
    dp_inv = 1.0d0/dp
    omega_p2 = dp/(beta*beta)
!    omega2 = omega_system*omega_system/dp
    omega2 = omega_system*omega_system

  !select case(Isimulation)
  !  case(0)
  !    simulation = "PIMD"
  !  case(1)
  !    simulation = "RPMD"
  !  case(2)
  !    simulation = "CMD"
  !  case(10)
  !    simulation = "MD"
  !end select

  Return
End Subroutine
