Subroutine Calc_Constant
  Use Parameters
  Implicit None
  integer :: iatom

!   twopi=2.d0*pi
  !! fs  -- > A.U.
  dt = dt * facttime
  !! amu -- > A.U.
  !do iatom=1,natom
  !   physmass(iatom) = physmass(iatom)*factmass
  !enddo
  physmass(1:natom) = physmass(1:natom)*factmass


  beta   = 1.d0/(boltz*temperature)
  !! Frequency of the system
  omega_system = 2.d0*pi/freq1
  omega_system = omega_system/facttime

  gnkt = 3.d0*dble(natom)/beta
  gkt = 1.d0/beta

  address0=trim(address2)//'/'
  laddress=len_trim(address2)+1
  addresstmp=trim(trim(address2)//'/')

!
!     /*   parameters for path integral simulation   */
!
    dp = dble(nbead)
    dp_inv = 1.0d0/dp
    omega_p = dsqrt(dp)/beta
    omega_p2 = omega_p*omega_p
!    omega2 = omega_system*omega_system/dp
    omega2 = omega_system*omega_system

  select case(Isimulation)
    case(0)
      simulation = "PIMD"
    case(1)
      simulation = "RPMD"
    case(2)
      simulation = "CMD"
    case(10)
      simulation = "MD"
  end select

  Return
End Subroutine