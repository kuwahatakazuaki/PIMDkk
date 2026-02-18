subroutine Calc_Constant
  use Parameters
  implicit none
  integer :: iatom
  real(8) :: omega_p, dp

  dt = dt*fs2AU   ! fs -> AU
  physmass(1:natom) = physmass(1:natom)*amu2AU  ! amu -> AU

  beta   = 1.d0/(temperature*K2AU)
  !! Frequency of the system
  omega_system = 2.d0*pi/(freq1*fs2AU)

  select case(Ndim)
    case(3)

      !if ( Natom >= 3 ) then
      !  Ndof = 3*Natom - 6
      !else if ( Natom == 2 ) then
      !  Ndof = 3*Natom - 5
      !else if ( Natom == 1 ) then
      !  Ndof = 3*Natom
      !end if
      Ndof = 3*Natom

    case(1)

      if ( Natom >= 2 ) then
        !Ndof = Natom
        Ndof = Natom - 1
      else if ( Natom == 1 ) then
        Ndof = Natom
      end if

    case default
      stop 'ERROR!! Wrong "Ndim"'

  end select

  gnkt = dble(Ndof)/beta
  gkt = 1.d0/beta

  !gnkt = 3.d0*dble(Natom)/beta
  !gkt = 1.d0/beta

  address0=trim(dir_scr)//'/'
  laddress=len_trim(dir_scr)+1
  addresstmp=trim(trim(dir_scr)//'/')

!
!     /*   parameters for path integral simulation   */
!
    dp = dble(Nbead)
    dp_inv = 1.0d0/dp
    omega_p2 = dp/(beta*beta)
!    omega2 = omega_system*omega_system/dp
    omega2 = omega_system*omega_system

  return
end subroutine Calc_Constant
