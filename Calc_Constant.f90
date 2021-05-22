Subroutine Calc_Constant

  Use Parameters
  Implicit None

  pi=3.14159265358979d0
!   twopi=2.d0*pi
  !! fs  -- > A.U.
  dt = dt * facttime
  !! amu -- > A.U.
  do iatom=1,natom
     physmass(iatom) = physmass(iatom)*factmass
  enddo
!   physmass(1:natom) = physmass(1:natom)*factmass
  !!
  beta   = 1.d0/(boltz*temperature)
  !! Frequency of the system
  omega_system = 2.d0*pi/freq1
  omega_system = omega_system/facttime

  If (Simulation==1 ) Then
     omega_system = 1.d0
     beta         = 5.d0
     physmass (1:natom) = 1.d0
     temperature  = 1.d0 / (boltz*beta)
!      dt = 0.003d0 * facttime
  Endif

  Return
End Subroutine
