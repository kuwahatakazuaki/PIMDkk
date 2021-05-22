Subroutine Temp_ctr

  Use Parameters
  Implicit None
  Double Precision :: temp_scale, tempi

  tempi=0.0D+00
  call kinetic_energy(tempi)
  tempi = 0.5d0*tempi
  tempi = 2.d0*tempi/(3.d0*dble(natom))/boltz
  tempi = tempi/dble(nbead)
!    Write(*,*) 'Temperature calculated from Initial Velocities = ',tempi      
!    Write(*,*) 'Rescale the Velocities to Set to Temperature   = ',temperature

  temp_scale = dsqrt(temperature/tempi)
!      Write(*,*) 'Rescale Paramter   = ',temp_scale
  do imode = 1, nbead
     do iatom = 1, natom
        vux(iatom,imode) = vux(iatom,imode)*temp_scale
        vuy(iatom,imode) = vuy(iatom,imode)*temp_scale
        vuz(iatom,imode) = vuz(iatom,imode)*temp_scale
     enddo
  enddo

Return
End Subroutine
