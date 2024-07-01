subroutine Ham_Temp
  use Parameters
  use utility, only: norm_seq
  Implicit None
  Double Precision                   :: qdummy, factqk, dkin
  integer :: i, j, k, iatom, imode, inhc, icolor
  real(8) :: get_kinetic_ene!, temp

  !dkin  = 0.d0
  !call Kinetic_Energy(dkin)
  !dkinetic = dkin
  !dkinetic = 0.5d0*dkinetic
  dkinetic = get_kinetic_ene()
  temp = 2.d0*dkinetic/dble(natom)/KtoAU/3.d0
  temp = temp/dble(nbead)

!
!   pressure (PV) from virial and kinetic energy 
!   Be careful, we use f not du/dr
!

  virial = 0.0d0
  do j = 1, nbead
    do i = 1, natom
      virial = virial + dot_product( fr(:,i,j) , r(:,i,j) )
    end do
  end do
  PV  = (dkin + virial)/3.0d0 * AUtoJ / dble(nbead)

!
!   /*  quantum kinetic energy (harmonic interaction):  *
!    *  primitive estimator                             */
!
  qkinetic = 0.d0
  do imode = 2, nbead
    do iatom = 1, natom
      factqk = 0.5d0*dnmmass(iatom,imode)*omega_p2
      qkinetic = qkinetic + factqk * norm_seq( ur(:,iatom,imode) )
    enddo
  enddo
!
!   /*  calculate the total hamiltonian  */
!
    hamiltonian = dkinetic + potential + qkinetic
!
!   /*  sum bath variables  */
!
  ebath = 0.d0
  ebath_cent=0.d0
  if ( Ncent > 0 ) then
  !YK removed i=1 since this is not used
  ! and the first qmass_cent differs from others do i = 1, nbead
    do imode = 2, nbead
      qdummy = qmass(imode)
      do inhc = 1, nnhc
        do iatom = 1, natom
          ebath = ebath                                                   &
                + 0.5d0 * qdummy * norm_seq( vrbath(:,iatom,inhc,imode) ) &
                + gkt * sum( rbath(:,iatom,inhc,imode) )
        enddo
      enddo
    enddo
  end if
!
!     /*  centroid thermostat  */
!
!YK First qmass_cent differs from others!!!

  select case(Ncent)
    case(1)
      ebath_cent=ebath_cent + 0.5d0*qmcent11(1)*vbc11(1)*vbc11(1) + gnkt*rbc11(1)
      do imode=2,nnhc
        ebath_cent=ebath_cent + 0.5d0*qmcent11(imode)*vbc11(imode)*vbc11(imode) + gkt*rbc11(imode)
      enddo
    case(3)
    do inhc=1,nnhc
      do iatom=1,natom
        ebath_cent=ebath_cent + 0.5d0 * qmcent31(inhc) * norm_seq( vrbc31(:,iatom,inhc) ) &
                              + gkt * sum( rbc31(:,iatom,inhc) )
      enddo
    enddo
  end select
!YK include new centroid thermostat energies
    hamiltonian = hamiltonian + ebath + ebath_cent

  Call Virial_Estimator

  return
end subroutine Ham_Temp
