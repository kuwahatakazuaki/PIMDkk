Subroutine Ham_Temp_Classical
  use Parameters
  use utility, only: norm_seq
  implicit none
  integer :: iatom, icolor, Inhc, i, j
  real(8) :: get_kinetic_ene

  dkinetic = get_kinetic_ene()
  temp = 2.d0*dkinetic*AU2K/dble(Ndim*Natom*Nbead)
  !temp = 2.d0*dkinetic*AU2K/dble(3*Natom*Nbead)

!
!   /*  calculate the total hamiltonian  */
!
  hamiltonian = dkinetic + potential
!
!   /*  sum bath variables  */
!   /*  centroid thermostat  */
!
!YK First qmass_cent differs from others!!!
  ebath_cent=0.d0
  select case(Ncent)
    case(0)
      continue
    case(1)
      ebath_cent=ebath_cent + 0.5d0*qmcent11(1)*vbc11(1)*vbc11(1) + gnkt*rbc11(1)
      do Inhc=2,nnhc
        ebath_cent=ebath_cent + 0.5d0*qmcent11(Inhc)*vbc11(Inhc)*vbc11(Inhc) &
                              + gkt*rbc11(Inhc)
      enddo
    case(3)
      do inhc=1,Nnhc
         do iatom=1,Natom
           ebath_cent=ebath_cent + 0.5d0 * qmcent31(inhc) * norm_seq( vrbc31(:,iatom,inhc) ) &
                                 + gkt * sum( rbc31(:,iatom,inhc) )
         enddo
      enddo
  end select
  hamiltonian = hamiltonian + ebath_cent

!     /*  total pseudo-hamiltonian  */
 !100  Continue
!YK include new centroid thermostat energies
!    hamiltonian = hamiltonian + ebath_cent
!YK

  return
end subroutine Ham_Temp_Classical

