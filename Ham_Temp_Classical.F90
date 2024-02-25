Subroutine Ham_Temp_Classical
  use Parameters
  use utility, only: norm_seq, outer_product
  implicit none
  integer :: iatom, imode, icolor, inhc, i, j
  real(8) :: get_kinetic_ene!, temp


  !dkinetic = 0.0d0
  !call Kinetic_Energy(dkinetic)
  !dkinetic = 0.5d0*dkinetic
  dkinetic = get_kinetic_ene()
  temp = 2.d0*dkinetic/dble(natom)/boltz/3.d0
  temp = temp/dble(nbead)


!block
!  real(8) :: moment(3)
!  moment(:) = 0.0d0
!  do i = 1, Natom
!    do j = 1, Nbead
!      moment(:) = moment(:) + physmass(i) * outer_product(ur(:,i,j),vur(:,i,j))
!    end do
!  end do
!  print *, istepsv, moment(:)
!end block


!
!   /*  calculate the total hamiltonian  */
!
  hamiltonian = dkinetic + potential
!
!   /*  sum bath variables  */
!
  !If(Iensemble==0) goto 100
!
!   /*  centroid thermostat  */
!
!YK First qmass_cent differs from others!!!
  !If(Iensemble==0.and.Simulation==3) goto 100

  !if(NCent==1) then
  select case(Ncent)
    case(0)
      continue
    case(1)
      ebath_cent=0.d0
      ebath_cent=ebath_cent + 0.5d0*qmcent11(1)*vbc11(1)*vbc11(1) + gnkt*rbc11(1)
      do imode=2,nnhc
        ebath_cent=ebath_cent + 0.5d0*qmcent11(imode)*vbc11(imode)*vbc11(imode) &
                              + gkt*rbc11(imode)
      enddo
      hamiltonian = hamiltonian + ebath_cent
  !endif
  !if(NCent==3) then
    case(3)
      ebath_cent=0.d0
      do inhc=1,nnhc
         do iatom=1,natom
           ebath_cent=ebath_cent + 0.5d0 * qmcent31(inhc) * norm_seq( vrbc31(:,iatom,inhc) ) &
                                 + gkt * sum( rbc31(:,iatom,inhc) )
         enddo
      enddo
      hamiltonian = hamiltonian + ebath_cent
  !endif
  end select

!     /*  total pseudo-hamiltonian  */
 !100  Continue
!YK include new centroid thermostat energies
!    hamiltonian = hamiltonian + ebath_cent
!YK

  Return
End Subroutine

