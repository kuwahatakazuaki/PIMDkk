Subroutine Ham_Temp_Classical

  Use Parameters
  Implicit None

  Integer                            :: i, j, k, ij = 0
  Double Precision                   :: qdummy, factqk, dkin


    dkin  = 0.d0
    call Kinetic_Energy(dkin)
    dkinetic = dkin
    dkinetic = 0.5d0*dkinetic
    temp = 2.d0*dkinetic/dble(natom)/boltz/3.d0
    temp = temp/dble(nbead)
!
!     /*  calculate the total hamiltonian  */
!
      hamiltonian = dkinetic + potential
!
!     /*  sum bath variables  */
!
      ebath_cent=0.d0
      If(NEnsemble==0) goto 100
!
!     /*  centroid thermostat  */
!
!YK First qmass_cent differs from others!!!

    If(NEnsemble==0.and.Simulation==3) goto 100
    if(NCent==1) then

       If(NColor==1) Then

          ebath_cent=ebath_cent + 0.5d0*qmcent11(1)*vbc11(1)*vbc11(1) &
                                + gnkt*rbc11(1)
          do imode=2,nnhc
             ebath_cent=ebath_cent + 0.5d0*qmcent11(imode)*vbc11(imode)*vbc11(imode) &
                                   + gkt*rbc11(imode)
          enddo

       Else
          do icolor=1,ncolor
             ebath_cent=ebath_cent + 0.5d0*qmcent1(1,icolor)*vbc1(1,icolor)*vbc1(1,icolor) &
                                   + gnkt*rbc1(1,icolor)
             do imode=2,nnhc
                ebath_cent=ebath_cent + 0.5d0*qmcent1(imode,icolor)*vbc1(imode,icolor)*vbc1(imode,icolor) &
                                      + gkt*rbc1(imode,icolor)
             enddo
          enddo
       EndIf
    endif

    if(NCent==3) then
       If(NColor==1) Then
          do inhc=1,nnhc
             do iatom=1,natom
                ebath_cent=ebath_cent + 0.5d0*qmcent31(inhc)*vxbc31(iatom,inhc)*vxbc31(iatom,inhc) &
                                      + 0.5d0*qmcent31(inhc)*vybc31(iatom,inhc)*vybc31(iatom,inhc) &
                                      + 0.5d0*qmcent31(inhc)*vzbc31(iatom,inhc)*vzbc31(iatom,inhc) &
                                      + gkt*xbc31(iatom,inhc)                                      &
                                      + gkt*ybc31(iatom,inhc)                                      &
                                      + gkt*zbc31(iatom,inhc)                                        
             enddo
          enddo
       Else

          do icolor=1,ncolor
             do inhc=1,nnhc
                do iatom=1,natom
                   ebath_cent=ebath_cent &
                             + 0.5d0*qmcent3(inhc,icolor)*vxbc3(iatom,inhc,icolor)*vxbc3(iatom,inhc,iatom) &
                             + 0.5d0*qmcent3(inhc,icolor)*vybc3(iatom,inhc,icolor)*vybc3(iatom,inhc,iatom) &
                             + 0.5d0*qmcent3(inhc,icolor)*vzbc3(iatom,inhc,icolor)*vzbc3(iatom,inhc,iatom) &
                             + gkt*xbc3(iatom,inhc,icolor)                              &
                             + gkt*ybc3(iatom,inhc,icolor)                              &
                             + gkt*zbc3(iatom,inhc,icolor)                                        
                enddo
             enddo
          enddo

       EndIf

    endif
!!    dummy = 0.5d0*qmass(1)
!!    ebath_cent = dummy*vbathcent(1)*vbathcent(1) + gnkt*rbathcent(1)
!!    do i = 2, nnhc
!!      ebath_cent = ebath_cent + dummy*vbathcent(i)*vbathcent(i) + gkt*rbathcent(i)
!!    enddo
!
!     /*  total pseudo-hamiltonian  */
!
!!    hamiltonian = hamiltonian + ebath + ebath_cent 
 100  Continue
!YK include new centroid thermostat energies
    hamiltonian = hamiltonian + ebath_cent
!      hamiltonian = hamiltonian + ebath 
!YK

    Return
  End Subroutine

