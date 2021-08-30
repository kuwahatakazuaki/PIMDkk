subroutine calc_hamil
  use global_variable
  use utility
  implicit none
!  Double Precision                   :: qdummy, factqk, dkin
  integer :: i, j, inhc, icolor
  real(8) :: dkinetic, qkinetic, hamiltonian, potential
  real(8) :: ebath, ebath_cent, virial, primitive, steptemp

  dkinetic = kinetic_energy()

  steptemp = 2.0d0 * dkinetic / ( dble(3*Natom*Nbead)*boltz )

  qkinetic = 0.0d0
  do j = 1, Nbead
    do i = 1, Natom
!      factqk = 0.5d0 * dnmmass(i,j) * omega_p2
!      qkinetic = qkinetic + factqk * dot_product(u(:,i,j),u(:,i,j))
      qkinetic = qkinetic + 0.5d0*dnmmass(i,j)*omega_p2 * dot_product(u(:,i,j),u(:,i,j))
    end do
  end do

  potential = sum(energy(:)) / dble(Nbead)


  ebath      = 0.0d0
  ebath_cent = 0.0d0


!  gnkt is 3.0d0 * dble(Natom) / beta
!  gkt  is 1.0d0 / beta
  if ( Nensemble == 0 ) goto 100
  do j = 2, Nbead
    do inhc = 1, Nnhc
      do i = 1, Natom
        ebath = ebath + 0.5d0 * qmass(j) * dot_product(vbath(:,i,inhc,j),vbath(:,i,inhc,j)) &
                      + gkt * dot_product(bath(:,i,inhc,j),bath(:,i,inhc,j))
      end do
    end do
  end do

  if ( Simulation == 3 ) goto 100   ! skip in CMD ??

  if ( Ncent == 1) then

    if ( Ncolor == 1 ) then
      ebath_cent = ebath_cent + 0.5 * qmcent11(1)*vbc11(1)*vbc11(1) + gkt*rbc11(1)
      do inhc = 2, Nnhc
        ebath_cent = ebath_cent + 0.5 * qmcent11(inhc)*vbc11(inhc)*vbc11(inhc) +gkt*rbc11(inhc)
      end do
    else if ( Ncolor >= 2 ) then
      do icolor = 1, Ncolor
        ebath_cent = ebath_cent + 0.5 *qmcent1(1,icolor)*vbc1(1,icolor)*vbc1(1,icolor) + gnkt*rbc1(1,icolor)
        do inhc = 2, Nnhc
          ebath_cent = ebath_cent + 0.5 *qmcent1(inhc,icolor)*vbc1(inhc,icolor)*vbc1(inhc,icolor) + gnkt*rbc1(inhc,icolor)
        end do
      end do
    end if

  else if ( Ncent == 3) then
    if ( Ncolor == 1 ) then
      do inhc = 1, Nnhc
        do i = 1, Natom
          ebath_cent = ebath_cent + 0.5*qmcent31(inhc)*dot_product(vbc31(:,i,inhc),vbc31(:,i,inhc)) &
                                  + gkt * sum(rbc31(:,i,inhc))
        end do
      end do
    else if ( Ncolor >= 2 ) then
      do icolor = 1, Ncolor
        do inhc = 1, Nnhc
          do i = 1, Natom
            ebath_cent = ebath_cent &
                        + 0.5*qmcent3(inhc,icolor)*dot_product(vbc3(:,i,inhc,icolor),vbc3(:,i,inhc,icolor)) &
                        + gkt * sum(rbc3(:,i,inhc,icolor))
          end do
        end do
      end do
    end if
  end if

100 continue

  hamiltonian = dkinetic + qkinetic + potential + ebath + ebath_cent
!  print *, dkinetic
!  print *, qkinetic
!  print *, potential
!  print *, ebath
!  print *, ebath_cent


!    Call Virial_Estimator
    call calc_estimator
    call print_ham

return
contains

  subroutine calc_estimator
    real(8) :: temp
    temp = 0.0d0
    do j = 1, Nbead
      do i = 1, Natom
        temp = temp + dot_product( f(:,i,j), (r(:,i,j)-u(:,i,1)) )
      end do
    end do
    virial = 1.5d0 * dble(Natom)/beta - 0.5d0 * temp

    primitive = 0.0d0
  end subroutine calc_estimator

  subroutine print_ham
    use global_variable
    integer :: iham
    open(newunit=iham,file=path_result//'/ham.dat',position='append')
      write(iham,9998) Istep, hamiltonian, potential, dkinetic, qkinetic, ebath, ebath_cent, steptemp, virial
    close(iham)

    if ( mod(Istep,100) == 0 ) then
      open(newunit=Uout,file=Oname,status='old',position='append')
        write(Uout,9999) Istep, hamiltonian, potential, dkinetic, qkinetic, ebath, ebath_cent, steptemp, virial
      close(Uout)
    end if
    return
9999 format(i7,8e17.9)
9998 format(i7,8e23.15)
  end subroutine print_ham

end subroutine calc_hamil

!
!    dkin  = 0.d0
!    call Kinetic_Energy(dkin)
!    dkinetic = dkin
!    dkinetic = 0.5d0*dkinetic
!    temp = 2.d0*dkinetic/dble(natom)/boltz/3.d0
!    temp = temp/dble(nbead)
!!
!!   /*  quantum kinetic energy (harmonic interaction):  *
!!    *  primitive estimator                             */
!!
!    qkinetic = 0.d0
!    do imode = 1, nbead
!       do iatom = 1, natom
!          factqk = 0.5d0*dnmmass(iatom,imode)*omega_p2
!          qkinetic = qkinetic                               &
!                   + factqk*ux(iatom,imode)*ux(iatom,imode) &
!                   + factqk*uy(iatom,imode)*uy(iatom,imode) &
!                   + factqk*uz(iatom,imode)*uz(iatom,imode)
!       enddo
!    enddo
!!
!!   /*  calculate the total hamiltonian  */
!!
!    hamiltonian = dkinetic + potential + qkinetic
!!!
!!   /*  sum bath variables  */
!!
!    ebath = 0.d0
!    ebath_cent=0.d0
!    If(NEnsemble==0.and.Simulation==5) goto 100
!!YK removed i=1 since this is not used
!!   and the first qmass_cent differs from others
!!   do i = 1, nbead
!    do imode = 2, nbead
!!YK
!       qdummy = qmass(imode)
!       do inhc = 1, nnhc
!          do iatom = 1, natom
!             ebath = ebath                                                           &
!                   + 0.5d0*qdummy*vxbath(iatom,inhc,imode)*vxbath(iatom,inhc,imode)  &
!                   + 0.5d0*qdummy*vybath(iatom,inhc,imode)*vybath(iatom,inhc,imode)  &
!                   + 0.5d0*qdummy*vzbath(iatom,inhc,imode)*vzbath(iatom,inhc,imode)  &
!                   + gkt*xbath(iatom,inhc,imode)                                     &
!                   + gkt*ybath(iatom,inhc,imode)                                     &
!                   + gkt*zbath(iatom,inhc,imode)
!          enddo
!       enddo
!    enddo

!!
!!     /*  centroid thermostat  */
!!
!!YK First qmass_cent differs from others!!!
!
!    If(NEnsemble==0.and.Simulation==3) goto 100
!
!    if(NCent==1) then
!
!       If(NColor==1) Then
!! qmcent11(1) = 3.d0*dble(natom)/beta/omega2
!          ebath_cent=ebath_cent + 0.5d0*qmcent11(1)*vbc11(1)*vbc11(1) &
!                                + gnkt*rbc11(1)
!          do imode=2,nnhc
!             ebath_cent=ebath_cent + 0.5d0*qmcent11(imode)*vbc11(imode)*vbc11(imode) &
!                                   + gkt*rbc11(imode)
!          enddo
!
!       Else
!
!          do icolor=1,ncolor
!             ebath_cent=ebath_cent + 0.5d0*qmcent1(1,icolor)*vbc1(1,icolor)*vbc1(1,icolor) &
!                                   + gnkt*rbc1(1,icolor)
!             do imode=2,nnhc
!                ebath_cent=ebath_cent + 0.5d0*qmcent1(imode,icolor)*vbc1(imode,icolor)*vbc1(imode,icolor) &
!                                      + gkt*rbc1(imode,icolor)
!             enddo
!          enddo
!
!       EndIf
!
!    endif
!
!    if(NCent==3) then
!
!       If(NColor==1) Then
!
!          do inhc=1,nnhc
!             do iatom=1,natom
!                ebath_cent=ebath_cent + 0.5d0*qmcent31(inhc)*vxbc31(iatom,inhc)*vxbc31(iatom,inhc) &
!                                      + 0.5d0*qmcent31(inhc)*vybc31(iatom,inhc)*vybc31(iatom,inhc) &
!                                      + 0.5d0*qmcent31(inhc)*vzbc31(iatom,inhc)*vzbc31(iatom,inhc) &
!                                      + gkt*xbc31(iatom,inhc)                                      &
!                                      + gkt*ybc31(iatom,inhc)                                      &
!                                      + gkt*zbc31(iatom,inhc)                                        
!             enddo
!          enddo
!
!       Else
!
!          do icolor=1,ncolor
!             do inhc=1,nnhc
!                do iatom=1,natom
!                   ebath_cent=ebath_cent &
!                             + 0.5d0*qmcent3(inhc,icolor)*vxbc3(iatom,inhc,icolor)*vxbc3(iatom,inhc,iatom) &
!                             + 0.5d0*qmcent3(inhc,icolor)*vybc3(iatom,inhc,icolor)*vybc3(iatom,inhc,iatom) &
!                             + 0.5d0*qmcent3(inhc,icolor)*vzbc3(iatom,inhc,icolor)*vzbc3(iatom,inhc,iatom) &
!                             + gkt*xbc3(iatom,inhc,icolor)                              &
!                             + gkt*ybc3(iatom,inhc,icolor)                              &
!                             + gkt*zbc3(iatom,inhc,icolor)                                        
!                enddo
!             enddo
!          enddo
!
!       EndIf
!
!    endif

