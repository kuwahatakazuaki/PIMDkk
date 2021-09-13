subroutine integ_nhc_cent
  use global_variable, only: Natom, Ncent, Ncolor, Nnhc, Nys, &
      r, vu, fictmass, &
      vbc11, fbc11, &
      qmcent11, &
      gnkt
  use utility

    select case(Ncent)
      case(1)
        call integ_nhc_cent1
      case(3)
        call integ_nhc_cent3
      case default
        call program_abort('Erro!!! Wrong "Ncent"')
    end select
  contains

! *** Ncent = 1 ***
  subroutine integ_nhc_cent1
    real(8) :: skin
    integer :: i, inhc, iys

    skin = 0.0d0
    do i = 1, Natom
      skin = skin + fictmass(i,1) * dot_product(vu(:,i,1),vu(:,i,1))
    end do

    if ( Ncolor == 1 ) then
      fbc11(1) = (skin-gnkt) / qmcent11(1)
      do inhc = 2, Nnhc
        fbc11(inhc) = (qmcent11(inhc-1)*vbc11(inhc-1)*vbc11(inhc-1) - gkt)/qmcent11(inhc)
      end do

      do iys = 1, Nys
      end do
    else if ( Ncolor > 1 ) then
    end if

  end subroutine integ_nhc_cent1
! *** End Ncent = 1 ***

! *** Ncent = 3 ***
  subroutine integ_nhc_cent3
  end subroutine integ_nhc_cent3
! *** End Ncent = 3 ***

end subroutine integ_nhc_cent


!  Subroutine Nhc_Integrate_Cent
!
!    Use Parameters
!    Implicit None
!
!    Double Precision                   :: skin, scale, dt_ys, vfact, pvfact
!    Integer                            :: i, iys
!
!    skin = 0.d0
!    do iatom = 1, natom
!       skin = skin                                        & 
!            + fictmass(iatom,1)*vux(iatom,1)*vux(iatom,1) &
!            + fictmass(iatom,1)*vuy(iatom,1)*vuy(iatom,1) &
!            + fictmass(iatom,1)*vuz(iatom,1)*vuz(iatom,1)
!    enddo
!
!    scale = 1.d0
!
!!     /* update the force */
!!YK Rewrote using qmass_cent
!
!    If(NColor==1) Then
!
!       fbc11(1) = (skin - gnkt)/qmcent11(1)
!       do inhc = 2, nnhc
!          fbc11(inhc) = (qmcent11(inhc-1)*vbc11(inhc-1)*vbc11(inhc-1) &
!                      - gkt)/qmcent11(inhc)
!       enddo
!
!!     /*  start multiple time step integration  */
!       do iys = 1, nys 
!
!!     /*  set time increment at this loop  */
!          dt_ys = dt*ysweight(iys)
!
!!     /* update the thermostat velocities */
!          vbc11(nnhc) = vbc11(nnhc) + 0.25d0*fbc11(nnhc)*dt_ys
!          do inhc = 1, nnhc-1
!             vfact=dexp(-0.125d0*   vbc11(nnhc-inhc+1)*dt_ys)
!             vbc11(nnhc-inhc) = vbc11(nnhc-inhc)*vfact*vfact &
!                              + 0.25d0*fbc11(nnhc-inhc)*vfact*dt_ys
!          enddo
!
!!     /* update the particle velocities */
!          pvfact = dexp(-0.5d0*vbc11(1)*dt_ys)
!          scale = scale*pvfact
!!     /* update the force */
!!YK changed to qmass_cent
!          fbc11(1)=(scale*scale*skin - gnkt)/qmcent11(1)
!!YK
!!     /* update the thermostat position */
!!          write(*,*) 'RBC'
!          do inhc = 1, nnhc
!             rbc11(inhc) = rbc11(inhc) &
!                         + 0.5d0*vbc11(inhc)*dt_ys
!          enddo
!!
!!     /* update the thermostat velocities */
!          do inhc = 1, nnhc-1
!             vfact = dexp(-0.125d0*vbc11(inhc+1)*dt_ys)
!             vbc11(inhc) = vbc11(inhc)*vfact*vfact &
!                         + 0.25d0*fbc11(inhc)*vfact*dt_ys
!!YK changed to qmass_cent
!             fbc11(inhc+1) = (qmcent11(inhc)*vbc11(inhc)*vbc11(inhc) &
!                           - gkt)/qmcent11(inhc+1)
!          enddo
!
!          vbc11(nnhc) = vbc11(nnhc) + 0.25d0*fbc11(nnhc)*dt_ys
!       enddo
!    
!    Else
!
!       do icolor=1,ncolor
!          fbc1(1,icolor) = (skin - gnkt)/qmcent1(1,icolor)
!          do inhc = 2, nnhc
!             fbc1(inhc,icolor) = (qmcent1(inhc-1,icolor)*vbc1(inhc-1,icolor)*vbc1(inhc-1,icolor) &
!                               - gkt)/qmcent1(inhc,icolor)
!          enddo
!       enddo   
!
!!     /*  start multiple time step integration  */
!       do icolor=1,ncolor
!          do iys = 1, nys 
!
!!     /*  set time increment at this loop  */
!             dt_ys = dt*ysweight(iys)
!
!!     /* update the thermostat velocities */
!             vbc1(nnhc,icolor) = vbc1(nnhc,icolor) + 0.25d0*fbc1(nnhc,icolor)*dt_ys
!             do inhc = 1, nnhc-1
!                vfact=dexp(-0.125d0*   vbc1(nnhc-inhc+1,icolor)*dt_ys)
!                vbc1(nnhc-inhc,icolor) = vbc1(nnhc-inhc,icolor)*vfact*vfact &
!                                       + 0.25d0*fbc1(nnhc-inhc,icolor)*vfact*dt_ys
!             enddo
!
!!     /* update the particle velocities */
!             pvfact = dexp(-0.5d0*vbc1(1,icolor)*dt_ys)
!             scale = scale*pvfact
!
!!     /* update the force */
!!YK changed to qmass_cent
!             fbc1(1,icolor)=(scale*scale*skin - gnkt)/qmcent1(1,icolor)
!!YK
!!     /* update the thermostat position */
!             do inhc = 1, nnhc
!                rbc1(inhc,icolor) = rbc1(inhc,icolor) &
!                                  + 0.5d0*vbc1(inhc,icolor)*dt_ys
!             enddo
!
!!     /* update the thermostat velocities */
!             do inhc = 1, nnhc-1
!                vfact = dexp(-0.125d0*vbc1(inhc+1,icolor)*dt_ys)
!                vbc1(inhc,icolor) = vbc1(inhc,icolor)*vfact*vfact &
!                                  + 0.25d0*fbc1(inhc,icolor)*vfact*dt_ys
!!YK changed to qmcent
!                fbc1(inhc+1,icolor) = &
!                                    + (qmcent1(inhc,icolor)*vbc1(inhc,icolor)*vbc1(inhc,icolor) &
!                                    -  gkt)/qmcent1(inhc+1,icolor)
!             enddo
!             vbc1(nnhc,icolor) = vbc1(nnhc,icolor) + 0.25d0*fbc1(nnhc,icolor)*dt_ys
!          enddo
!       enddo
!    EndIf
!
!!     /* update the paricle velocities */
!    do iatom = 1, natom
!       vux(iatom,1) = vux(iatom,1)*scale
!       vuy(iatom,1) = vuy(iatom,1)*scale
!       vuz(iatom,1) = vuz(iatom,1)*scale
!    enddo
!
!    Return
!  End Subroutine
!
!
!  Subroutine nhc_integrate_cent3
!
!    Use Parameters
!    Implicit None
!    Double Precision                   :: skin, scale, dt_ys, vfact, pvfact
!    Double Precision                   :: vxfact,vyfact,vzfact
!    Double Precision                   :: pvxfact,pvyfact,pvzfact
!    Double Precision                   :: scalex,scaley,scalez
!    Double Precision                   :: dkinx,dkiny,dkinz
!    Integer                            :: i, iys
!
!
!    If(NColor==1) Then
!
!!     /*  start multiple time step integration  */
!       do iys = 1, NYS
!!     /*  set time increment at this loop  */
!          dt_ys = dt*ysweight(iys)
!          do iatom = 1, NATOM
!
!!     /*  calculate kinetic energy of centroid coordiates  */
!!     /*  for each atom                                    */
!
!             dkinx = fictmass(iatom,1)*vux(iatom,1)*vux(iatom,1)
!             dkiny = fictmass(iatom,1)*vuy(iatom,1)*vuy(iatom,1)
!             dkinz = fictmass(iatom,1)*vuz(iatom,1)*vuz(iatom,1)
!             scalex = 1.d0
!             scaley = 1.d0
!             scalez = 1.d0
!
!!     /* update the force */
!
!             fxbc31(iatom,1) = (dkinx - gkt)/qmcent31(1)
!             fybc31(iatom,1) = (dkiny - gkt)/qmcent31(1)
!             fzbc31(iatom,1) = (dkinz - gkt)/qmcent31(1)
!   
!             do inhc = 2, NNHC
!                fxbc31(iatom,inhc) = &
!                                    (qmcent31(inhc-1)*vxbc31(iatom,inhc-1)*vxbc31(iatom,inhc-1) &
!                                   - gkt)/qmcent31(inhc)
!                fybc31(iatom,inhc) = &
!                                    (qmcent31(inhc-1)*vybc31(iatom,inhc-1)*vybc31(iatom,inhc-1) &
!                                   - gkt)/qmcent31(inhc)
!                fzbc31(iatom,inhc) = &
!                                    (qmcent31(inhc-1)*vzbc31(iatom,inhc-1)*vzbc31(iatom,inhc-1) &
!                                   - gkt)/qmcent31(inhc)
!             enddo
!
!!     /* update the thermostat velocities */
!
!             vxbc31(iatom,NNHC) = vxbc31(iatom,NNHC) &
!                                + 0.25d0*fxbc31(iatom,NNHC)*dt_ys
!             vybc31(iatom,NNHC) = vybc31(iatom,NNHC) &
!                                + 0.25d0*fybc31(iatom,NNHC)*dt_ys
!             vzbc31(iatom,NNHC) = vzbc31(iatom,NNHC) &
!                                + 0.25d0*fzbc31(iatom,NNHC)*dt_ys
!
!            do inhc = 1, NNHC-1
!               vxfact=dexp(-0.125d0*vxbc31(iatom,NNHC-inhc+1)*dt_ys)
!               vyfact=dexp(-0.125d0*vybc31(iatom,NNHC-inhc+1)*dt_ys)
!               vzfact=dexp(-0.125d0*vzbc31(iatom,NNHC-inhc+1)*dt_ys)
!               vxbc31(iatom,NNHC-inhc) = vxbc31(iatom,NNHC-inhc)*vxfact*vxfact &
!                                       + 0.25d0*fxbc31(iatom,NNHC-inhc)*vxfact*dt_ys   
!               vybc31(iatom,NNHC-inhc) = vybc31(iatom,NNHC-inhc)*vyfact*vyfact &
!                                       + 0.25d0*fybc31(iatom,NNHC-inhc)*vyfact*dt_ys
!               vzbc31(iatom,NNHC-inhc) = vzbc31(iatom,NNHC-inhc)*vzfact*vzfact &
!                                       + 0.25d0*fzbc31(iatom,NNHC-inhc)*vzfact*dt_ys
!            enddo
!
!!     /* update the particle velocities */
!
!             pvxfact = dexp(-0.5d0*vxbc31(iatom,1)*dt_ys)
!             pvyfact = dexp(-0.5d0*vybc31(iatom,1)*dt_ys)
!             pvzfact = dexp(-0.5d0*vzbc31(iatom,1)*dt_ys)
!   
!             scalex = scalex*pvxfact
!             scaley = scaley*pvyfact
!             scalez = scalez*pvzfact
!
!!     /* update the force */
!             fxbc31(iatom,1)=(scalex*scalex*dkinx - gkt)/qmcent31(1)
!             fybc31(iatom,1)=(scaley*scaley*dkiny - gkt)/qmcent31(1)
!             fzbc31(iatom,1)=(scalez*scalez*dkinz - gkt)/qmcent31(1)
!
!!     /* update the thermostat position */
!             do inhc = 1, NNHC
!                xbc31(iatom,inhc) = xbc31(iatom,inhc) &
!                                  + 0.5d0*vxbc31(iatom,inhc)*dt_ys
!                ybc31(iatom,inhc) = ybc31(iatom,inhc) &
!                                  + 0.5d0*vybc31(iatom,inhc)*dt_ys
!                zbc31(iatom,inhc) = zbc31(iatom,inhc) &
!                                  + 0.5d0*vzbc31(iatom,inhc)*dt_ys
!             enddo
!
!!     /* update the thermostat velocities */
!             do inhc = 1, NNHC-1
!                vxfact = dexp(-0.125d0*vxbc31(iatom,inhc+1)*dt_ys)
!                vyfact = dexp(-0.125d0*vybc31(iatom,inhc+1)*dt_ys)
!                vzfact = dexp(-0.125d0*vzbc31(iatom,inhc+1)*dt_ys)
!   
!                vxbc31(iatom,inhc) = vxbc31(iatom,inhc)*vxfact*vxfact &
!                                   + 0.25d0*fxbc31(iatom,inhc)*vxfact*dt_ys
!                vybc31(iatom,inhc) = vybc31(iatom,inhc)*vyfact*vyfact &
!                                   + 0.25d0*fybc31(iatom,inhc)*vyfact*dt_ys
!                vzbc31(iatom,inhc) = vzbc31(iatom,inhc)*vzfact*vzfact &
!                                   + 0.25d0*fzbc31(iatom,inhc)*vzfact*dt_ys
!   
!                fxbc31(iatom,inhc+1) = (qmcent31(inhc)*vxbc31(iatom,inhc)*vxbc31(iatom,inhc) &
!                                     - gkt)/qmcent31(inhc+1) 
!                fybc31(iatom,inhc+1) = (qmcent31(inhc)*vybc31(iatom,inhc)*vybc31(iatom,inhc) &
!                                     - gkt)/qmcent31(inhc+1)
!                fzbc31(iatom,inhc+1) = (qmcent31(inhc)*vzbc31(iatom,inhc)*vzbc31(iatom,inhc) &
!                                     - gkt)/qmcent31(inhc+1) 
!             enddo
!   
!             vxbc31(iatom,NNHC) = vxbc31(iatom,NNHC)  &
!                                + 0.25d0*fxbc31(iatom,NNHC)*dt_ys
!             vybc31(iatom,NNHC) = vybc31(iatom,NNHC)  &
!                                + 0.25d0*fybc31(iatom,NNHC)*dt_ys
!             vzbc31(iatom,NNHC) = vzbc31(iatom,NNHC)  &
!                                + 0.25d0*fzbc31(iatom,NNHC)*dt_ys
!
!!     /* update the paricle velocities */
!             vux(iatom,1) = vux(iatom,1)*scalex
!             vuy(iatom,1) = vuy(iatom,1)*scaley
!             vuz(iatom,1) = vuz(iatom,1)*scalez
!          enddo
!       enddo
!
!    Else
!
!!     /*  start multiple time step integration  */
!       do iys = 1, NYS
!
!!     /*  set time increment at this loop  */
!          dt_ys = dt*ysweight(iys)
!          do iatom = 1, NATOM
!             do icolor = 1, NCOLOR
!
!!     /*  calculate kinetic energy of centroid coordiates  */
!!     /*  for each atom                                    */
!
!                dkinx = fictmass(iatom,1)*vux(iatom,1)*vux(iatom,1)
!                dkiny = fictmass(iatom,1)*vuy(iatom,1)*vuy(iatom,1)
!                dkinz = fictmass(iatom,1)*vuz(iatom,1)*vuz(iatom,1)
!
!                scalex = 1.d0
!                scaley = 1.d0
!                scalez = 1.d0
!
!!     /* update the force */
!
!                fxbc3(iatom,1,icolor) = (dkinx - gkt)/qmcent3(1,icolor)
!                fybc3(iatom,1,icolor) = (dkiny - gkt)/qmcent3(1,icolor)
!                fzbc3(iatom,1,icolor) = (dkinz - gkt)/qmcent3(1,icolor)
!
!                do inhc = 2, NNHC
!                   fxbc3(iatom,inhc,icolor) = (qmcent3(inhc-1,icolor)    &
!                                            * vxbc3(iatom,inhc-1,icolor) &
!                                            * vxbc3(iatom,inhc-1,icolor) &
!                                            - gkt)/qmcent3(inhc,icolor)
!                   fybc3(iatom,inhc,icolor) = (qmcent3(inhc-1,icolor)    &
!                                            * vybc3(iatom,inhc-1,icolor) &
!                                            * vybc3(iatom,inhc-1,icolor) &
!                                            - gkt)/qmcent3(inhc,icolor)
!                   fzbc3(iatom,inhc,icolor) = (qmcent3(inhc-1,icolor)    &
!                                            * vzbc3(iatom,inhc-1,icolor) &
!                                            * vzbc3(iatom,inhc-1,icolor) &
!                                            - gkt)/qmcent3(inhc,icolor)
!                enddo
!
!!     /* update the thermostat velocities */
!
!                vxbc3(iatom,NNHC,icolor) = vxbc3(iatom,NNHC,icolor) &
!                                         + 0.25d0*fxbc3(iatom,NNHC,icolor)*dt_ys
!                vybc3(iatom,NNHC,icolor) = vybc3(iatom,NNHC,icolor) &
!                                         + 0.25d0*fybc3(iatom,NNHC,icolor)*dt_ys
!                vzbc3(iatom,NNHC,icolor) = vzbc3(iatom,NNHC,icolor) &
!                                         + 0.25d0*fzbc3(iatom,NNHC,icolor)*dt_ys
! 
!                do inhc = 1, NNHC-1
!                   vxfact=dexp(-0.125d0*vxbc3(iatom,NNHC-inhc+1,icolor)*dt_ys)
!                   vyfact=dexp(-0.125d0*vybc3(iatom,NNHC-inhc+1,icolor)*dt_ys)
!                   vzfact=dexp(-0.125d0*vzbc3(iatom,NNHC-inhc+1,icolor)*dt_ys)
!
!                   vxbc3(iatom,NNHC-inhc,icolor) = vxbc3(iatom,NNHC-inhc,icolor)        &
!                                                 * vxfact*vxfact                        &
!                                                 + 0.25d0*fxbc3(iatom,NNHC-inhc,icolor) &
!                                                 * vxfact*dt_ys   
!                   vybc3(iatom,NNHC-inhc,icolor) = vybc3(iatom,NNHC-inhc,icolor)        &
!                                                 * vyfact*vyfact                        &
!                                                 + 0.25d0*fybc3(iatom,NNHC-inhc,icolor) &
!                                                 * vyfact*dt_ys
!                   vzbc3(iatom,NNHC-inhc,icolor) = vzbc3(iatom,NNHC-inhc,icolor)        &
!                                                 * vzfact*vzfact                        &
!                                                 + 0.25d0*fzbc3(iatom,NNHC-inhc,icolor) &
!                                                 * vzfact*dt_ys
!                enddo
!
!!     /* update the particle velocities */
!
!                pvxfact = dexp(-0.5d0*vxbc3(iatom,1,icolor)*dt_ys)
!                pvyfact = dexp(-0.5d0*vybc3(iatom,1,icolor)*dt_ys)
!                pvzfact = dexp(-0.5d0*vzbc3(iatom,1,icolor)*dt_ys)
!        
!                scalex = scalex*pvxfact
!                scaley = scaley*pvyfact
!                scalez = scalez*pvzfact
!
!!     /* update the force */
!                fxbc3(iatom,1,icolor)=(scalex*scalex*dkinx - gkt)/qmcent3(1,icolor)
!                fybc3(iatom,1,icolor)=(scaley*scaley*dkiny - gkt)/qmcent3(1,icolor)
!                fzbc3(iatom,1,icolor)=(scalez*scalez*dkinz - gkt)/qmcent3(1,icolor)
!
!!     /* update the thermostat position */
!                do inhc = 1, NNHC
!                   xbc3(iatom,inhc,icolor) = xbc3(iatom,inhc,icolor)               &
!                                           + 0.5d0*vxbc3(iatom,inhc,icolor)*dt_ys
!                   ybc3(iatom,inhc,icolor) = ybc3(iatom,inhc,icolor)               &
!                                           + 0.5d0*vybc3(iatom,inhc,icolor)*dt_ys
!                   zbc3(iatom,inhc,icolor) = zbc3(iatom,inhc,icolor)               &
!                                           + 0.5d0*vzbc3(iatom,inhc,icolor)*dt_ys
!                enddo
!
!!     /* update the thermostat velocities */
!                do inhc = 1, NNHC-1
!                   vxfact = dexp(-0.125d0*vxbc3(iatom,inhc+1,icolor)*dt_ys)
!                   vyfact = dexp(-0.125d0*vybc3(iatom,inhc+1,icolor)*dt_ys)
!                   vzfact = dexp(-0.125d0*vzbc3(iatom,inhc+1,icolor)*dt_ys)
!
!                   vxbc3(iatom,inhc,icolor) = vxbc3(iatom,inhc,icolor)*vxfact*vxfact &
!                                            + 0.25d0*fxbc3(iatom,inhc,icolor)*vxfact*dt_ys
!                   vybc3(iatom,inhc,icolor) = vybc3(iatom,inhc,icolor)*vyfact*vyfact &
!                                            + 0.25d0*fybc3(iatom,inhc,icolor)*vyfact*dt_ys
!                   vzbc3(iatom,inhc,icolor) = vzbc3(iatom,inhc,icolor)*vzfact*vzfact &
!                                            + 0.25d0*fzbc3(iatom,inhc,icolor)*vzfact*dt_ys
!      
!                   fxbc3(iatom,inhc+1,icolor) = (qmcent3(inhc,icolor)      &
!                                              * vxbc3(iatom,inhc,icolor)   &
!                                              * vxbc3(iatom,inhc,icolor)   &
!                                              - gkt)/qmcent3(inhc+1,icolor) 
!                   fybc3(iatom,inhc+1,icolor) = (qmcent3(inhc,icolor)      &
!                                              * vybc3(iatom,inhc,icolor)   &
!                                              * vybc3(iatom,inhc,icolor)   &
!                                              - gkt)/qmcent3(inhc+1,icolor)
!                   fzbc3(iatom,inhc+1,icolor) = (qmcent3(inhc,icolor)      &
!                                              * vzbc3(iatom,inhc,icolor)   &
!                                              * vzbc3(iatom,inhc,icolor)   &
!                                              - gkt)/qmcent3(inhc+1,icolor) 
!                enddo
!
!                vxbc3(iatom,NNHC,icolor) = vxbc3(iatom,NNHC,icolor)  &
!                                         + 0.25d0*fxbc3(iatom,NNHC,icolor)*dt_ys
!                vybc3(iatom,NNHC,icolor) = vybc3(iatom,NNHC,icolor)  &
!                                         + 0.25d0*fybc3(iatom,NNHC,icolor)*dt_ys
!                vzbc3(iatom,NNHC,icolor) = vzbc3(iatom,NNHC,icolor)  &
!                                         + 0.25d0*fzbc3(iatom,NNHC,icolor)*dt_ys
!
!!     /* update the paricle velocities */
!                vux(iatom,1) = vux(iatom,1)*scalex
!                vuy(iatom,1) = vuy(iatom,1)*scaley
!                vuz(iatom,1) = vuz(iatom,1)*scalez
!             enddo
!          enddo
!       enddo
!
!    EndIf   
!
!    return
!    End Subroutine 
