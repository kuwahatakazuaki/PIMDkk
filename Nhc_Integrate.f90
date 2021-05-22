  Subroutine Nhc_Integrate

    Use Parameters
    Implicit None

    Integer                            :: iys, i
    Double Precision                   :: dkinx,  dkiny, dkinz 
    Double Precision                   :: scalex, scaley, scalez 
    Double Precision                   :: pvxfact, pvyfact, pvzfact 
    Double Precision                   ::  vxfact,  vyfact,  vzfact 
    Double Precision                   ::  dt_ys 

    do iys = 1, nys

       dt_ys = dt_ref*ysweight(iys)

!YK Excluded imode=1, which should do nothing
       do imode = 2, nbead
!YK
          do iatom = 1, natom

!     /*  calculate kinetic energy of normal mode coordinate  *
!      *  for each atom                                       */
             dkinx = fictmass(iatom,imode)*vux(iatom,imode)*vux(iatom,imode)
             dkiny = fictmass(iatom,imode)*vuy(iatom,imode)*vuy(iatom,imode)
             dkinz = fictmass(iatom,imode)*vuz(iatom,imode)*vuz(iatom,imode)
             scalex = 1.d0
             scaley = 1.d0
             scalez = 1.d0

!     /* update the force */
             fxbath(iatom,1,imode) = (dkinx - gkt)/qmass(imode)
             fybath(iatom,1,imode) = (dkiny - gkt)/qmass(imode)
             fzbath(iatom,1,imode) = (dkinz - gkt)/qmass(imode)
             do inhc = 2, nnhc
                fxbath(iatom,inhc,imode) = (qmass(imode) &
                                         * vxbath(iatom,inhc-1,imode) &
                                         * vxbath(iatom,inhc-1,imode) &
                                         - gkt)/qmass(imode)
                fybath(iatom,inhc,imode) = (qmass(imode) &
                                         * vybath(iatom,inhc-1,imode) &
                                         * vybath(iatom,inhc-1,imode) &
                                         - gkt)/qmass(imode)
                fzbath(iatom,inhc,imode) = (qmass(imode) &
                                         * vzbath(iatom,inhc-1,imode) &
                                         * vzbath(iatom,inhc-1,imode) &
                                         - gkt)/qmass(imode)
             enddo

!     /* update the thermostat velocities */
             vxbath(iatom,nnhc,imode) = vxbath(iatom,nnhc,imode) &
                                      + 0.25d0*fxbath(iatom,nnhc,imode)*dt_ys
             vybath(iatom,nnhc,imode) = vybath(iatom,nnhc,imode) &
                                      + 0.25d0*fybath(iatom,nnhc,imode)*dt_ys
             vzbath(iatom,nnhc,imode) = vzbath(iatom,nnhc,imode) &
                                      + 0.25d0*fzbath(iatom,nnhc,imode)*dt_ys

             do inhc = 1, nnhc-1

                vxfact=dexp(-0.125d0*vxbath(iatom,nnhc-inhc+1,imode)*dt_ys)
                vyfact=dexp(-0.125d0*vybath(iatom,nnhc-inhc+1,imode)*dt_ys)
                vzfact=dexp(-0.125d0*vzbath(iatom,nnhc-inhc+1,imode)*dt_ys)
              
                vxbath(iatom,nnhc-inhc,imode) = vxbath(iatom,nnhc-inhc,imode)*vxfact*vxfact &
                                              + 0.25d0*fxbath(iatom,nnhc-inhc,imode)*vxfact*dt_ys
                vybath(iatom,nnhc-inhc,imode) = vybath(iatom,nnhc-inhc,imode)*vyfact*vyfact &
                                              + 0.25d0*fybath(iatom,nnhc-inhc,imode)*vyfact*dt_ys
                vzbath(iatom,nnhc-inhc,imode) = vzbath(iatom,nnhc-inhc,imode)*vzfact*vzfact &
                                              + 0.25d0*fzbath(iatom,nnhc-inhc,imode)*vzfact*dt_ys
             enddo

!     /* update the particle velocities */
             pvxfact = dexp(-0.5d0*vxbath(iatom,1,imode)*dt_ys)
             pvyfact = dexp(-0.5d0*vybath(iatom,1,imode)*dt_ys)
             pvzfact = dexp(-0.5d0*vzbath(iatom,1,imode)*dt_ys)
     
             scalex = scalex*pvxfact
             scaley = scaley*pvyfact
             scalez = scalez*pvzfact

!     /* update the force */
             fxbath(iatom,1,imode)=(scalex*scalex*dkinx - gkt)/qmass(imode)
             fybath(iatom,1,imode)=(scaley*scaley*dkiny - gkt)/qmass(imode)
             fzbath(iatom,1,imode)=(scalez*scalez*dkinz - gkt)/qmass(imode)

!     /* update the thermostat position */
             do inhc = 1, nnhc
                xbath(iatom,inhc,imode) = xbath(iatom,inhc,imode) &
                                        + 0.5d0*vxbath(iatom,inhc,imode)*dt_ys
                ybath(iatom,inhc,imode) = ybath(iatom,inhc,imode) &
                                        + 0.5d0*vybath(iatom,inhc,imode)*dt_ys
                zbath(iatom,inhc,imode) = zbath(iatom,inhc,imode) &
                                        + 0.5d0*vzbath(iatom,inhc,imode)*dt_ys
             enddo

!     /* update the tehrmostat velocities */
             do inhc = 1, nnhc-1

                vxfact = dexp(-0.125d0*vxbath(iatom,inhc+1,imode)*dt_ys)
                vyfact = dexp(-0.125d0*vybath(iatom,inhc+1,imode)*dt_ys)
                vzfact = dexp(-0.125d0*vzbath(iatom,inhc+1,imode)*dt_ys)
    
                vxbath(iatom,inhc,imode) = vxbath(iatom,inhc,imode)*vxfact*vxfact &
                                         + 0.25d0*fxbath(iatom,inhc,imode)*vxfact*dt_ys
                vybath(iatom,inhc,imode) = vybath(iatom,inhc,imode)*vyfact*vyfact &
                                         + 0.25d0*fybath(iatom,inhc,imode)*vyfact*dt_ys
                vzbath(iatom,inhc,imode) = vzbath(iatom,inhc,imode)*vzfact*vzfact &
                                         + 0.25d0*fzbath(iatom,inhc,imode)*vzfact*dt_ys

                fxbath(iatom,inhc+1,imode) = (qmass(imode)*vxbath(iatom,inhc,imode) &
                                           * vxbath(iatom,inhc,imode)- gkt)/qmass(imode)
                fybath(iatom,inhc+1,imode) = (qmass(imode)*vybath(iatom,inhc,imode) &
                                           * vybath(iatom,inhc,imode)- gkt)/qmass(imode)
                fzbath(iatom,inhc+1,imode) = (qmass(imode)*vzbath(iatom,inhc,imode) &
                                           * vzbath(iatom,inhc,imode)- gkt)/qmass(imode)
             enddo

             vxbath(iatom,nnhc,imode) = vxbath(iatom,nnhc,imode) + 0.25d0*fxbath(iatom,nnhc,imode)*dt_ys
             vybath(iatom,nnhc,imode) = vybath(iatom,nnhc,imode) + 0.25d0*fybath(iatom,nnhc,imode)*dt_ys
             vzbath(iatom,nnhc,imode) = vzbath(iatom,nnhc,imode) + 0.25d0*fzbath(iatom,nnhc,imode)*dt_ys

!     /* update the paricle velocities */
             vux(iatom,imode) = vux(iatom,imode)*scalex
             vuy(iatom,imode) = vuy(iatom,imode)*scaley
             vuz(iatom,imode) = vuz(iatom,imode)*scalez
          enddo
       enddo
    enddo

    Return

  End Subroutine
