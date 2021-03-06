subroutine integ_nhc_cent
  use global_variable, only: Natom, Ncent, Ncolor, Nnhc, Nys, &
      r, vu, fictmass, ysweight, &
      rbc11, vbc11, fbc11, qmcent11, &
      rbc1,  vbc1,  fbc1,  qmcent1, &
      rbc31, vbc31, fbc31, qmcent31, &
      rbc3,  vbc3,  fbc3,  qmcent3, &
      gnkt, gkt
  use utility
  implicit none
  integer :: i, inhc, iys, icolor

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
    real(8) :: skin, vfact, pvfact
    real(8) :: dt, dt_ys

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
        dt_ys = dt*ysweight(iys)

        vbc11(Nnhc) = vbc11(Nnhc) + 0.25d0 * fbc11(Nnhc) * dt_ys
        do inhc = 1, Nnhc-1
          vfact = dexp( -0.125d0 * vbc11(Nnhc-inhc+1) * dt_ys )
          vbc11(Nnhc-inhc) = vbc11(Nnhc-inhc)*vfact*vfact + 0.25d0*fbc11(Nnhc-inhc)*vfact*dt_ys
        end do

        pvfact = dexp( -0.5d0*vbc11(1)*dt_ys )
        fbc11(1) = ( pvfact * pvfact * skin - gnkt ) / qmcent11(1)

        rbc11(:) = rbc11(:) + 0.5d0 * vbc11(:) * dt_ys

        do inhc = 1, Nnhc-1
          vfact = dexp( -0.125d0 * vbc11(inhc+1) * dt_ys )
          vbc11(inhc) = vbc11(inhc) * vfact*vfact + 0.25d0 * fbc11(inhc) * vfact * dt_ys
          fbc11(inhc+1) = (qmcent11(inhc)*vbc11(inhc)*vbc11(inhc) - gkt ) / qmcent11(inhc+1)
        end do
        vbc11(Nnhc) = vbc11(Nnhc) + 0.25d0 * fbc11(Nnhc) * dt_ys
      end do
!      qmcent11(1) = 3.0d0 * dble(Natom) / beta / omega2
!      qmcent11(2:Nnhc) = 1.0d0 / beta / omega2

    else if ( Ncolor > 1 ) then
      do icolor = 1, Ncolor
        fbc1(1,icolor) = (skin - gnkt)/qmcent1(1,icolor)
        do inhc = 2, Nnhc
          fbc1(inhc,icolor) = ( qmcent1(inhc-1,icolor) * vbc1(inhc-1,icolor) * vbc1(inhc-1,icolor) - gkt ) / qmcent1(inhc,icolor)
        end do
      end do

      do icolor = 1, Ncolor
        do iys = 1, Nys
          dt_ys = dt*ysweight(iys)
          vbc1(Nnhc,icolor) = vbc1(Nnhc,icolor) + 0.25d0 * fbc1(Nnhc,icolor) * dt_ys

          do inhc = 1, Nnhc-1
            vfact = dexp( -0.125d0 * vbc1(Nnhc-inhc+1,icolor) * dt_ys )
            vbc1(Nnhc-inhc,icolor) = vbc1(Nnhc-inhc,icolor) * vfact * vfact + 0.25d0 * fbc1(Nnhc-inhc,icolor) * vfact * dt_ys
          end do

          pvfact = dexp(-0.5d0 * vbc1(1,icolor) * dt_ys )
          fbc1(1,icolor) = (pvfact*pvfact*skin - gnkt) / qmcent1(1,icolor)

          rbc1(:,icolor) = rbc1(:,icolor) + 0.5d0 * vbc1(:,icolor) * dt_ys

          do inhc = 1, Nnhc-1
            vfact = dexp( -0.125d0 * vbc1(inhc+1,icolor) * dt_ys )
            vbc1(inhc,icolor) = vbc1(inhc,icolor) * vfact * vfact + 0.25d0 * fbc1(inhc,icolor) * vfact * dt_ys

            fbc1(inhc+1,icolor) = (qmcent1(inhc,icolor) * vbc1(inhc,icolor) * vbc1(inhc,icolor) - gkt ) / qmcent1(inhc+1,icolor)
          end do
          vbc1(Nnhc,icolor) = vbc1(Nnhc,icolor) + 0.25d0 + fbc1(Nnhc,icolor) * dt_ys
        end do
      end do
    end if

    vu(:,:,1) = vu(:,:,1) * pvfact

  end subroutine integ_nhc_cent1
! *** End Ncent = 1 ***

! *** Ncent = 3 ***
  subroutine integ_nhc_cent3
    real(8) :: dt, dt_ys, vfact(3), pvfact(3)

    real(8) :: dkin(3)

    if (Ncolor == 1) then
      do iys = 1, Nys
        dt_ys = dt*ysweight(iys)

        do i = 1, Natom
          dkin(:) = fictmass(i,1) * vu(:,i,1) * vu(:,i,1)
          fbc31(:,i,1) = (dkin(:) - gkt) / qmcent31(1)

          do inhc = 2, Nnhc
            fbc31(:,i,inhc) = (qmcent31(inhc-1) * vbc31(:,i,inhc-1) * vbc31(:,i,inhc-1) - gkt) / qmcent31(inhc)
          end do

          vbc31(:,i,Nnhc) = vbc31(:,i,inhc) + 0.25d0 * fbc31(:,i,Nnhc) * dt_ys

          do inhc = 1, Nnhc-1
            vfact(:) = dexp( -0.125d0 * vbc31(:,i,Nnhc-inhc+1) * dt_ys )
            vbc31(:,i,Nnhc-inhc) = vbc31(:,i,Nnhc-inhc) + 0.25d0 * fbc31(:,i,Nnhc-inhc) * vfact(:) * dt_ys
          end do

          pvfact(:) = dexp( -0.5d0 * vbc31(:,i,1) * dt_ys )
          fbc31(:,i,1) = ( pvfact(:) * pvfact(:) * dkin(:) - gkt ) / qmcent31(1)

          !rbc31(:,i,inhc) = rbc31(:,i,inhc) + 0.5d0 * vbc31(:,i,inhc) * dt_ys
          do inhc = 1, Nnhc
            rbc31(:,i,inhc) = rbc31(:,i,inhc) + 0.5d0 * vbc31(:,i,inhc) * dt_ys
          end do

          do inhc = 1, Nnhc-1
            vfact(:) = dexp( -0.125d0 * vbc31(:,i,inhc+1) * dt_ys)
            vbc31(:,i,inhc) = vbc31(:,i,inhc) * vfact * vfact + 0.25d0 * fbc31(:,i,inhc) * vfact * dt_ys
            fbc31(:,i,inhc+1) = ( qmcent31(inhc) * vbc31(:,i,inhc) * vbc31(:,i,inhc) - gkt) / qmcent31(inhc+1)
          end do

          vbc31(:,i,Nnhc) = vbc31(:,i,Nnhc) + 0.25d0 * fbc31(:,i,Nnhc) * dt_ys
          vu(:,i,1) = vu(:,i,1) * pvfact(:)

        end do

      end do

    else if (Ncolor > 1) then
      do iys = 1, Nys
        dt_ys = dt*ysweight(iys)
        do i = 1, Natom
          do icolor = 1, Ncolor
            dkin(:) = fictmass(i,1) * vu(:,i,1) * vu(:,i,1)
            fbc3(:,i,1,icolor) = (dkin(:) - gkt) / qmcent3(1,icolor)

            do inhc = 2, Nnhc
              ! +++ HERE +++ !
            end do

          end do
        end do
      end do
          !do inhc = 1, Nnhc-1
          !  vfact(:) = dexp( -0.125d0 * vbc3(:,i,inhc+1,icolor) * dt_ys )
          !  vbc3(:,i,inhc,icolor) = vbc3(:,i,inhc,icolor) * vfact(:) * vfact(:) &
          !                          + 0.25d0 * fbc3(:,i,inhc,icolor) * vfact(:) * dt_ys
          !  fbc3(:,i,inhc+1,icolor) = &
          !      (qmcent3(inhc,icolor) * vbc3(:,i,inhc,icolor) * vbc3(:,i,inhc,icolor) - gkt) / qmcent3(inhc+1,icolor) 
          !end do
    end if
  end subroutine integ_nhc_cent3
! *** End Ncent = 3 ***

end subroutine integ_nhc_cent

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
