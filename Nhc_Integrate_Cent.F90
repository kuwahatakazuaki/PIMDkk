Subroutine Nhc_Integrate_Cent
  use Parameters
  use utility, only: norm_seq
  implicit none
  Double Precision :: skin, scale, dt_ys, vfact, pvfact
  integer :: i, iys, iatom, Inhc, icolor

  skin = 0.d0
  do iatom = 1, Natom
    skin = skin + fictmass(iatom,1) * norm_seq( vur(:,iatom,1) )
  enddo

  scale = 1.d0

!   /* update the force */

  fbc11(1) = (skin - gnkt)/qmcent11(1)
  do Inhc = 2, nnhc
     fbc11(Inhc) = (qmcent11(Inhc-1)*vbc11(Inhc-1)*vbc11(Inhc-1) &
                 - gkt)/qmcent11(Inhc)
  enddo

!  /*  start multiple time step integration  */
  do iys = 1, nys

!  /*  set time increment at this loop  */
     dt_ys = dt*ysweight(iys)

!  /* update the thermostat velocities */
     vbc11(nnhc) = vbc11(nnhc) + 0.25d0*fbc11(nnhc)*dt_ys
     do Inhc = 1, nnhc-1
        vfact=dexp(-0.125d0*   vbc11(nnhc-Inhc+1)*dt_ys)
        vbc11(nnhc-Inhc) = vbc11(nnhc-Inhc)*vfact*vfact &
                         + 0.25d0*fbc11(nnhc-Inhc)*vfact*dt_ys
     enddo

!  /* update the particle velocities */
     pvfact = dexp(-0.5d0*vbc11(1)*dt_ys)
     scale = scale*pvfact
!  /* update the force */
     fbc11(1)=(scale*scale*skin - gnkt)/qmcent11(1)
!  /* update the thermostat position */
     do Inhc = 1, nnhc
        rbc11(Inhc) = rbc11(Inhc) &
                    + 0.5d0*vbc11(Inhc)*dt_ys
     enddo

!  /* update the thermostat velocities */
     do Inhc = 1, nnhc-1
        vfact = dexp(-0.125d0*vbc11(Inhc+1)*dt_ys)
        vbc11(Inhc) = vbc11(Inhc)*vfact*vfact &
                    + 0.25d0*fbc11(Inhc)*vfact*dt_ys
        fbc11(Inhc+1) = (qmcent11(Inhc)*vbc11(Inhc)*vbc11(Inhc) &
                      - gkt)/qmcent11(Inhc+1)
     enddo

     vbc11(nnhc) = vbc11(nnhc) + 0.25d0*fbc11(nnhc)*dt_ys
  enddo


!     /* update the paricle velocities */
  !do iatom = 1, natom
  !  vur(:,iatom,1) = vur(:,iatom,1)*scale
  !enddo
  vur(:,:,1) = vur(:,:,1) * scale

  Return
End Subroutine


Subroutine nhc_integrate_cent3
  Use Parameters
  Implicit None
  Double Precision  :: skin, scale, dt_ys, vfact, pvfact
  real(8) :: vrfact(3), pvrfact(3), scaler(3), dkinr(3)
  integer :: i, iys, iatom, Inhc, icolor

!  /*  start multiple time step integration  */
  do iys = 1, Nys
!  /*  set time increment at this loop  */
    dt_ys = dt*ysweight(iys)
    do iatom = 1, Natom

!  /*  calculate kinetic energy of centroid coordiates  */
!  /*  for each atom                                    */

      dkinr(:) = fictmass(iatom,1)*vur(:,iatom,1)*vur(:,iatom,1)
      scaler(:) = 1.d0

!  /* update the force */

      frbc31(:,iatom,1) = (dkinr(:) - gkt)/qmcent31(1)

      do Inhc = 2, Nnhc
        frbc31(:,iatom,Inhc) &
          = (qmcent31(Inhc-1)*vrbc31(:,iatom,Inhc-1)*vrbc31(:,iatom,Inhc-1) - gkt)/qmcent31(Inhc)
      enddo

!  /* update the thermostat velocities */

      vrbc31(:,iatom,Nnhc) = vrbc31(:,iatom,Nnhc) + 0.25d0*frbc31(:,iatom,Nnhc)*dt_ys

      do Inhc = 1, Nnhc-1
        vrfact = dexp( -0.125d0*vrbc31(:,iatom,Nnhc-Inhc+1)*dt_ys )
        vrbc31(:,iatom,Nnhc-Inhc) &
          = vrbc31(:,iatom,Nnhc-Inhc)*vrfact*vrfact + 0.25d0*frbc31(:,iatom,Nnhc-Inhc)*vrfact*dt_ys
      enddo

!  /* update the particle velocities */

      pvrfact(:) = dexp(-0.5d0*vrbc31(:,iatom,1)*dt_ys)

      scaler(:) = scaler(:)*pvrfact

!  /* update the force */
      frbc31(:,iatom,1)=(scaler(:)*scaler(:)*dkinr(:) - gkt)/qmcent31(1)

!  /* update the thermostat position */
      do Inhc = 1, Nnhc
        rbc31(:,iatom,Inhc) = rbc31(:,iatom,Inhc) + 0.5d0*vrbc31(:,iatom,Inhc)*dt_ys
      enddo

!  /* update the thermostat velocities */
      do Inhc = 1, Nnhc-1
        vrfact(:) = dexp(-0.125d0*vrbc31(:,iatom,Inhc+1)*dt_ys)

        vrbc31(:,iatom,Inhc) &
          = vrbc31(:,iatom,Inhc)*vrfact(:)*vrfact(:) + 0.25d0*frbc31(:,iatom,Inhc)*vrfact*dt_ys

        frbc31(:,iatom,Inhc+1) &
          = (qmcent31(Inhc)*vrbc31(:,iatom,Inhc)*vrbc31(:,iatom,Inhc) - gkt)/qmcent31(Inhc+1) 
      enddo

      vrbc31(:,iatom,Nnhc) = vrbc31(:,iatom,Nnhc) + 0.25d0*frbc31(:,iatom,Nnhc)*dt_ys

!  /* update the paricle velocities */
      vur(:,iatom,1) = vur(:,iatom,1)*scaler(:)
    enddo
  enddo

return
End Subroutine
