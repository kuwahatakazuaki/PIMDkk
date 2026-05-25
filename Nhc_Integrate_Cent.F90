subroutine nhc_Integrate_Cent
  use Parameters
  use utility, only: norm_seq
  implicit none
  Double Precision :: skin, scale, dt_ys, vfact, pvfact
  integer :: i, iys, Iatom, Inhc

  skin = 0.d0
  do Iatom = 1, Natom
    skin = skin + fictmass(Iatom,1) * norm_seq( vur(:,Iatom,1) )
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
  vur(:,:,1) = vur(:,:,1) * scale

  Return
end subroutine nhc_Integrate_Cent


subroutine nhc_Integrate_Cent3
  Use Parameters
  implicit none
  real(8) :: skin, scale, dt_ys, vfact, pvfact
  real(8) :: vrfact(Ndim), pvrfact(Ndim), dkinr(Ndim)!, scaler(3)
  integer :: i, iys, Iatom, Inhc

!  /*  start multiple time step integration  */
  do iys = 1, Nys
!  /*  set time increment at this loop  */
    dt_ys = dt*ysweight(iys)
    do Iatom = 1, Natom

!  /*  calculate kinetic energy of centroid coordiates  */
!  /*  for each atom                                    */

      dkinr(:) = fictmass(Iatom,1)*vur(:,Iatom,1)*vur(:,Iatom,1)

!  /* update the force */

      frbc31(:,Iatom,1) = (dkinr(:) - gkt)/qmcent31(1)

      do Inhc = 2, Nnhc
        frbc31(:,Iatom,Inhc) &
          = (qmcent31(Inhc-1)*vrbc31(:,Iatom,Inhc-1)*vrbc31(:,Iatom,Inhc-1) - gkt)/qmcent31(Inhc)
      enddo

!  /* update the thermostat velocities */

      vrbc31(:,Iatom,Nnhc) = vrbc31(:,Iatom,Nnhc) + 0.25d0*frbc31(:,Iatom,Nnhc)*dt_ys

      do Inhc = 1, Nnhc-1
        vrfact = dexp( -0.125d0*vrbc31(:,Iatom,Nnhc-Inhc+1)*dt_ys )
        vrbc31(:,Iatom,Nnhc-Inhc) &
          = vrbc31(:,Iatom,Nnhc-Inhc)*vrfact*vrfact + 0.25d0*frbc31(:,Iatom,Nnhc-Inhc)*vrfact*dt_ys
      enddo

!  /* update the particle velocities */

      pvrfact(:) = dexp(-0.5d0*vrbc31(:,Iatom,1)*dt_ys)

!  /* update the force */
      frbc31(:,Iatom,1)=(pvrfact(:)*pvrfact(:)*dkinr(:) - gkt)/qmcent31(1)

!  /* update the thermostat position */
      do Inhc = 1, Nnhc
        rbc31(:,Iatom,Inhc) = rbc31(:,Iatom,Inhc) + 0.5d0*vrbc31(:,Iatom,Inhc)*dt_ys
      enddo

!  /* update the thermostat velocities */
      do Inhc = 1, Nnhc-1
        vrfact(:) = dexp(-0.125d0*vrbc31(:,Iatom,Inhc+1)*dt_ys)

        vrbc31(:,Iatom,Inhc) &
          = vrbc31(:,Iatom,Inhc)*vrfact(:)*vrfact(:) + 0.25d0*frbc31(:,Iatom,Inhc)*vrfact*dt_ys

        frbc31(:,Iatom,Inhc+1) &
          = (qmcent31(Inhc)*vrbc31(:,Iatom,Inhc)*vrbc31(:,Iatom,Inhc) - gkt)/qmcent31(Inhc+1) 
      enddo

      vrbc31(:,Iatom,Nnhc) = vrbc31(:,Iatom,Nnhc) + 0.25d0*frbc31(:,Iatom,Nnhc)*dt_ys

!  /* update the paricle velocities */
      vur(:,Iatom,1) = vur(:,Iatom,1)*pvrfact(:)
    enddo
  enddo

return
end subroutine nhc_Integrate_Cent3
