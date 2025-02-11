Subroutine Nhc_Integrate
  Use Parameters
  use utility, only: norm_seq
  Implicit None

  Integer           :: iys, i
  real(8) :: dkinr(3), pvrfact(3), vrfact(3)!, scaler(3)
  Double Precision  :: dt_ys
  integer :: iatom, imode, inhc

  do iys = 1, nys

    dt_ys = dt_ref*ysweight(iys)

!YK Excluded imode=1, which should do nothing
    do imode = 2, nbead

      do iatom = 1, natom

!   /*  calculate kinetic energy of normal mode coordinate  *
!    *  for each atom                                       */
        dkinr(:) = fictmass(iatom,imode)*vur(:,iatom,imode)*vur(:,iatom,imode)

!   /* update the force */
        frbath(:,iatom,1,imode) = (dkinr(:) - gkt)/qmass(imode)
        do inhc = 2, nnhc
          frbath(:,iatom,inhc,imode) &
            = (qmass(imode)*vrbath(:,iatom,inhc-1,imode)*vrbath(:,iatom,inhc-1,imode) - gkt)/qmass(imode)
        enddo

!   /* update the thermostat velocities */
        vrbath(:,iatom,nnhc,imode) = vrbath(:,iatom,nnhc,imode) + 0.25d0*frbath(:,iatom,nnhc,imode)*dt_ys

        do inhc = 1, nnhc-1
          vrfact(:)=dexp(-0.125d0*vrbath(:,iatom,nnhc-inhc+1,imode)*dt_ys)

          vrbath(:,iatom,nnhc-inhc,imode) &
            = vrbath(:,iatom,nnhc-inhc,imode)*vrfact(:)*vrfact(:) &
              + 0.25d0*frbath(:,iatom,nnhc-inhc,imode)*vrfact(:)*dt_ys
        enddo

!   /* update the particle velocities */
        pvrfact(:) = dexp(-0.5d0*vrbath(:,iatom,1,imode)*dt_ys)

!   /* update the force */
        frbath(:,iatom,1,imode)=(pvrfact(:)*pvrfact(:)*dkinr(:) - gkt)/qmass(imode)

!   /* update the thermostat position */
        do inhc = 1, nnhc
          rbath(:,iatom,inhc,imode) &
            = rbath(:,iatom,inhc,imode) + 0.5d0*vrbath(:,iatom,inhc,imode)*dt_ys
        enddo

!   /* update the tehrmostat velocities */
        do inhc = 1, nnhc-1
          vrfact(:) = dexp(-0.125d0*vrbath(:,iatom,inhc+1,imode)*dt_ys)

          vrbath(:,iatom,inhc,imode) &
            = vrbath(:,iatom,inhc,imode)*vrfact(:)*vrfact(:) &
              + 0.25d0*frbath(:,iatom,inhc,imode)*vrfact*dt_ys

          frbath(:,iatom,inhc+1,imode) &
            = (qmass(imode)*vrbath(:,iatom,inhc,imode)*vrbath(:,iatom,inhc,imode) - gkt) &
              / qmass(imode)
        enddo

        vrbath(:,iatom,nnhc,imode) &
          = vrbath(:,iatom,nnhc,imode) + 0.25d0*frbath(:,iatom,nnhc,imode)*dt_ys

!   /* update the paricle velocities */
        vur(:,iatom,imode) = vur(:,iatom,imode)*pvrfact(:)
      enddo
    enddo
  enddo

  Return
End Subroutine
