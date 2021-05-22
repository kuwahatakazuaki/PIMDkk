  Subroutine Remove_TnR_All

    Use Parameters
    Implicit None

    Double Precision :: sumvx,sumvy,sumvz,totmas
    Double Precision :: comx,comy,comz
    Double Precision :: avelx,avely,avelz
    Double Precision :: cixx,ciyy,cizz,cixy,ciyz,cizx
    Double Precision :: detmat,gradx,grady,gradz
    Double Precision, Dimension(1:3,1:3) :: dinvmat
!
!   Remove Translation and Rotation Velocity of 
!   COM for Centroids and Noncentroids
!
    do imode=1,nbead
       comx  = 0.d0
       comy  = 0.d0
       comz  = 0.d0
       sumvx = 0.d0
       sumvy = 0.d0
       sumvz = 0.d0
       totmas= 0.0d0
       cixx  = 0.0d0
       ciyy  = 0.0d0
       cizz  = 0.0d0
       cixy  = 0.0d0
       ciyz  = 0.0d0
       cizx  = 0.0d0
       avelx = 0.0d0
       avely = 0.0d0
       avelz = 0.0d0
!
! Calculate Center of Mass and Translational Velocity of COM
! 
       do iatom = 1, natom
          sumvx  = sumvx  + vux(iatom,imode)*fictmass(iatom,imode)
          sumvy  = sumvy  + vuy(iatom,imode)*fictmass(iatom,imode)
          sumvz  = sumvz  + vuz(iatom,imode)*fictmass(iatom,imode)
          comx   = comx   + ux(iatom,imode)*fictmass(iatom,imode)
          comy   = comy   + uy(iatom,imode)*fictmass(iatom,imode)
          comx   = comz   + uz(iatom,imode)*fictmass(iatom,imode)
          totmas = totmas + fictmass(iatom,imode)
       enddo
       sumvx = sumvx / totmas
       sumvy = sumvy / totmas
       sumvz = sumvz / totmas
       comx  = comx  / totmas
       comy  = comy  / totmas
       comz  = comz  / totmas
!
! Calculate Rotational Velocity of COM
! Construct the Momenta 
!
       do iatom = 1, natom
          cixx   = cixx   + fictmass(iatom,imode) &
                          * ((uy(iatom,imode)-comy)**2+(uz(iatom,imode)-comz)**2)
          ciyy   = ciyy   + fictmass(iatom,imode) &
                          * ((ux(iatom,imode)-comx)**2+(uz(iatom,imode)-comz)**2)
          cizz   = cizz   + fictmass(iatom,imode) &
                          * ((uy(iatom,imode)-comy)**2+(ux(iatom,imode)-comx)**2)
          cixy   = cixy   - fictmass(iatom,imode) &
                          * (ux(iatom,imode)-comx)*(uy(iatom,imode)-comy)
          ciyz   = ciyz   - fictmass(iatom,imode) &
                          * (uz(iatom,imode)-comz)*(uy(iatom,imode)-comy)
          cizx   = cizx   - fictmass(iatom,imode) &
                          * (ux(iatom,imode)-comx)*(uz(iatom,imode)-comz)
       enddo
!
! Obtain the Inverse Matrix of Momenta
!
       detmat= cixx*ciyy*cizz + 2.0d0*cixy*ciyz*cizx &
              -ciyz*ciyz*cixx - cizx*cizx*ciyy - cixy*cixy*cizz
       if(detmat<=0.00001d0) goto 100
       detmat=1.0d0/detmat
       dinvmat(1,1)=(ciyy*cizz-ciyz*ciyz)*detmat
       dinvmat(1,2)=(ciyz*cizx-cizz*cixy)*detmat
       dinvmat(1,3)=(cixy*ciyz-cizx*ciyy)*detmat
       dinvmat(2,1)=dinvmat(1,2)
       dinvmat(2,2)=(cizz*cixx-cizx*cizx)*detmat
       dinvmat(2,3)=(cizx*cixy-cixx*ciyz)*detmat
       dinvmat(3,1)=dinvmat(1,3)
       dinvmat(3,2)=dinvmat(2,3)
       dinvmat(3,3)=(cixx*ciyy-cixy*cixy)*detmat
!
!  Calculate Angular Velocities from Angular Momentum
!
       do iatom = 1, natom
          gradx = (uy(iatom,imode)-comy)*vuz(iatom,imode) &
                - (uz(iatom,imode)-comz)*vuy(iatom,imode)
          grady = (uz(iatom,imode)-comz)*vux(iatom,imode) &
                - (ux(iatom,imode)-comx)*vuz(iatom,imode)
          gradz = (ux(iatom,imode)-comx)*vuy(iatom,imode) &
                - (uy(iatom,imode)-comy)*vux(iatom,imode)
          avelx = avelx + fictmass(iatom,imode) &
                * (dinvmat(1,1)*gradx + dinvmat(1,2)*grady + dinvmat(1,3)*gradz)
          avely = avely + fictmass(iatom,imode) &
                * (dinvmat(2,1)*gradx + dinvmat(2,2)*grady + dinvmat(2,3)*gradz)
          avelz = avelz + fictmass(iatom,imode) &
                * (dinvmat(3,1)*gradx + dinvmat(3,2)*grady + dinvmat(3,3)*gradz)
       enddo

 100   continue
!
! Subtract Translation and Rotation Velocities of COM from 
! Velocities of Atoms
!
       do iatom = 1, natom
          vux(iatom,imode) = vux(iatom,imode) - sumvx &
                       - avely*(uz(iatom,imode)-comz) &
                       + avelz*(uy(iatom,imode)-comy)
          vuy(iatom,imode) = vuy(iatom,imode) - sumvy &
                       - avelz*(ux(iatom,imode)-comx) &
                       + avelx*(uz(iatom,imode)-comz)
          vuz(iatom,imode) = vuz(iatom,imode) - sumvz &
                       - avelx*(uy(iatom,imode)-comy) &
                       + avely*(ux(iatom,imode)-comx)
       enddo


    enddo

    Return
  End Subroutine
