Subroutine Init_Bath
  use Parameters
  Implicit None
  Double Precision :: gasd, vsigma
  integer :: inhc, imode, iatom, icolor

!     /*  massive Nose-Hoover chain  */
!YK Remove the initial velocity for the centroid 
! which is not used
  do inhc = 1, nnhc
    do iatom= 1, natom
      rbath(:,iatom,inhc,1) = 0.d0
      vrbath(:,iatom,inhc,1) = 0.d0
    enddo
  enddo

  do imode = 2, nbead
    vsigma = dsqrt(1.d0/beta/qmass(imode))
    do inhc = 1, nnhc
      do iatom = 1, natom
        call gasdev(gasd); vrbath(1,iatom,inhc,imode) = vsigma*gasd
        call gasdev(gasd); vrbath(2,iatom,inhc,imode) = vsigma*gasd
        call gasdev(gasd); vrbath(3,iatom,inhc,imode) = vsigma*gasd

        rbath(:,iatom,inhc,imode) = 0.d0
      enddo
    enddo
  enddo
!
!   /*  single thermostat attached to centroids  */

!YK Commented out original since qmass_cent differs when i=1 or not
!     vsigma = dsqrt(1.d0/beta/qmass(1))
!     do i = 1, nnhc
!       vbathcent(i) = vsigma*gasdev()
!       rbathcent(i) = 0.d0
!     enddo
!YK
!     /*  single thermostat attached to centroids  */
!YK how about one NH chain?
  select case(Ncent)
    case(1)
      do inhc = 1, nnhc
        vsigma = dsqrt(1.d0/beta/qmcent11(inhc))
        call gasdev(gasd)
        vbc11(inhc) = vsigma*gasd
        rbc11(inhc) = 0.d0
      enddo
    case(3)
      do inhc = 1, nnhc
        do iatom=1,natom
          vsigma = dsqrt(1.d0/beta/qmcent31(inhc))
          call gasdev(gasd); vrbc31(1,iatom,inhc) = vsigma*gasd
          call gasdev(gasd); vrbc31(2,iatom,inhc) = vsigma*gasd
          call gasdev(gasd); vrbc31(3,iatom,inhc) = vsigma*gasd
          rbc31(:,iatom,inhc) = 0.d0
        enddo
      enddo
  end select

  Return
End Subroutine
