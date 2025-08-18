subroutine Init_Bath
  use Parameters
  use utility, only: gasdev
  Implicit None
  real(8) :: vsigma
  integer :: inhc, imode, iatom, icolor

!     /*  massive Nose-Hoover chain  */
!YK Remove the initial velocity for the centroid 
  rbath(:, :,:,1) = 0.d0
  vrbath(:,:,:,1) = 0.d0

  do imode = 2, nbead
    vsigma = dsqrt(1.d0/beta/qmass(imode))
    do inhc = 1, nnhc
      do iatom = 1, natom
        !call gasdev(gasd); vrbath(1,iatom,inhc,imode) = vsigma*gasd
        !call gasdev(gasd); vrbath(2,iatom,inhc,imode) = vsigma*gasd
        !call gasdev(gasd); vrbath(3,iatom,inhc,imode) = vsigma*gasd
        vrbath(1,iatom,inhc,imode) = vsigma*gasdev()
        vrbath(2,iatom,inhc,imode) = vsigma*gasdev()
        vrbath(3,iatom,inhc,imode) = vsigma*gasdev()

        rbath(:,iatom,inhc,imode) = 0.d0
      enddo
    enddo
  enddo

!YK how about one NH chain?
  select case(Ncent)
    case(1)
      do inhc = 1, nnhc
        vsigma = dsqrt(1.d0/beta/qmcent11(inhc))
        !call gasdev(gasd)
        vbc11(inhc) = vsigma*gasdev()
        rbc11(inhc) = 0.d0
      enddo
    case(3)
      do inhc = 1, nnhc
        do iatom=1,natom
          vsigma = dsqrt(1.d0/beta/qmcent31(inhc))
          !call gasdev(gasd); vrbc31(1,iatom,inhc) = vsigma*gasd
          !call gasdev(gasd); vrbc31(2,iatom,inhc) = vsigma*gasd
          !call gasdev(gasd); vrbc31(3,iatom,inhc) = vsigma*gasd
          vrbc31(1,iatom,inhc) = vsigma*gasdev()
          vrbc31(2,iatom,inhc) = vsigma*gasdev()
          vrbc31(3,iatom,inhc) = vsigma*gasdev()
          rbc31(:,iatom,inhc) = 0.d0
        enddo
      enddo
  end select

  return
end subroutine
