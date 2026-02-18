subroutine Init_Bath
  use Parameters
  use utility, only: gasdev
  Implicit None
  real(8) :: vsigma
  integer :: inhc, imode, iatom, icolor, Idim_

!     /*  massive Nose-Hoover chain  */
!YK Remove the initial velocity for the centroid 
  rbath(:, :,:,1) = 0.d0
  vrbath(:,:,:,1) = 0.d0

  do imode = 2, Nbead
    vsigma = dsqrt(1.d0/beta/qmass(imode))
    do inhc = 1, Nnhc
      do iatom = 1, Natom
        !vrbath(1,iatom,inhc,imode) = vsigma*gasdev()
        !vrbath(2,iatom,inhc,imode) = vsigma*gasdev()
        !vrbath(3,iatom,inhc,imode) = vsigma*gasdev()
        do Idim_ = 1, Ndim
          vrbath(Idim_,iatom,inhc,imode) = vsigma*gasdev()
        end do

        rbath(:,iatom,inhc,imode) = 0.d0
      enddo
    enddo
  enddo

!YK how about one NH chain?
  select case(Ncent)
    case(1)
      do inhc = 1, nnhc
        vsigma = dsqrt(1.d0/beta/qmcent11(inhc))
        vbc11(inhc) = vsigma*gasdev()
        rbc11(inhc) = 0.d0
      enddo
    case(3)
      do inhc = 1, nnhc
        do iatom=1,natom
          vsigma = dsqrt(1.d0/beta/qmcent31(inhc))
          !vrbc31(1,iatom,inhc) = vsigma*gasdev()
          !vrbc31(2,iatom,inhc) = vsigma*gasdev()
          !vrbc31(3,iatom,inhc) = vsigma*gasdev()
          do Idim_ = 1, Ndim
            vrbc31(Idim_,iatom,inhc) = vsigma*gasdev()
          end do
          rbc31(:,iatom,inhc) = 0.d0
        enddo
      enddo
  end select

  return
end subroutine
