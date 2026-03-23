subroutine Init_Bath
  use Parameters
  use utility, only: fill_gaussian_random
  Implicit None
  real(8) :: vsigma
  integer :: inhc, imode

!     /*  massive Nose-Hoover chain  */
!YK Remove the initial velocity for the centroid 
  rbath(:,:,:,:) = 0.d0
  vrbath(:,:,:,1) = 0.d0

  do imode = 2, Nbead
    vsigma = dsqrt(1.d0/beta/qmass(imode))
    call fill_gaussian_random(vrbath(:,:,:,imode))
    vrbath(:,:,:,imode) = vsigma * vrbath(:,:,:,imode)
  enddo

!YK how about one NH chain?
  select case(Ncent)
    case(1)
      call fill_gaussian_random(vbc11)
      vbc11(:) = dsqrt(1.d0/beta/qmcent11(:)) * vbc11(:)
      rbc11(:) = 0.d0
    case(3)
      do inhc = 1, nnhc
        vsigma = dsqrt(1.d0/beta/qmcent31(inhc))
        call fill_gaussian_random(vrbc31(:,:,inhc))
        vrbc31(:,:,inhc) = vsigma * vrbc31(:,:,inhc)
      enddo
      rbc31(:,:,:) = 0.d0
  end select

  return
end subroutine
