Subroutine Init_Bath_Classical
  use Parameters
  use utility, only: fill_gaussian_random
  Implicit None
  real(8) :: vsigma
  integer :: inhc

!     /*  single thermostat attached to centroids  */
!YK how about one NH chain?
  select case(Ncent)
    case(0)
      continue
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

  Return
End Subroutine
