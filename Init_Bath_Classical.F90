Subroutine Init_Bath_Classical
  Use Parameters
  Implicit None
  Double Precision :: gasd, vsigma
  integer :: inhc, imode, iatom, icolor

!     /*  single thermostat attached to centroids  */
!YK how about one NH chain?
  select case(Ncent)
    case(0)
      continue
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
