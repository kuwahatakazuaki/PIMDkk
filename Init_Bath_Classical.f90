Subroutine Init_Bath_Classical
  Use Parameters
  Implicit None
  Double Precision :: gasd
  integer :: inhc, imode, iatom, icolor

!     /*  single thermostat attached to centroids  */
!YK how about one NH chain?
  !If(NCent==1) Then
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
  !EndIf ! NCent==1
  !If(NCent==3) Then
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
  !EndIf ! NCent==3
  end select

  Return
End Subroutine
