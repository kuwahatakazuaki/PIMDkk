  Subroutine Init_Bath_Classical

    Use Parameters
    Implicit None
    Double Precision :: gasd

!     /*  single thermostat attached to centroids  */
!YK how about one NH chain?
      If(NCent==1) Then

         If(NColor==1) Then

            do inhc = 1, nnhc
               vsigma = dsqrt(1.d0/beta/qmcent11(inhc))
               call gasdev(gasd)
               vbc11(inhc) = vsigma*gasd
               rbc11(inhc) = 0.d0
            enddo
            close(55)
         Else ! Ncolor==1

            do icolor=1,ncolor
               do inhc = 1, nnhc
                  vsigma = dsqrt(1.d0/beta/qmcent1(inhc,icolor))
                  call gasdev(gasd)
                  vbc1(inhc,icolor) = vsigma*gasd
                  rbc1(inhc,icolor) = 0.d0
               enddo
            enddo

         EndIf ! NColor==1

      EndIf ! NCent==1

      If(NCent==3) Then

         If(NColor==1) Then

            do inhc = 1, nnhc
               do iatom=1,natom
                  vsigma = dsqrt(1.d0/beta/qmcent31(inhc))
                  call gasdev(gasd)
                  vxbc31(iatom,inhc) = vsigma*gasd
                  call gasdev(gasd)
                  vybc31(iatom,inhc) = vsigma*gasd
                  call gasdev(gasd)
                  vzbc31(iatom,inhc) = vsigma*gasd
                  xbc31(iatom,inhc) = 0.d0
                  ybc31(iatom,inhc) = 0.d0
                  zbc31(iatom,inhc) = 0.d0
               enddo
            enddo

         Else ! NColor==1

            do icolor=1,ncolor
               do inhc = 1, nnhc
                  do iatom=1,natom
                     vsigma = dsqrt(1.d0/beta/qmcent3(inhc,icolor))
                     call gasdev(gasd)
                     vxbc3(iatom,inhc,icolor) = vsigma*gasd
                     call gasdev(gasd)
                     vybc3(iatom,inhc,icolor) = vsigma*gasd
                     call gasdev(gasd)
                     vzbc3(iatom,inhc,icolor) = vsigma*gasd
                     xbc3(iatom,inhc,icolor) = 0.d0
                     ybc3(iatom,inhc,icolor) = 0.d0
                     zbc3(iatom,inhc,icolor) = 0.d0
                  enddo
               enddo
            enddo

         EndIf ! NColor==1

      EndIf ! NCent==3

      Return
  End Subroutine
