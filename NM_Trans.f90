Subroutine NM_Trans(NTmp)

  Use Parameters
  Implicit None
  Integer       :: NTmp, jmode
  integer :: i, j


  if (NTmp ==  0) then
!     /*  from normal mode variables to real variables  *
!      *  x(i) = x(i) + sum_j tnm(i,j)*u(j)             */
!     /*  initialize array  */
     do imode = 1, nbead
        do iatom = 1, natom
           x(iatom,imode) = 0.d0
           y(iatom,imode) = 0.d0
           z(iatom,imode) = 0.d0
        enddo
     enddo

     do iatom = 1, natom
        do imode = 1, nbead
           do jmode = 1, nbead
              x(iatom,imode) = x(iatom,imode) + tnm(imode,jmode)*ux(iatom,jmode)
              y(iatom,imode) = y(iatom,imode) + tnm(imode,jmode)*uy(iatom,jmode)
              z(iatom,imode) = z(iatom,imode) + tnm(imode,jmode)*uz(iatom,jmode)
           enddo
        enddo
     enddo

  elseif (NTmp == 1) then  ! We don't use this
!     /*  from real variables to normal mode variables  *
!      *  u(i) = u(i) + sum_j tnminv(i,j)*x(j)          */
!     /*  initialize array  */
     do imode = 1, nbead
        do iatom = 1, natom
           ux(iatom,imode) = 0.d0
           uy(iatom,imode) = 0.d0
           uz(iatom,imode) = 0.d0
        enddo
     enddo

     do iatom = 1, natom
        do imode = 1, nbead
           do jmode = 1, nbead
              ux(iatom,imode) = ux(iatom,imode) + tnminv(imode,jmode)*x(iatom,jmode)
              uy(iatom,imode) = uy(iatom,imode) + tnminv(imode,jmode)*y(iatom,jmode)
              uz(iatom,imode) = uz(iatom,imode) + tnminv(imode,jmode)*z(iatom,jmode)
           enddo
        enddo
     enddo
  endif

  Return
End Subroutine
