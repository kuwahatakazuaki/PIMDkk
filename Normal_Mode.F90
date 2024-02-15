Subroutine Normal_Mode
  Use Parameters
  use utility, only: program_abort
  Implicit None

  Double Precision     :: sqp, sqpinv, dnorm, dum
  Double Precision     :: di, dj
  integer :: i, j, jmode, imode

!     /* parameters */
  dp = dble(nbead)
  sqp = dsqrt(dp)
  sqpinv = 1.d0/dsqrt(dp)
  dnorm = dsqrt(2.d0/dp)

!     /*  making unitary matrix for diagnalizing the spring matrix *
!      *  using analytical expressions                             */

  dum = -1.d0
  do imode = 1, nbead
     u(imode,1) = sqpinv
     u(imode,nbead) = dum*sqpinv
     dum = dum*(-1.d0)
  enddo

  do imode = 1, (nbead-2)/2
     di = dble(imode)
     do jmode = 1, nbead
        dj = dble(jmode)
        u(jmode,2*imode) = dnorm*dcos(2.d0*pi*di*dj/dp)
        u(jmode,2*imode+1) = dnorm*dsin(2.d0*pi*di*dj/dp)
     enddo
  enddo

  do imode = 1, nbead
     do jmode = 1, nbead
        uinv(jmode,imode) = u(imode,jmode)
     enddo
  enddo

  do imode = 1, nbead
     do jmode = 1, nbead
        tnm(imode,jmode) = sqp*u(imode,jmode)
        tnminv(imode,jmode) = sqpinv*uinv(imode,jmode)
     enddo
  enddo

Return
End Subroutine
