!subroutine normal_mode_matrix
!  use global_variable, only: Nbead, unitary, uni_inv, pi, &
!                             tnm, tnminv, myrank
!  implicit none
!  integer :: i, j
!  real(8) :: signf
!  real(8) :: sdbinv, pref
!  real(8) :: di, dj
!  real(8) :: temp
!  real(8) :: sqp, sqpinv
!
!  sdbinv = 1.0d0 / dsqrt( dble(Nbead) )
!  pref   = dsqrt(2.0d0 / dble(Nbead))
!
!  signf = -1.0d0
!  do i = 1, Nbead
!    unitary(i,1) = sdbinv
!    unitary(i,Nbead) = signf * sdbinv
!    signf = signf * (-1.0d0)
!  end do
!
!  do i = 1, (Nbead-2) / 2
!    di = dble(i)
!    do j = 1, Nbead
!      dj = dble(j)
!      temp = 2.0d0 * pi * di * dj / dble(Nbead)
!      unitary(j, 2*i  ) = pref * dcos(temp)
!      unitary(j, 2*i+1) = pref * dsin(temp)
!    end do
!  end do
!
!  do i = 1, Nbead
!    do j = 1, Nbead
!      uni_inv(j,i) = unitary(i,j)
!    end do
!  end do
!
!  sqp    = sqrt(dble(Nbead))
!  sqpinv = 1.0d0 / dsqrt(dble(Nbead))
!  tnm(:,:) = sqp * unitary(:,:)
!  tnminv(:,:) = sqpinv * uni_inv(:,:)
!!  unitary(:,:) = sqp    * unitary(:,:)
!!  uni_inv(:,:) - sqpinv * uni_inv(:,:)
!
!!if (myrank == 0) then
!!  do i = 1, Nbead
!!    do j = 1, Nbead
!!      print '(2I3,F10.5)', i,j, unitary(i,j)/pref
!!    end do
!!  end do
!!  print *, "pref is ", pref
!!end if
!
!end subroutine normal_mode_matrix



!  Double Precision     :: sqp, sqpinv, dnorm, dum
!
!!     /* parameters */
!  dp = dble(nbead)
!  sqp = dsqrt(dp)
!  sqpinv = 1.d0/dsqrt(dp)
!  dnorm = dsqrt(2.d0/dp)
!
!!     /*  making unitary matrix for diagnalizing the spring matrix *
!!      *  using analytical expressions                             */
!
!  dum = -1.d0
!  do imode = 1, nbead
!     u(imode,1) = sqpinv
!     u(imode,nbead) = dum*sqpinv
!     dum = dum*(-1.d0)
!  enddo
!
!  do imode = 1, (nbead-2)/2
!     di = dble(imode)
!     do jmode = 1, nbead
!        dj = dble(jmode)
!        u(jmode,2*imode)   = dnorm*dcos(2.d0*pi*di*dj/dp)
!        u(jmode,2*imode+1) = dnorm*dsin(2.d0*pi*di*dj/dp)
!     enddo
!  enddo
!
!  do imode = 1, nbead
!     do jmode = 1, nbead
!        uinv(jmode,imode) = u(imode,jmode)
!     enddo
!  enddo
!
!  do imode = 1, nbead
!     do jmode = 1, nbead
!        tnm(imode,jmode) = sqp*u(imode,jmode)
!        tnminv(imode,jmode) = sqpinv*uinv(imode,jmode)
!     enddo
!  enddo

