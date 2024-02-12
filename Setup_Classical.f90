Subroutine Setup_Classical
  use Parameters
  implicit none
  double precision :: rndnumber
  integer :: inhc, icolor, imode, iatom

!YK Initiate Random Number Generator
  Call RandomG(0,rndnumber)

!     /*   for multiple time step   */
!
!   We must be careful in RPMD!
!
  dt_ref = dt
  ista = 1
  iend = 1
!
!     /*   bath parameters for path integral MD   */
!
!YK different qmass required for centroid MD
   If(NCent==1) Then

      !If(NColor==1) Then
         qmcent11(1) = 3.d0*dble(natom)/beta/omega2
         do inhc=2,nnhc
            qmcent11(inhc) = 1.d0/beta/omega2
         enddo
      !EndIf

      !If(NColor>=2) Then
      !   do icolor=1,NColor
      !      factor_color=dble(NColor-icolor+1)/dble(ncolor)
      !      omega2=(omega_system*factor_color)**2
      !      qmcent1(1,icolor) = 3.d0*dble(natom)/beta/omega2
      !      do inhc=2,nnhc
      !         qmcent1(inhc,icolor) = 1.d0/beta/omega2
      !      enddo
      !   enddo
      !   qmcent1(1,ncolor) = 3.d0*dble(natom)/beta/omega_p2
      !   do inhc=2,nnhc
      !      qmcent1(inhc,ncolor) = 1.d0/beta/omega_p2
      !   enddo
      !EndIf

   EndIf

   If(NCent==3) Then

      !If(NColor==1) Then
         qmcent31(1) = 3.d0*dble(natom)/beta/omega2
         do inhc=2,nnhc
            qmcent31(inhc) = 1.d0/beta/omega_p2
         enddo
      !EndIf

      !If(NColor>=2) Then
      !   do icolor=1,NColor
      !      factor_color=dble(NColor-icolor+1)/dble(ncolor)
      !      omega2=(omega_system*factor_color)**2
      !      qmcent3(1,icolor) = 3.d0*dble(natom)/beta/omega2
      !      do inhc=2,nnhc
      !         qmcent3(inhc,icolor) = 1.d0/beta/omega2
      !      enddo
      !   enddo
      !   qmcent3(1,ncolor) = 3.d0*dble(natom)/beta/omega_p2
      !   do inhc=2,nnhc
      !      qmcent3(inhc,ncolor) = 1.d0/beta/omega_p2
      !   enddo
      !EndIf

   EndIf


!     /*   set parameters for higher order decomposition   *
!      *   of propagator (Yoshida and Suzuki parameters)   */

  if(nys==1) then
     ysweight(1) = 1.d0
  else if (nys==3) then
     ysweight(1) = 1.d0/(2.d0 - 2.d0**(1.d0/3.d0))
     ysweight(2) = 1.d0 - 2.d0*ysweight(1)
     ysweight(3) = ysweight(1)
  else if (nys==5) then
     ysweight(1) = 1.d0/(4.d0 - 4.d0**(1.d0/3.d0))
     ysweight(2) = ysweight(1)
     ysweight(3) = 1.d0 - 4.d0*ysweight(1)
     ysweight(4) = ysweight(1)
     ysweight(5) = ysweight(1)
  end if

  do iatom = 1, natom
     fictmass(iatom,1) = physmass(iatom)
  enddo

  Return
End Subroutine
