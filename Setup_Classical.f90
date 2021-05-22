  Subroutine Setup_Classical

    Use Parameters
    Implicit None

    double precision :: rndnumber

!      call system('mkdir '//trim(address))
!      call system('cp *.* '//trim(address))
!
!YK Initiate Random Number Generator
    Call RandomG(0,rndnumber)

!     /*   for multiple time step   */
!
!   We must be careful in RPMD!
!
    dt_ref = dt
!
!     /*   parameters for path integral simulation   */
!
     dp = 1.0D+00
     omega_p = dsqrt(dp)/beta
     omega_p2 = omega_p*omega_p
     omega2 = omega_system*omega_system/dp
     omega2 = omega_system*omega_system
!
!     /*   bath parameters for path integral MD   */
!
!YK different qmass required for centroid MD
     If(NCent==1) Then
 
        If(NColor==1) Then
           qmcent11(1) = 3.d0*dble(natom)/beta/omega2
           do inhc=2,nnhc
              qmcent11(inhc) = 1.d0/beta/omega2
           enddo
        EndIf
  
        If(NColor>=2) Then
           do icolor=1,NColor
              factor_color=dble(NColor-icolor+1)/dble(ncolor)
              omega2=(omega_system*factor_color)**2
              qmcent1(1,icolor) = 3.d0*dble(natom)/beta/omega2
              do inhc=2,nnhc
                 qmcent1(inhc,icolor) = 1.d0/beta/omega2
              enddo
           enddo
           qmcent1(1,ncolor) = 3.d0*dble(natom)/beta/omega_p2
           do inhc=2,nnhc
              qmcent1(inhc,ncolor) = 1.d0/beta/omega_p2
           enddo
        EndIf
 
     EndIf
 
     If(NCent==3) Then
 
        If(NColor==1) Then
           qmcent31(1) = 3.d0*dble(natom)/beta/omega2
           do inhc=2,nnhc
              qmcent31(inhc) = 1.d0/beta/omega_p2
           enddo
        EndIf
 
        If(NColor>=2) Then
           do icolor=1,NColor
              factor_color=dble(NColor-icolor+1)/dble(ncolor)
              omega2=(omega_system*factor_color)**2
              qmcent3(1,icolor) = 3.d0*dble(natom)/beta/omega2
              do inhc=2,nnhc
                 qmcent3(inhc,icolor) = 1.d0/beta/omega2
              enddo
           enddo
           qmcent3(1,ncolor) = 3.d0*dble(natom)/beta/omega_p2
           do inhc=2,nnhc
              qmcent3(inhc,ncolor) = 1.d0/beta/omega_p2
           enddo
        EndIf
 
     EndIf
!YK
!     /*   For centroid MD, qmass should be scaled by gamma2, since   *
!      *   natural frequencies of modes are omega_p**2/gamma**2       */
!
    gnkt = 3.d0*dble(natom)/beta
    gkt = 1.d0/beta
!
!     /*   set parameters for higher order decomposition   *
!      *   of propagator (Yoshida and Suzuki parameters)   */
!
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
!YK Print Some Constants Used in Calculation
    write(*,*)
    write(*,*) ' +++++   Calculate Some Constants   +++++ '
    Write(*,'(a20,D15.9)') ' Random Number    = ',rndnumber
    write(*,'(a20,D15.9)') ' dt               = ',dt
    write(*,'(a20,D15.9)') ' dt_ref           = ',dt_ref
    write(*,'(a20,D15.9)') ' omega_p          = ',omega_p
    write(*,'(a20,D15.9)') ' omega_p2         = ',omega_p2
    write(*,'(a20,D15.9)') ' omega_system     = ',omega_system
    write(*,'(a20,D15.9)') ' omega2           = ',omega2
    write(*,'(a20,D15.9)') ' gnkt             = ',gnkt
    write(*,'(a20,D15.9)') ' gkt              = ',gkt
    write(*,'(a20,D15.9)') ' ysweight1        = ',ysweight(1)
    write(*,'(a20,D15.9)') ' ysweight2        = ',ysweight(2)
    write(*,'(a20,D15.9)') ' ysweight3        = ',ysweight(3)
    write(*,'(a20,D15.9)') ' ysweight4        = ',ysweight(4)
    write(*,'(a20,D15.9)') ' ysweight5        = ',ysweight(5)
    write(*,*) ' ++++++++++++++++++++++++++++++++++++++++ '
    write(*,*)
 
    Return
  End Subroutine
