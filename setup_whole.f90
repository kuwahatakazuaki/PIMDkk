subroutine setup_whole
  use global_variable
  use mpi
  implicit none
  integer :: i, j, icolor, inhc
  real(8) :: rndnumber
  real(8) :: omega
!  real(8) :: factor_color

  call random_generator(1,rndnumber)

  beta = 1.0d0 / (KtoAU * temperature)
  dt_ref   = dt / dble(Nref)
!  dp       = dble(Nbead)
  omega_p2 = dble(Nbead) / (beta*beta)
  omega_system = 2.0d0 * pi / freq1 / facttime
  omega2   = omega_system * omega_system

  gnkt  = 3.0d0 * dble(Natom) / beta
  gkt   = 1.0d0 / beta


  if      ( Ncent == 1 ) then

    if      ( Ncolor == 1) then

      qmcent11(1) = 3.0d0 * dble(Natom) / beta / omega2
      qmcent11(2:Nnhc) = 1.0d0 / beta / omega2

    else if ( Ncolor >= 2 ) then

      do icolor = 1, Ncolor
        omega        = (omega_system * dble(Ncolor) / dble(Ncolor+1-icolor) )
        omega2       = omega * omega
        qmcent1(1     ,icolor) = 3.0d0 * dble(Natom) / (beta * omega2)
        qmcent1(2:Nnhc,icolor) = 1.0d0 / (beta * omega2)
        qmcent1(1     ,Ncolor) = 3.0d0 * dble(Natom) / (beta * omega_p2)
        qmcent1(2:Nnhc,Ncolor) = 1.0d0 / (beta * omega_p2)
      end do

    end if

  else if (Ncent == 3) then

    if      ( Ncolor == 1) then


!      qmcent31(1) = 3.0d0 * dble("Nbead"?) / beta / omega2 Kuwahata 2020/02/05
      qmcent31(1) = 3.0d0 * dble(Natom) / beta / omega2
      qmcent31(2:Nnhc) =  1.0d0 / beta / omega2

    else if ( Ncolor >= 2 ) then

      do icolor = 1, Ncolor
        omega        = (omega_system * dble(Ncolor) / dble(Ncolor+1-icolor))
        omega2       = omega * omega
        qmcent3(1     ,icolor) = 3.0d0 * dble(Natom) / (beta * omega2)
        qmcent3(2:Nnhc,icolor) = 1.0d0 / (beta * omega2)
      end do
      qmcent3(1     ,Ncolor) = 3.0d0 * dble(Natom) / (beta * omega_p2)
      qmcent3(2:Nnhc,Ncolor) = 1.0d0 / (beta * omega_p2)

    end if
  end if

! qmass is the mass of thermostat. we don't use qmass(1) ??
  qmass(1) = 0.0d0
  if (Simulation == 3) then  ! CMD or RPMD
    qmass(2:Nbead) = 1.0d0 / beta / omega_p2 * gamma * gamma
  else
    qmass(2:Nbead) = 1.0d0 / beta / omega_p2
  end if

  if      (Nys == 1) then
    ysweight(1) = 1.d0
  else if (Nys == 3) then
    ysweight(1) = 1.d0/(2.d0 - 2.d0**(1.d0/3.d0))
    ysweight(2) = 1.d0 - 2.d0*ysweight(1)
    ysweight(3) = ysweight(1)
  else if (Nys == 5) then
    ysweight(1) = 1.d0/(4.d0 - 4.d0**(1.d0/3.d0))
    ysweight(2) = ysweight(1)
    ysweight(3) = 1.d0 - 4.d0*ysweight(1)
    ysweight(4) = ysweight(1)
    ysweight(5) = ysweight(1)
  end if


end subroutine setup_whole

!
!!YK Print Some Constants Used in Calculation
!if (MyRank == 0) then
!    write(*,*) 
!!    write(*,*) 'MyRank = ', MyRank
!    write(*,*) ' +++++   Calculate Some Constants   +++++ '
!    Write(*,'(a20,D15.9)') ' Random Number    = ',rndnumber  ! Now we don't use rndnumber
!    write(*,'(a20,D15.9)') ' dp               = ',dp
!    write(*,'(a20,D15.9)') ' dt               = ',dt
!    write(*,'(a20,D15.9)') ' dt_ref           = ',dt_ref
!    write(*,'(a20,D15.9)') ' omega_p          = ',omega_p
!    write(*,'(a20,D15.9)') ' omega_p2         = ',omega_p2
!    write(*,'(a20,D15.9)') ' omega_system     = ',omega_system
!    write(*,'(a20,D15.9)') ' omega2           = ',omega2
!    write(*,'(a20,D15.9)') ' gnkt             = ',gnkt
!    write(*,'(a20,D15.9)') ' gkt              = ',gkt
!    write(*,'(a20,D15.9)') ' ysweight1        = ',ysweight(1)
!    write(*,'(a20,D15.9)') ' ysweight2        = ',ysweight(2)
!    write(*,'(a20,D15.9)') ' ysweight3        = ',ysweight(3)
!    write(*,'(a20,D15.9)') ' ysweight4        = ',ysweight(4)
!    write(*,'(a20,D15.9)') ' ysweight5        = ',ysweight(5)
!    write(*,*) ' ++++++++++++++++++++++++++++++++++++++++ '
!    write(*,*)
!end if
!
!Return
!End Subroutine
