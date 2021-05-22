subroutine check_input
  use global_variable
  implicit none
  integer :: i

!  write(*,'(a)')  ' +++++ Input Check +++++   '
  write(*,'(a)')  ' +++++ Input Check +++++++++++++++++++++++++++'
  write(*,10000)  ' +++++ Process ID          ', getpid()
  write(*,10000)  ' +++++ Simulation type     ', simulation
  write(*,10000)  ' +++++ Number of Atoms     ', Natom
  write(*,10000)  ' +++++ Number of Beads     ', Nbead
  write(*,10000)  ' +++++ Number of Steps     ', Nstep
  write(*,10001)  ' +++++ Given Temperature   ', temperature
  write(*,10001)  ' +++++ Given time step     ', dt
  write(*,10000)  ' +++++ Type of Ensemble    ', Nensemble
  write(*,10000)  ' +++++ Method of QM calc.  ', Ntheory
  write(*,10000)  ' +++++ Method of Center NHC', Ncent
  write(*,10000)  ' +++++ Color of Center NHC ', Ncolor
  write(*,10000)  ' +++++ Length of Center NHC', Nnhc
  write(*,10002)  ' +++++ Flag for Restart    ', Lrestart
  write(*,10000)  ' +++++ Flag for Force Calc.', Nforce
  write(*,10002)  ' +++++ QM/MM calculation?  ', Lqmmm
  write(*,10002)  ' +++++ Additional Basis?   ', Lgengau
  write(*,10000)  ' +++++ Seed for Random No.1', iseed(1)
  write(*,10000)  ' +++++ Seed for Random No.2', iseed(2)
  write(*,10000)  ' +++++ Seed for Random No.3', iseed(3)
  write(*,10000)  ' +++++ Seed for Random No.4', iseed(4)
  write(*,10003)  ' +++++ Address of Result   ', path_result
  write(*,10003)  ' +++++ Address of Scratch  ', path_scr
  write(*,'()')
  write(*,'(a)')  ' +++++   Calculate Some Constants   ++++++++++ '
!  write(*,10001)  ' +++++ Random Number  ', rndnumber
!  write(*,10001)  ' +++++ dp            ', dp
  write(*,10001)  ' +++++ dt            ', dt
  write(*,10001)  ' +++++ dt_ref        ', dt_ref
!  write(*,10001)  ' +++++ omega_p       ', omega_p
  write(*,10005)  ' +++++ omega_system  ', omega_system
  write(*,10005)  ' +++++ omega2        ', omega2
  write(*,10005)  ' +++++ omega_p2      ', omega_p2
  write(*,10001)  ' +++++ gnkt          ', gnkt
  write(*,10001)  ' +++++ gkt           ', gkt
  write(*,10001)  ' +++++ ysweight1     ', ysweight(1)
  write(*,10001)  ' +++++ ysweight2     ', ysweight(2)
  write(*,10001)  ' +++++ ysweight3     ', ysweight(3)
  write(*,10001)  ' +++++ ysweight4     ', ysweight(4)
  write(*,10001)  ' +++++ ysweight5     ', ysweight(5)
  write(*,'()')
  write(*,'(a)')  ' +++++ Given Atomic Label, Mass, and Coords +++'
  do i = 1, Natom
    write(*,10004) alabel(i), physmass(i)/factmass, u(:,i,1) * AUtoAng
  end do
  write(*,'()')
  write(*,'(a)')  ' +++++ End Input Check ++++++++++++++++++++++++'
  write(*,'()')

return
10000 format(a,I7)
10002 format(a,L7)
10001 format(a,F11.5)
10005 format(a,G15.5)
10003 format(a,5x,a)
10004 format(8x,a,2x,4F11.6)

end subroutine check_input

!  Write (*,'(a)')' +++++ Input Check +++++   '
!  Write (*,9999) ' +++++ Flag for Force Calc ', nforce
!!  Write (*,9997) ' +++++ Address of Scratch  ', address
!!  Write (*,9997) ' +++++ Address of Scratch2 ', address2
!! Kuwahata 2019/11/16
!  Write (*,9997) ' +++++ Address of Result   ', trim(address)
!  Write (*,9997) ' +++++ Address of Scratch  ', trim(address2)
!! End Kuwahata 2019/11/16
!  Write (*,9999) ' +++++ proc for extra out  ', kproc
!  Write (*,*)
!  Write (*,*) ' +++++ Given Atomic Label, Mass, and Coords +++++'
!  Do iatom=1, natom
!     Write(*,9996) alabel(iatom),PhysMass(iatom),Ux(iatom,1), Uy(iatom,1), Uz(iatom,1)
!  Enddo
!  Write (*,9998) ' +++++ Beta                ', beta
!  Write (*,9998) ' +++++ Adiabaticity param. ', gamma
!  Write (*,9999) ' +++++ Frequency of vel sv ', nsavevel
!  Write (*,9999) ' +++++ Frequency of out sv ', nsavelog
!  Write (*,9999) ' +++++ Frequency of chk sv ', nsavechk
!
!  Return
! 9999 Format(A26,I15)
! 9998 Format(A26,F15.9)
! 9997 Format(A26,8X,A)
!! 9997 Format(A26,A80)
! 9996 Format(A3,1X,4F15.9)
!End Subroutine
