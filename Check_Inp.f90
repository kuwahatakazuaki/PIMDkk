Subroutine Check_Inp

  Use Parameters
  Use Parameter_tk
  Implicit None

!YK added some new parameters for CMD
  Write (*,'(a)')' +++++ Input Check +++++   '
  Write (*,9999) ' +++++ Simulation Typ      ', Simulation
  Write (*,9999) ' +++++ Number of Atoms     ', natom
  Write (*,9999) ' +++++ Number of Beads     ', nbead
  Write (*,9999) ' +++++ Number of Steps     ', nstep
  Write (*,9998) ' +++++ Given Temperature   ', Temperature
  Write (*,9998) ' +++++ Given time step     ', dt
  Write (*,9998) ' +++++ Beta                ', beta
  Write (*,9998) ' +++++ Adiabaticity param. ', gamma
  Write (*,9999) ' +++++ Type of Ensemble    ', nensemble
  Write (*,9999) ' +++++ Method of QM calc.  ', theory
  Write (*,9999) ' +++++ QM/MM calculation?  ', nqmmm
  Write (*,9999) ' +++++ Additional Basis?   ', ngengau
  Write (*,9999) ' +++++ Frequency of vel sv ', nsavevel
  Write (*,9999) ' +++++ Frequency of out sv ', nsavelog
  Write (*,9999) ' +++++ Frequency of chk sv ', nsavechk
  Write (*,9999) ' +++++ Method of Centr NHC ', ncent
  Write (*,9999) ' +++++ Color of Centr NHC  ', ncolor
  Write (*,9999) ' +++++ Length of Centr NHC ', nnhc
  Write (*,9999) ' +++++ Flag for Restart    ', nrestart
  Write (*,9999) ' +++++ Flag for Force Calc ', nforce
  Write (*,9999) ' +++++ Seed for Random No.1', iseed1
  Write (*,9999) ' +++++ Seed for Random No.2', iseed2
  Write (*,9999) ' +++++ Seed for Random No.3', iseed3
  Write (*,9999) ' +++++ Seed for Random No.4', iseed4
!  Write (*,9997) ' +++++ Address of Scratch  ', address
!  Write (*,9997) ' +++++ Address of Scratch2 ', address2
! Kuwahata 2019/11/16
  Write (*,9997) ' +++++ Address of Result   ', trim(address)
  Write (*,9997) ' +++++ Address of Scratch  ', trim(address2)
! End Kuwahata 2019/11/16
  Write (*,9999) ' +++++ proc for extra out  ', kproc
  Write (*,*)
  Write (*,*) ' +++++ Given Atomic Label, Mass, and Coords +++++'
  Do iatom=1, natom
     Write(*,9996) alabel(iatom),PhysMass(iatom),Ux(iatom,1), Uy(iatom,1), Uz(iatom,1)
  Enddo

  Return
 9999 Format(A26,I15)
 9998 Format(A26,F15.9)
 9997 Format(A26,8X,A)
! 9997 Format(A26,A80)
 9996 Format(A3,1X,4F15.9)
End Subroutine
