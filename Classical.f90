Subroutine Classical

  Use Parameters
! koba add (2017.11.16) for address0
  use Parameter_tk, only: laddress, addresstmp, address0
!  Use Parameter_tk

  Implicit None
  Integer        :: istep, iref

! Kuwahata 2020/01/11
!  address0=trim(address)//'/'
!  laddress=len_trim(address)+1
!  addresstmp=trim(address0)
  address0=trim(address2)//'/'
  laddress=len_trim(address2)+1
  addresstmp=trim(address0)
! End Kuwahata 2020/01/11

  Call Calc_Constant
  Call Check_Inp
  Call Setup_Classical
!tkawatsu
!    IF(NForce>=5) Call Set_etc_MPI_tk
! koba imported from Set_etc_MPI_tk.F90 (2017.11.16)
!    bohr_inv = 1.0d0/bohr
    dp_inv = 1.0d0/dp

! Kuwahata 2019/08/01 We no longer use open file.
!       The files are open in 'Force_Gaussian.f90'
!    IF(NForce>=5) Call Open_files_tk
! End Kuwahata

!    If(NForce==5) Call Set_Gamess_MPI_tk
!    If(NForce==6) Call Set_Gaussian_MPI_tk
    !If(NForce==7) Call Set_Turbomole_tk
!
  If(NRestart==1) Then
     Call restart_read_classical
  Else
! Kuwahata 2020/01/12
!     call print_ham_classical_ini
     call print_ini_cl
! End Kuwahata 2020/01/12
     Call Init_Velocity
     Call Init_Bath_Classical
     Call Temp_ctr
!   Call Getforce_Ref
     istepsv=0
     Call Force_Classical
!    Call Vupdate 
!    If(NCent==1) Then
!       Call Nhc_Integrate_Cent 
!    EndIf
!    If(NCent==3) Then
!       Call Nhc_Integrate_Cent3
!    EndIf
     Call Ham_Temp_Classical
     Call Print_Ham_Classical(nrstep)
!       Call Print_Ham_Classical(0,nrstep)
  EndIf

!    Path Integral Loop
  Do istep=nrstep+1, nstep
     istepsv=istep
     If(NCent==1) Then
        Call Nhc_Integrate_Cent
     EndIf
     If(NCent==3) Then
        Call Nhc_Integrate_Cent3
     EndIf
     Call Vupdate
     Call Uupdate

!    Call Getforce_Ref
     Call Force_Classical
     Call Vupdate
     If(NCent==1) Then
        Call Nhc_Integrate_Cent
     EndIf
     If(NCent==3) Then
        Call Nhc_Integrate_Cent3
     EndIf
     Call Ham_Temp_Classical
     Call Print_Ham_Classical(istep)
     Call Restart_Write_Classical(istep)
  EndDo

!tkawatsu
!    If(NForce==5) Call Unset_Gamess_MPI_tk
!    IF(NForce>=5) Call Unset_etc_MPI_tk
! Kuwahata 20190801 We no longer use open file
!    IF(NForce>=5) Call Close_files_tk
! kuwahata


Return
End Subroutine
