Subroutine CMD_NEW_MPI

  Use Parameters
  Use MPI
  Implicit None
  Integer  :: istep, iref

  Call Calc_Constant
  If(MyRank==0) Then
    Call Check_Inp
  EndIf
  Call Setup_MPI
  Call Normal_Mode
  Call Init_Mass
  If(NRestart==1) Then
     If(MyRank==0) Then
        Call Restart_Read
     EndIf
     Call Broad3
  Else
     Call NM_Position
     Call Init_Velocity
     Call Init_Bath
     Call Broad3
     Call Temp_ctr
!   Call Getforce_Ref
     Call NM_Trans (0)
     istepsv=0
     Call Force_New_MPI_tk
     Call Getfnm
!   Call Vupdate

!    If(Nensemble==1) Then
!       Call Nhc_Integrate_Cent 
!    EndIf
     Call Ham_Temp
     If(MyRank==0) Then
        Call Print_Ham(0,nrstep)
     EndIf
  EndIf

  Call Getforce_Ref
  Do istep=nrstep+1, nstep 
     istepsv=istep
     If(Nensemble==1) Then
        Call Nhc_Integrate_Cent 
     EndIf
     Call Vupdate 

     Do iref=1, nref
        Call Nhc_Integrate
        Call Vupdate_Ref 
        Call Uupdate 
        Call Getforce_Ref
        Call Vupdate_Ref 
        Call Nhc_Integrate

     Enddo

     Call Getforce_Ref
     Call NM_Trans (0)
     Call Force_NEW_MPI_tk
     Call Getfnm 
     Call Vupdate 
     If(Nensemble==1) Then
        Call Nhc_Integrate_Cent 
     EndIf
     Call Ham_Temp 
     If(MyRank==0) Then
        Call Print_Ham(1,istep)
        Call Restart_Write(istep)
     EndIf

  Enddo

Return
End Subroutine
