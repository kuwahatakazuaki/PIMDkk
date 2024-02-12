  Subroutine RPMD_NEW

    Use Parameters
    Implicit None

    Integer                            :: istep, iref 

!
    Call Calc_Constant
    Call Check_Inp
    Call Setup 
    Call Normal_Mode 
    Call Init_Mass 
    If(NRestart==1) Then
       Call restart_read
    Else
       Call NM_Position 
       Call Init_Velocity
       If(NEnsemble==1) Then
          Call Init_Bath
       EndIf
       Call Temp_ctr 
!    Call Getforce_Ref
       Call NM_Trans (0)
       istepsv=0
       Call Force_New
       Call Getfnm 
!    Call Vupdate 

!    If(Nensemble==1) Then
!       Call Nhc_Integrate_Cent 
!    EndIf
       Call Ham_Temp 
       Call Print_Ham(0,nrstep)
    EndIf
!    Call Print_Out
!
!    Path Integral Loop
!!!HERE!!!
    Call Getforce_Ref
!!!HERE!!!
    Do istep=nrstep+1, nstep 
       istepsv=istep
       If(Nensemble==1) Then
          Call Nhc_Integrate_Cent 
       EndIf
       Call Vupdate 
!
!       Nref Loop 
       Do iref=1, nref
          If(NEnsemble==1) Then
             Call Nhc_Integrate
          EndIf
          Call Vupdate_Ref 
          Call Uupdate 
          Call Getforce_Ref
          Call Vupdate_Ref 
          If(NEnsemble==1) Then
             Call Nhc_Integrate
          EndIf
       Enddo
       Call Getforce_Ref
!       Nref Loop End
!
       Call NM_Trans (0)
       Call Force_NEW
       Call Getfnm 
       Call Vupdate 
       If(Nensemble==1) Then
          Call Nhc_Integrate_Cent 
       EndIf
       Call Ham_Temp 
       Call Print_Ham(1,istep)
       Call Restart_Write(istep)
!      Call Analyze(istep)

    Enddo

    Return
  End Subroutine 
