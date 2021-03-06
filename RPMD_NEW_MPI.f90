  Subroutine RPMD_NEW_MPI

    Use Parameters
    Use MPI
    Implicit None

    Integer                            :: istep, iref 

!
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
!    Call Getforce_Ref
       Call NM_Trans (0)
       istepsv=0
       Call Force_New_MPI
       Call Getfnm 
!    Call Vupdate 

!    If(Nensemble==1) Then
!       Call Nhc_Integrate_Cent 
!    EndIf
       Call Ham_Temp 
       If(MyRank==0) Then
          Call Print_Ham(0,nrstep)
       EndIf
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
          If(Nensemble==1) Then
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
       Call NM_Trans (0)
       Call Force_NEW_MPI
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
