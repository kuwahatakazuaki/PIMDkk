Subroutine PI_NEW_MPI

  Use Parameters
  use Parameter_tk, only: laddress, addresstmp, address0
  Use MPI
  Implicit None
  Integer :: istep, iref

  Call Calc_Constant
  If(MyRank==0) Then
     Call Check_Inp
  EndIf
  Call Setup_MPI

! Kuwahata 2019/11/17 address2 is used for Scratch finle
!  laddress=len_trim(address)+1
!  addresstmp=trim(address0)

  address0=trim(address2)//'/'
  laddress=len_trim(address2)+1
  addresstmp=trim(trim(address2)//'/')
! End Kuwahata 2019/11/17

  IF(NForce>=5 .and. NForce<=10) Call Set_etc_MPI_tk

  select case(NForce)
    case(5)
      call Set_Gamess_MPI_tk
    case(6)
      call Set_Gaussian_MPI_tk  ! Set for chk, rwf, etc.
!    case(7)
!      call Set_Turbomole_MPI_tk
!    case(8)
!      call Set_VASP  ! We are making
!    case(9)
!      call Set_siesta
    case(11:12)
      print *, "Do I need something prepare??"
    case default
      stop 'ERROR!!! Wrong "NForce" option'
  end select


  Call Normal_Mode
  Call Init_Mass
  If(NRestart==1) Then
     If(MyRank==0) Then
        Call restart_read_old
     EndIf
     Call Broad3
  Else
!     if (NForce>=5) call print_ini_qm
     Call NM_Position
     If(MyRank==0) Then
        Call Init_Velocity
        Call Init_Bath
     EndIf
     Call Broad3
     Call Temp_ctr
     Call NM_Trans (0) ! x(i) = x(i) + sum_j tnm(i,j)*u(j)

     istepsv=0

     Call Force_New_MPI_tk
     Call Getfnm  ! fu(i) = fu(i) + sum_j fx(j)*tnm(j,i)

!    Call Vupdate 
!    If(NCent==1) Then
!       Call Nhc_Integrate_Cent 
!    EndIf
!    If(NCent==3) Then
!       Call Nhc_Integrate_Cent3
!    EndIf
     If(MyRank==0) Then
        Call Ham_Temp
        Call Print_Ham_tk(nrstep)
     EndIf
  EndIf

  Call Getforce_Ref
  if(myrank==0 .or. NForce<5) then

    Do istep=nrstep+1, nstep
      istepsv=istep
      If(NCent==1) Then
         Call Nhc_Integrate_Cent
      EndIf
      If(NCent==3) Then
         Call Nhc_Integrate_Cent3
      EndIf
      Call Vupdate

      Do iref=1, nref  ! nref = 5
         Call Nhc_Integrate
         Call Vupdate_Ref
         Call Uupdate
         Call Getforce_Ref
         Call Vupdate_Ref
         Call Nhc_Integrate
      Enddo

!     Call Getforce_Ref
      Call NM_Trans (0)     ! x(i) = x(i) + sum_j tnm(i,j)*u(j)
      Call Force_New_MPI_tk ! Obtaining fx
      if (umbrella_sampling .eqv. .True.) call calc_umbrella
      Call Getfnm           ! fu(i) = fu(i) + sum_j fx(j)*tnm(j,i)
      Call Vupdate
      If(NCent==1) Then
         Call Nhc_Integrate_Cent
      EndIf
      If(NCent==3) Then
         Call Nhc_Integrate_Cent3
      EndIf
      Call Ham_Temp
      If(MyRank==0) Then
         Call Print_Ham_tk(istep)
         Call Restart_Write(istep)
      EndIf
    Enddo
  else
    Do istep=nrstep+1, nstep
       istepsv=istep
       Call Force_New_MPI_tk
    Enddo
  endif

  If(NForce==5) Call Unset_Gamess_MPI_tk
  IF(NForce>=5 .and. NForce<=10) Call Unset_etc_MPI_tk
!  IF(NForce>=5 .and. NForce<=10) Call Close_files_MPI_tk

  if (MyRank == 0) then
    write(*,'(a)') &
      &'***********************************************************************&
      &************************************************************************'
  end if

  Return
  End Subroutine
