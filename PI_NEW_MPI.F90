subroutine PI_NEW_MPI
  use Parameters
  use utility, only: program_abort
  implicit none
  integer :: iref, Uout

  call Setup_time_mass
  call set_pallarel
  call Normal_Mode
  call Init_Mass

  !call set_each_exe

  select case(Iforce)
    case(1)
      call Set_mopac
    case(6)
      call Set_Gaussian_MPI_tk  ! Set for chk, rwf, etc.
    case(8)
      call Set_VASP
    case(9)
      call Set_siesta
    !case(10)
    !  call Set_Molcas_MPI
    case(21)
      call set_nnp_araidai
    case(22)
      call set_nnp_matlantis
  end select

  if ( Lrestart .eqv. .True. ) then
    if (MyRank == 0) Then
       call restart_read
    end if
    call Broad3
    !call print_ini_restart
  else
    call print_ini
    call NM_Position
    if (MyRank == 0) Then
      call Init_Velocity
      call Init_Bath
    end if
    call Broad3
    call Temp_ctr
    call nmtrans_ur2r ! x(i) = x(i) + sum_j tnm(i,j)*u(j)
    call Force_New_MPI_tk
    if ( mod(istepsv,out_step) == 0 ) call print_result_qm
    call nmtrans_fr2fur     !call Getfnm  ! fu(i) = fu(i) + sum_j fx(j)*tnm(j,i)

    if ( MyRank == 0 ) then
      call Ham_Temp
      call print_ham(Irestep)
    end if
  end if

  call Getforce_Ref
  if ( MyRank == 0 ) then

    main_loop: &
    do istepsv = Irestep+1, nstep
      select case(Ncent)
        case(0)
          continue
        case(1)
         call Nhc_Integrate_Cent
        case(3)
         call Nhc_Integrate_Cent3
      end select
      call Vupdate

      if ( Ncent == 0 ) then ! NVE simulation
        do iref=1, Nref   ! Nref = 5
          call Vupdate_Ref
          call Uupdate
          call Getforce_Ref
          call Vupdate_Ref
        end do
      else
        do iref=1, Nref   ! Nref = 5
          call Nhc_Integrate
          call Vupdate_Ref
          call Uupdate
          call Getforce_Ref
          call Vupdate_Ref
          call Nhc_Integrate
        end do
      end if

!     call Getforce_Ref
      call nmtrans_ur2r       ! x(i) = x(i) + sum_j tnm(i,j)*u(j)
      call Force_New_MPI_tk   ! Obtaining fx
      if ( mod(istepsv,out_step) == 0 ) call print_result_qm
      call nmtrans_fr2fur     ! fu(i) = fu(i) + sum_j fx(j)*tnm(j,i) !call Getfnm
      call Vupdate

      select case(Ncent)
        case(0)
          continue
        case(1)
         call Nhc_Integrate_Cent
        case(3)
         call Nhc_Integrate_Cent3
      end select
      call Ham_Temp

      if (MyRank == 0) Then
        call print_ham(istepsv)
        !call Print_Ham_tk(istepsv)
        call Restart_Write(istepsv)
      end if

      if (mod(istepsv,10) == 0) then
        call exit_program
      end if

    end do main_loop

  else
    do istepsv = Irestep+1, nstep
      call Force_New_MPI_tk
    end do
  end if

  if (MyRank == 0) then
    open(newunit=Uout,file=Fout,status='old',position='append')
      write(Uout,'(" ",a)')  repeat('*',121)
    close(Uout)
  end if

  return
end subroutine PI_NEW_MPI
