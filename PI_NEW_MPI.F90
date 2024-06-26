subroutine PI_NEW_MPI
  use Parameters
  use utility, only: program_abort
  implicit none
  integer :: istep, iref, Uout

  call Setup_time_mass
  call set_pallarel
  call Normal_Mode
  call Init_Mass

  select case(Iforce)
    case(1)
      call Set_mopac
!    case(5)
!      call Set_Gamess_MPI_tk
    case(6)
      call Set_Gaussian_MPI_tk  ! Set for chk, rwf, etc.
    case(8)
      call Set_VASP
    case(9)
      call Set_siesta
    !case(10)
    !  call Set_Molcas_MPI
    case(11:16)
    case(21)
      call set_nnp_araidai
    !case default
    !  call program_abort('ERROR!!! Wrong "Iforce" option')
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

    istepsv=0

    call Force_New_MPI_tk
    call nmtrans_fr2fur     !call Getfnm  ! fu(i) = fu(i) + sum_j fx(j)*tnm(j,i)

    if ( MyRank == 0 ) then
      call Ham_Temp
      call Print_Ham_tk(nrstep)
    end if
  end if

  call Getforce_Ref
  if ( MyRank == 0 ) then

    do istep = nrstep+1, nstep
      istepsv = istep
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
        call Print_Ham_tk(istep)
        call Restart_Write(istep)
      end if
    end do

    if (mod(istep,100) == 0) then
      call exit_program
    end if

  else
    do istep = nrstep+1, nstep
      istepsv = istep
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
