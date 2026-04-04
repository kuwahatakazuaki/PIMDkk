subroutine Classical
  use Parameters
  use utility, only: program_abort
  implicit none
  integer :: iref, Uout

  call Setup_time_mass
  call Init_Mass
  call set_Iforce

  if ( Lrestart .eqv. .True. ) then
    !call restart_read_classical
    call restart_read
  else
    call print_ini
    call Init_Velocity
    call Init_Bath
    call Temp_ctr
    call Force_Classical
    if ( mod(istepsv,out_step) == 0 ) call print_result
    call Ham_Temp_Classical
    call print_ham(Irestep)
  end if

  do istepsv = Irestep + 1, Nstep
    select case(Ncent)
      case(0)
        continue
      case(1)
        call Nhc_Integrate_Cent
      case(3)
        call Nhc_Integrate_Cent3
    end select
    call update_vel_nor
    call update_pos_nor
    call Force_Classical
    if ( mod(istepsv,out_step) == 0 ) call print_result
    call update_vel_nor
    select case(Ncent)
      case(0)
        continue
      case(1)
        call Nhc_Integrate_Cent
      case(3)
        call Nhc_Integrate_Cent3
    end select
    call Ham_Temp_Classical
    call print_ham(istepsv)
    !if ( mod(istepsv,out_step)==0 ) call restart_write_Classical(istepsv)
    if ( mod(istepsv,out_step)==0 ) call restart_write(istepsv)

    if (mod(istepsv,10) == 0) then
      call exit_program
    end if
  end do

  if (MyRank == 0) then
    open(newunit=Uout,file=Fout,status='old',position='append')
      write(Uout,'(" ",a)') repeat('*',95)
    close(Uout)
  end if

return
end subroutine Classical
