subroutine Classical
  use Parameters
  use utility, only: program_abort
  implicit none
  integer :: istep, iref, Uout

  Call Setup_Classical

  select case(Iforce)
    case(1)
      call Set_mopac
    case(6)
      call Set_Gaussian_MPI_tk
    case(8)
      call Set_VASP
    case(11:16)
!      print *, "Do I need something prepare??"
    case(21)
      call set_nnp_araidai
  end select


  if ( Lrestart .eqv. .True. ) then
    call restart_read_classical
  else
    call print_ini
    !call print_ini_cl
    call Init_Velocity
    call Init_Bath_Classical
    call Temp_ctr
    istepsv=0
    call Force_Classical
    call Ham_Temp_Classical
    call Print_Ham_Classical(nrstep)
  end if

  !Do istep=nrstep+1, nstep
  do istep = nrstep+1, nstep
    istepsv=istep
    select case(Ncent)
      case(0)
        continue
      case(1)
       call Nhc_Integrate_Cent
      case(3)
       call Nhc_Integrate_Cent3
    end select
    call Vupdate
    call Uupdate
    call Force_Classical
    call Vupdate
    select case(Ncent)
      case(0)
        continue
      case(1)
       call Nhc_Integrate_Cent
      case(3)
       call Nhc_Integrate_Cent3
    end select
    call Ham_Temp_Classical
    call Print_Ham_Classical(istep)
    call Restart_Write_Classical(istep)

    if (mod(istep,100) == 0) then
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
