subroutine print_ham(Istep)
  use Parameters
  use utility, only: get_time
  implicit none
  integer, intent(in) :: Istep
  integer :: Uout

  select case(Isimulation)
    case(0:2)
      call print_ham_qm
    case(10)
      call print_ham_cl
  end select

contains
  subroutine print_ham_qm
    if ( mod(Istep,out_step) == 0 ) then
      open(newunit=Uout,file=trim(dir_result)//'/ham.dat',status='unknown',form='formatted',position='append')
        write(Uout,1001) Istep, hamiltonian, temp, potential,dkinetic,qkinetic,ebath,ebath_cent, E_Virial
      close(Uout)
    end if

    if (mod(Istep,100)==0) Then
      open(newunit=Uout,file=Fout,status='old',position='append')
      write(Uout,1002) Istep, hamiltonian, temp, potential,dkinetic,qkinetic,ebath,ebath_cent, get_time()
      close(Uout)
    end if

    1001 format(i7,8e16.8)
    1002 format(i7,7e13.5,a20)
    return
  end subroutine print_ham_qm

  subroutine print_ham_cl
    if ( mod(Istep,out_step) == 0 ) then
      open(newunit=Uout,file=trim(dir_result)//'/ham.dat',status='unknown',form='formatted',position='append')
        write(Uout,2001) Istep, hamiltonian, temp, potential, dkinetic, ebath_cent
      close(Uout)
    end if

    If(mod(Istep,100)==0) Then
      open(newunit=Uout,file=Fout,status='old',position='append')
        write(Uout,2002)  Istep, hamiltonian, temp, potential, dkinetic, ebath_cent, get_time()
      close(Uout)
    EndIf

    2001 format(i7,5e16.8)
    2002 format(i7,5e13.5,a20)
    return
  end subroutine print_ham_cl
end subroutine print_ham



!subroutine Print_Ham_tk(Istep)
!  use Parameters
!  use utility, only: get_time
!  implicit none
!  integer, intent(in) :: Istep
!  integer :: Uout
!
!  if ( mod(Istep,out_step) == 0 ) then
!    open(newunit=Uout,file=trim(dir_result)//'/ham.dat',status='unknown',form='formatted',position='append')
!      write(Uout,1001) Istep, hamiltonian, temp, potential,dkinetic,qkinetic,ebath,ebath_cent, E_Virial
!    close(Uout)
!  end if
!
!  if (mod(Istep,100)==0) Then
!    open(newunit=Uout,file=Fout,status='old',position='append')
!    write(Uout,1002) Istep, hamiltonian, temp, potential,dkinetic,qkinetic,ebath,ebath_cent, get_time()
!    close(Uout)
!  end if
!
!!9999 format(i7,8e17.9)
!!9998 format(i7,8e23.15)
!!9997 format(a137)
!1001 format(i7,8e16.8)
!1002 format(i7,7e13.5,a20)
!return
!end subroutine Print_Ham_tk

