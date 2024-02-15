subroutine Print_Ham_Classical(istep)
  use Parameters
  use utility, only: get_time
  implicit none
  integer :: istep, Uout

  open(newunit=Uout,file=trim(address)//'/ham.dat',status='unknown',form='formatted',position='append')
    write(Uout,9996) istep, hamiltonian, temp, potential, dkinetic, ebath_cent
    !write(Uout,9998) istep, hamiltonian, potential, dkinetic, ebath_cent, temp
  close(Uout)

  If(mod(istep,100)==0) Then
    open(newunit=Uout,file=Fout,status='old',position='append')
      write(Uout,9995)  istep, hamiltonian, temp, potential, dkinetic, ebath_cent, get_time()
      !write(*,9999)  istep,hamiltonian,potential,dkinetic,ebath_cent,temp
    close(Uout)
  EndIf

9999 format(i7,5e17.9)
9998 format(i7,5e23.15)
9997 format(a95)
9996 format(i7,5e16.8)
9995 format(i7,5e13.5,a20)
Return
End Subroutine

