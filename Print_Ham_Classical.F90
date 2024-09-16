subroutine Print_Ham_Classical(Istep)
  use Parameters
  use utility, only: get_time
  implicit none
  integer, intent(in) :: Istep
  integer :: Uout

  if ( mod(Istep,out_step) == 0 ) then
    open(newunit=Uout,file=trim(address)//'/ham.dat',status='unknown',form='formatted',position='append')
      write(Uout,2001) Istep, hamiltonian, temp, potential, dkinetic, ebath_cent
    close(Uout)
  end if

  If(mod(Istep,100)==0) Then
    open(newunit=Uout,file=Fout,status='old',position='append')
      write(Uout,2002)  Istep, hamiltonian, temp, potential, dkinetic, ebath_cent, get_time()
    close(Uout)
  EndIf

!9999 format(i7,5e17.9)
!9998 format(i7,5e23.15)
!9997 format(a95)
2001 format(i7,5e16.8)
2002 format(i7,5e13.5,a20)
Return
End Subroutine

