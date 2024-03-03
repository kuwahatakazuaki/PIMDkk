subroutine Print_Ham_tk(ij)
  use Parameters
  use utility, only: get_time
  implicit none
  integer :: ij, Uham,Uout

  open(newunit=Uham,file=trim(address)//'/ham.dat',status='unknown',form='formatted',position='append')
    write(Uham,9996) ij, hamiltonian, temp, potential,dkinetic,qkinetic,ebath,ebath_cent, E_Virial
    !write(Uham,9998) ij,hamiltonian,potential,dkinetic,qkinetic,ebath,ebath_cent,temp,E_Virial
  close(Uham)


  If(mod(ij,100)==0) Then
    open(newunit=Uout,file=Fout,status='old',position='append')
    write(Uout,9995) ij, hamiltonian, temp, potential,dkinetic,qkinetic,ebath,ebath_cent, get_time()
    !Write(*,9999)  ij,hamiltonian,potential,dkinetic,qkinetic,ebath,ebath_cent,temp,E_Virial
    close(Uout)
  EndIf

9999 format(i7,8e17.9)
9998 format(i7,8e23.15)
9997 format(a137)
9996 format(i7,8e16.8)
9995 format(i7,7e13.5,a20)
Return
End Subroutine

