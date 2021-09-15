Subroutine Print_Ham_tk(ij)

Use Parameters
Implicit None

Integer  :: ij

Open(iham,file=trim(address)//'/ham.dat',status='unknown',form='formatted',position='append')
  write(iham,9998) ij,hamiltonian,potential,dkinetic,qkinetic,ebath,ebath_cent,temp,E_Virial
close(iham)


 If(mod(ij,100)==0) Then
     Write(*,9999)  ij,hamiltonian,potential,dkinetic,qkinetic,ebath,ebath_cent,temp,E_Virial
 EndIf

 9999 format(i7,8e17.9)
 9998 format(i7,8e23.15)
 9997 format(a137)
    Return
End Subroutine
