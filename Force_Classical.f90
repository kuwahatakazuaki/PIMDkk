Subroutine Force_Classical

  Use Parameters

  Implicit None

  Do iatom=1,natom
     x(iatom,1)=ux(iatom,1)
     y(iatom,1)=uy(iatom,1)
     z(iatom,1)=uz(iatom,1)
  EndDo

!    If(NForce==1) Call Force_MOPAC
    If(NForce==2) goto 900 ! call Force_Gaussian
    If(NForce==6) Call Force_Gaussian_classical
!    If(NForce==3) Call Force_DFTB
!    If(NForce==4) Call Force_Turbomole
!    If(NForce==7) Call Force_Turbomole_tk

  Do iatom=1,natom
     fux(iatom,1)=fx(iatom,1)
     fuy(iatom,1)=fy(iatom,1)
     fuz(iatom,1)=fz(iatom,1)
  EndDo

  Return
900 continue
print *, 'We no longer use "Nforce = 2"'
print *, 'Please use "Nforce = 6"'
stop
End Subroutine

