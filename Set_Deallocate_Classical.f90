  Subroutine Set_Deallocate_Classical

    Use Parameters
    Implicit None

    Deallocate  (physmass)
    Deallocate ( ux)
    Deallocate ( uy)
    Deallocate ( uz)
    Deallocate ( x)
    Deallocate ( y)
    Deallocate ( z)
    Deallocate (vux)
    Deallocate (vuy)
    Deallocate (vuz)
    Deallocate (fux)
    Deallocate (fuy)
    Deallocate (fuz)
    Deallocate (fx)
    Deallocate (fy)
    Deallocate (fz)
    Deallocate (alabel)
    Deallocate (dipole)
    Deallocate (dipolex)
    Deallocate (dipoley)
    Deallocate (dipolez)
    Deallocate (charge)
    Deallocate (nbo)
    Deallocate (Eenergy)
    Deallocate (homo)
    Deallocate (lumo)
!YK Changed for allowing colors
    If(NCent==1) Then
       If(NColor==1) Then
          Deallocate  (rbc11)
          Deallocate  (vbc11)
          Deallocate  (fbc11)
       Else
          Deallocate  (rbc1)
          Deallocate  (vbc1)
          Deallocate  (fbc1)
       EndIf
    EndIf
    If(NCent==3) Then
       If(NColor==1) Then
          Deallocate  (xbc31)
          Deallocate  (ybc31)
          Deallocate  (zbc31)
          Deallocate  (vxbc31)
          Deallocate  (vybc31)
          Deallocate  (vzbc31)
          Deallocate  (fxbc31)
          Deallocate  (fybc31)
          Deallocate  (fzbc31)
       Else
          Deallocate  (xbc3)
          Deallocate  (ybc3)
          Deallocate  (zbc3)
          Deallocate  (vxbc3)
          Deallocate  (vybc3)
          Deallocate  (vzbc3)
          Deallocate  (fxbc3)
          Deallocate  (fybc3)
          Deallocate  (fzbc3)
       EndIf
    EndIf
!YK
    Deallocate  (dnmmass)
    Deallocate  (fictmass)
!YK added for allowing colors
    If(NCent==1) Then
       If(NColor==1) Then
          Deallocate  (qmcent11)
       Else
          Deallocate  (qmcent1)
       EndIf
    EndIf
    If(NCent==3) Then
       If(NColor==1) Then
          Deallocate  (qmcent31)
       Else
          Deallocate  (qmcent3)
       EndIf
    EndIf
!YK
    Deallocate   (ysweight)
    If (Order == 4)Then
      Deallocate   (hess)
      Deallocate (pot_ti)
      Deallocate ( fx_ti)
      Deallocate ( fy_ti)
      Deallocate ( fz_ti)
      Deallocate ( fx_org)
      Deallocate ( fy_org)
      Deallocate ( fz_org)
    EndIf
    If(NForce==3) Then
      Deallocate(no_atom)
    EndIf

    Return

  End Subroutine



