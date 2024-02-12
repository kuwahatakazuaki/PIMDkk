Subroutine Set_Deallocate
  Use Parameters
  Implicit None

  deallocate  (physmass)
  deallocate(  r)
  deallocate( ur)
  deallocate(vur)
  deallocate( fr)
  deallocate(fur)
  deallocate (alabel)
  deallocate (dipoler)
  deallocate (charge)
  deallocate (hfcc)
  deallocate (nbo)
  deallocate (Eenergy)
  deallocate (homo)
  deallocate (lumo)

  if ( Isimulation /= 10 ) then
    deallocate(fur_ref)
    deallocate (tnm)
    deallocate (tnminv)
    deallocate (u)
    deallocate (uinv)
    deallocate(rbath)
    deallocate(vrbath)
    deallocate(frbath)
    deallocate  (qmass)
  end if

  select case(Ncent)
    case(1)
    !If(NCent==1) Then
        deallocate  (rbc11)
        deallocate  (vbc11)
        deallocate  (fbc11)
        deallocate  (qmcent11)
    !EndIf
    !If(NCent==3) Then
    case(3)
        deallocate(rbc31)
        deallocate(vrbc31)
        deallocate(frbc31)
        deallocate  (qmcent31)
    !EndIf
  end select
  deallocate  (dnmmass)
  deallocate  (fictmass)
!YK added for allowing colors
  deallocate   (ysweight)

  Return
End Subroutine Set_Deallocate


