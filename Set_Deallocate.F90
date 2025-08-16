Subroutine Set_Deallocate
  use Parameters
  implicit none

  deallocate(physmass)
  deallocate(  r)
  deallocate( ur)
  deallocate(vur)
  deallocate( fr)
  deallocate(fur)
  deallocate(alabel)
  deallocate(dipoler)
  deallocate(charge)
  deallocate(hfcc)
  deallocate(nbo)
  deallocate(pot_bead)

  if ( Isimulation /= 10 ) then
    deallocate(fur_ref)
    deallocate(tnm)
    deallocate(tnminv)
    deallocate(u)
    deallocate(uinv)
    deallocate(rbath)
    deallocate(vrbath)
    deallocate(frbath)
    deallocate(qmass)
  end if

  select case(Ncent)
    case(1)
      deallocate(rbc11)
      deallocate(vbc11)
      deallocate(fbc11)
      deallocate(qmcent11)
    case(3)
      deallocate(rbc31)
      deallocate(vrbc31)
      deallocate(frbc31)
      deallocate(qmcent31)
  end select
  deallocate(dnmmass)
  deallocate(fictmass)
  deallocate(ysweight)

  return
end subroutine Set_Deallocate


