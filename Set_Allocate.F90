subroutine Set_Allocate
  use Parameters
  implicit none

  allocate(  atom_num(natom))
  allocate( r(3,Natom,Nbead) )
  allocate( ur(3,Natom,Nbead) )
  allocate( vur(3,Natom,Nbead) )
  allocate( fr(3,Natom,Nbead) )
  allocate( fur(3,Natom,Nbead) )
  allocate  (physmass(natom))
  allocate  (dnmmass(natom, nbead))
  allocate  (fictmass(natom, nbead))

  if ( Isimulation /= 10 ) then
    allocate  (qmass(nbead))
    allocate( fur_ref(3,Natom,Nbead) )
    allocate (tnm(nbead, nbead))
    allocate (tnminv(nbead, nbead))
    allocate (u(nbead, nbead))
    allocate (uinv(nbead, nbead))
    allocate ( rbath(3,Natom,Nnhc,Nbead))
    allocate ( vrbath(3,Natom,Nnhc,Nbead))
    allocate ( frbath(3,Natom,Nnhc,Nbead))
  end if

  allocate (alabel(natom))
  allocate(dipoler(3,Nbead))
  allocate(pressure(nbead))
  allocate (charge(natom,nbead))
  allocate (hfcc(natom,nbead))
  allocate (nbo(natom,nbead))
  allocate (Eenergy(nbead))
  allocate (homo(nbead))
  allocate (lumo(nbead))

  select case(Ncent)
    case(1)
      allocate  (rbc11(nnhc))
      allocate  (vbc11(nnhc))
      allocate  (fbc11(nnhc))
      allocate  (qmcent11(nnhc))
    case(3)
      allocate ( rbc31(3,Natom,Nnhc) )
      allocate ( vrbc31(3,Natom,Nnhc) )
      allocate ( frbc31(3,Natom,Nnhc) )
      allocate  (qmcent31(nnhc))
  end select

  allocate   (ysweight(nys))

  return
end subroutine Set_Allocate

