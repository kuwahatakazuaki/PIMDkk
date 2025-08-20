subroutine Set_Allocate
  use Parameters
  implicit none

  !allocate( atom_num(Natom) )
  allocate( r(3,Natom,Nbead) )
  allocate( ur(3,Natom,Nbead) )
  allocate( vur(3,Natom,Nbead) )
  allocate( fr(3,Natom,Nbead) )
  allocate( fur(3,Natom,Nbead) )
  allocate( physmass(Natom))
  allocate( dnmmass(Natom, Nbead))
  allocate( fictmass(Natom, Nbead))

  if ( Isimulation == 3 ) then
    allocate( r_old(3,Natom,Nbead) )
    allocate( ur_old(3,Natom,Nbead) )
    allocate( vur_old(3,Natom,Nbead) )
    allocate( fr_old(3,Natom,Nbead) )
    allocate( fur_old(3,Natom,Nbead) )
    allocate( fur_ref_old(3,Natom,Nbead) )
    allocate( pot_old(Nbead))
  end if

  if ( Isimulation /= 10 ) then
    allocate( qmass(Nbead))
    allocate( fur_ref(3,Natom,Nbead) )
    allocate( tnm(Nbead, Nbead))
    allocate( tnminv(Nbead, Nbead))
    allocate( u(Nbead, Nbead))
    allocate( uinv(Nbead, Nbead))
    allocate( rbath(3,Natom,Nnhc,Nbead))
    allocate( vrbath(3,Natom,Nnhc,Nbead))
    allocate( frbath(3,Natom,Nnhc,Nbead))
  end if

  allocate( alabel(Natom))
  allocate( dipoler(3,Nbead))
  allocate( pressure(Nbead))
  allocate( charge(Natom,Nbead))
  allocate( hfcc(Natom,Nbead))
  allocate( nbo(Natom,Nbead))
  allocate( pot_bead(Nbead))

  select case(Ncent)
    case(1)
      allocate( rbc11(Nnhc))
      allocate( vbc11(Nnhc))
      allocate( fbc11(Nnhc))
      allocate( qmcent11(Nnhc))
    case(3)
      allocate( rbc31(3,Natom,Nnhc) )
      allocate( vrbc31(3,Natom,Nnhc) )
      allocate( frbc31(3,Natom,Nnhc) )
      allocate( qmcent31(Nnhc))
  end select

  allocate( ysweight(nys))

  return
end subroutine Set_Allocate

