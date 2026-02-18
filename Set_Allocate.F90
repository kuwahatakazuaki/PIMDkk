subroutine Set_Allocate
  use Parameters
  implicit none

  allocate( r(Ndim,Natom,Nbead) )
  allocate( ur(Ndim,Natom,Nbead) )
  allocate( vr(Ndim,Natom,Nbead) )
  allocate( vur(Ndim,Natom,Nbead) )
  allocate( fr(Ndim,Natom,Nbead) )
  allocate( fur(Ndim,Natom,Nbead) )
  allocate( physmass(Natom))
  allocate( dnmmass(Natom, Nbead))
  allocate( fictmass(Natom, Nbead))

  if ( Isimulation == 3 ) then
    allocate( r_old(Ndim,Natom,Nbead) )
    allocate( ur_old(Ndim,Natom,Nbead) )
    allocate( vur_old(Ndim,Natom,Nbead) )
    allocate( fr_old(Ndim,Natom,Nbead) )
    allocate( fur_old(Ndim,Natom,Nbead) )
    allocate( fur_ref_old(Ndim,Natom,Nbead) )
    allocate( pot_old(Nbead))
  end if

  if ( Isimulation /= 10 ) then
    allocate( qmass(Nbead))
    allocate( fur_ref(Ndim,Natom,Nbead) )
    allocate( tnm(Nbead, Nbead))
    allocate( tnminv(Nbead, Nbead))
    allocate( u(Nbead, Nbead))
    allocate( uinv(Nbead, Nbead))
    allocate( rbath(Ndim,Natom,Nnhc,Nbead))
    allocate( vrbath(Ndim,Natom,Nnhc,Nbead))
    allocate( frbath(Ndim,Natom,Nnhc,Nbead))
  end if

  allocate( alabel(Natom))
  allocate( dipoler(Ndim,Nbead))
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
      allocate( rbc31(Ndim,Natom,Nnhc) )
      allocate( vrbc31(Ndim,Natom,Nnhc) )
      allocate( frbc31(Ndim,Natom,Nnhc) )
      allocate( qmcent31(Nnhc))
  end select

  allocate( ysweight(nys))

  return
end subroutine Set_Allocate

