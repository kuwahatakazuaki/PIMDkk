subroutine Set_Allocate
  use Parameters
  implicit none
  integer    :: th

  allocate(  atom_num(natom))
  allocate( r(3,Natom,Nbead) )
  allocate( ur(3,Natom,Nbead) )
  allocate( vur(3,Natom,Nbead) )
  allocate( fr(3,Natom,Nbead) )
  allocate( fur(3,Natom,Nbead) )
  allocate  (physmass(natom))
    !physmass (1:natom) = 0.d0
  allocate  (dnmmass(natom, nbead))
        !dnmmass(1:natom,1:nbead) = 0.d0
  allocate  (fictmass(natom, nbead))
        !fictmass(1:natom,1:nbead) = 0.d0

  if ( Isimulation /= 10 ) then
    allocate  (qmass(nbead))
          !qmass(1:nbead) = 0.d0
    allocate( fur_ref(3,Natom,Nbead) )
    allocate (tnm(nbead, nbead))
       !tnm(1:nbead,1:nbead) = 0.d0
    allocate (tnminv(nbead, nbead))
       !tnminv(1:nbead,1:nbead) = 0.d0
    allocate (u(nbead, nbead))
       !u(1:nbead,1:nbead) = 0.d0
    allocate (uinv(nbead, nbead))
       !uinv(1:nbead,1:nbead) = 0.d0
    allocate ( rbath(3,Natom,Nnhc,Nbead))
    allocate ( vrbath(3,Natom,Nnhc,Nbead))
    allocate ( frbath(3,Natom,Nnhc,Nbead))
  end if

  allocate (alabel(natom))
     !alabel(1:natom) = '  '
  allocate(dipoler(3,Nbead))
  allocate(pressure(nbead))
  allocate (charge(natom,nbead))
     !charge(1:natom,1:nbead) = 0.d0
  allocate (hfcc(natom,nbead))
     !hfcc(1:natom,1:nbead) = 0.d0
  allocate (nbo(natom,nbead))
     !nbo(1:natom,1:nbead) = 0.d0
  allocate (Eenergy(nbead))
     !Eenergy(1:nbead) = 0.d0
  allocate (homo(nbead))
     !homo(1:nbead) = 0.d0
  allocate (lumo(nbead))
     !lumo(1:nbead) = 0.d0

  select case(Ncent)
    case(1)
      allocate  (rbc11(nnhc))
            !rbc11(1:nnhc) = 0.d0
      allocate  (vbc11(nnhc))
            !vbc11(1:nnhc) = 0.d0
      allocate  (fbc11(nnhc))
            ! fbc11(1:nnhc) = 0.d0
      allocate  (qmcent11(nnhc))
            !qmcent11(1:nnhc) = 0.d0
    case(3)
      allocate ( rbc31(3,Natom,Nnhc) )
      allocate ( vrbc31(3,Natom,Nnhc) )
      allocate ( frbc31(3,Natom,Nnhc) )
      allocate  (qmcent31(nnhc))
            !qmcent31(1:nnhc) = 0.d0
  end select

  allocate   (ysweight(nys))
    !ysweight (1:nys) = 0.d0

  return
end subroutine Set_Allocate

