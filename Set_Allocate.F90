subroutine Set_Allocate
  use Parameters
  implicit none
  integer    :: th

!YK Changed for allowing colors
  allocate (  atom_num(natom))
  allocate( r(3,Natom,Nbead) )
  allocate( ur(3,Natom,Nbead) )
  allocate( vur(3,Natom,Nbead) )
  allocate( fr(3,Natom,Nbead) )
  allocate( fur(3,Natom,Nbead) )
  allocate  (physmass(natom), Stat =th)
    physmass (1:natom) = 0.d0
  allocate  (dnmmass(natom, nbead), Stat=th)
        dnmmass(1:natom,1:nbead) = 0.d0
  allocate  (fictmass(natom, nbead), Stat=th)
        fictmass(1:natom,1:nbead) = 0.d0

  if ( Isimulation /= 10 ) then
    allocate  (qmass(nbead), Stat=th)
          qmass(1:nbead) = 0.d0
    allocate( fur_ref(3,Natom,Nbead) )
    allocate (tnm(nbead, nbead), Stat=th)
       tnm(1:nbead,1:nbead) = 0.d0
    allocate (tnminv(nbead, nbead), Stat=th)
       tnminv(1:nbead,1:nbead) = 0.d0
    allocate (u(nbead, nbead), Stat=th)
       u(1:nbead,1:nbead) = 0.d0
    allocate (uinv(nbead, nbead), Stat=th)
       uinv(1:nbead,1:nbead) = 0.d0
    allocate ( rbath(3,Natom,Nnhc,Nbead))
    allocate ( vrbath(3,Natom,Nnhc,Nbead))
    allocate ( frbath(3,Natom,Nnhc,Nbead))
  end if

  allocate (alabel(natom), Stat=th)
     alabel(1:natom) = '  '
  allocate(dipoler(3,Natom))
  allocate(pressure(nbead))
  allocate (charge(natom,nbead), Stat=th)
     charge(1:natom,1:nbead) = 0.d0
  allocate (hfcc(natom,nbead), Stat=th)
     hfcc(1:natom,1:nbead) = 0.d0
  allocate (nbo(natom,nbead), Stat=th)
     nbo(1:natom,1:nbead) = 0.d0
  allocate (Eenergy(nbead), Stat=th)
     Eenergy(1:nbead) = 0.d0
  allocate (homo(nbead), Stat=th)
     homo(1:nbead) = 0.d0
  allocate (lumo(nbead), Stat=th)
     lumo(1:nbead) = 0.d0

  select case(Ncent)
    case(1)
      allocate  (rbc11(nnhc), Stat=th)
            rbc11(1:nnhc) = 0.d0
      allocate  (vbc11(nnhc), Stat=th)
            vbc11(1:nnhc) = 0.d0
      allocate  (fbc11(nnhc), Stat=th)
             fbc11(1:nnhc) = 0.d0
      allocate  (qmcent11(nnhc), Stat=th)
            qmcent11(1:nnhc) = 0.d0
    case(3)
      allocate ( rbc31(3,Natom,Nnhc) )
      allocate ( vrbc31(3,Natom,Nnhc) )
      allocate ( frbc31(3,Natom,Nnhc) )
      allocate  (qmcent31(nnhc), Stat=th)
            qmcent31(1:nnhc) = 0.d0
  end select

  allocate   (ysweight(nys), Stat = th)
    ysweight (1:nys) = 0.d0

  return
end subroutine

