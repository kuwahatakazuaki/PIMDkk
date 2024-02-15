subroutine Set_Allocate_Classical
  use Parameters
  implicit none
  integer    :: th


!YK Changed for allowing colors
  Allocate (  atom_num(natom))
  allocate( r(3,Natom,Nbead) )
  allocate( ur(3,Natom,Nbead) )
  allocate( vur(3,Natom,Nbead) )
  allocate( fr(3,Natom,Nbead) )
  allocate( fur(3,Natom,Nbead) )
  Allocate  (physmass(natom), Stat =th)
    physmass (1:natom) = 0.d0
  Allocate  (dnmmass(natom, nbead), Stat=th)
        dnmmass(1:natom,1:nbead) = 0.d0
  Allocate  (fictmass(natom, nbead), Stat=th)
        fictmass(1:natom,1:nbead) = 0.d0

  Allocate (alabel(natom), Stat=th)
     alabel(1:natom) = '  '
  allocate(dipoler(3,Natom))
  allocate(pressure(nbead))
  Allocate (charge(natom,nbead), Stat=th)
     charge(1:natom,1:nbead) = 0.d0
  Allocate (hfcc(natom,nbead), Stat=th)
     hfcc(1:natom,1:nbead) = 0.d0
  Allocate (nbo(natom,nbead), Stat=th)
     nbo(1:natom,1:nbead) = 0.d0
  Allocate (Eenergy(nbead), Stat=th)
     Eenergy(1:nbead) = 0.d0
  Allocate (homo(nbead), Stat=th)
     homo(1:nbead) = 0.d0
  Allocate (lumo(nbead), Stat=th)
     lumo(1:nbead) = 0.d0

  If(NCent==3) Then
    allocate ( rbc31(3,Natom,Nnhc) )
    allocate ( vrbc31(3,Natom,Nnhc) )
    allocate ( frbc31(3,Natom,Nnhc) )
  EndIf
  If(NCent==1) Then
    Allocate  (rbc11(nnhc), Stat=th)
          rbc11(1:nnhc) = 0.d0
    Allocate  (vbc11(nnhc), Stat=th)
          vbc11(1:nnhc) = 0.d0
    Allocate  (fbc11(nnhc), Stat=th)
           fbc11(1:nnhc) = 0.d0
  EndIf

  If(NCent==1) Then
        Allocate  (qmcent11(nnhc), Stat=th)
              qmcent11(1:nnhc) = 0.d0
  EndIf
  If(NCent==3) Then
        Allocate  (qmcent31(nnhc), Stat=th)
              qmcent31(1:nnhc) = 0.d0
  EndIf

  Allocate   (ysweight(nys), Stat = th)
    ysweight (1:nys) = 0.d0

  Return
End Subroutine

