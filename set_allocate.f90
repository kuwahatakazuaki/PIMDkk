subroutine set_allocate
  use global_variable
  use mpi
  implicit none
  integer :: ios, ierr



  allocate( r(3,Natom, Nbead),           source=0.0d0 )
  allocate( f(3,Natom, Nbead),           source=0.0d0 )
  allocate( u(3,Natom, Nbead),           source=0.0d0 )
  allocate(vu(3,Natom, Nbead),           source=0.0d0 )
  allocate( fu(3,Natom, Nbead),          source=0.0d0 )
  allocate( fu_ref(3,Natom, Nbead),      source=0.0d0 )
  allocate( tnm(Nbead, Nbead),           source=0.0d0 )
  allocate( tnminv(Nbead, Nbead),        source=0.0d0 )
  allocate( physmass(Natom),             source=0.0d0 )
  allocate( energy(Nbead),               source=0.0d0 )
  allocate( bath(3,Natom, Nnhc, Nbead),  source=0.0d0 )
  allocate(vbath(3,Natom, Nnhc, Nbead),  source=0.0d0 )
  allocate(fbath(3,Natom, Nnhc, Nbead),  source=0.0d0 )
  allocate( dipole(4,Nbead),             source=0.0d0 )
  allocate( charge(Natom,Nbead),         source=0.0d0 )
  allocate( alabel(Natom))

!  allocate(path_result,  source=trim(path_result_temp))
!  path_result = trim(path_result)
!  allocate(path_scr,     source=trim(path_scr_temp)   )
!  path_scr = trim(path_scr)

  allocate(qmass(Nbead),                 source=0.0d0 )
  allocate(dnmmass(Natom,Nbead),         source=0.0d0 )
  allocate(fictmass(Natom,Nbead),        source=0.0d0 )
  allocate(ysweight(Nys),                source=0.0d0 )

  if (Ncent == 1) then
    if (Ncolor == 1) then

      allocate( rbc11(Nnhc),                 source=0.0d0 )
      allocate( vbc11(Nnhc),                 source=0.0d0 )
      allocate( fbc11(Nnhc),                 source=0.0d0 )
      allocate( qmcent11(Nnhc),              source=0.0d0 )

    else

      allocate( rbc1(Nnhc,Ncolor),           source=0.0d0 )
      allocate( vbc1(Nnhc,Ncolor),           source=0.0d0 )
      allocate( fbc1(Nnhc,Ncolor),           source=0.0d0 )
      allocate( qmcent1(Nnhc,Ncolor),        source=0.0d0 )

    end if

  else if (Ncent == 3) then
    if (Ncolor == 1) then

      allocate( rbc31(3,Natom,Nnhc),         source=0.0d0 )
      allocate( vbc31(3,Natom,Nnhc),         source=0.0d0 )
      allocate( fbc31(3,Natom,Nnhc),         source=0.0d0 )
      allocate( qmcent31(Nnhc),              source=0.0d0 )

    else

      allocate( rbc3(3,Natom,Nnhc,Ncolor),   source=0.0d0 )
      allocate( vbc3(3,Natom,Nnhc,Ncolor),   source=0.0d0 )
      allocate( fbc3(3,Natom,Nnhc,Ncolor),   source=0.0d0 )
      allocate( qmcent3(Nnhc,Ncolor),        source=0.0d0 )

    end if
  end if



!        allocate(path_result, source=trim(line))
end subroutine set_allocate
!  Subroutine Set_Allocate
!
!!   Use MPI
!    Use Parameters
!    Implicit None
!    Integer    :: th


!!YK Changed for allowing colors
!    Allocate (  x(natom, nbead), Stat=th)
!      If (th /= 0) Stop " Not Enough Memory "
!         x(1:natom,1:nbead) = 0.d0
!    Allocate (  y(natom, nbead), Stat=th)
!      If (th /= 0) Stop " Not Enough Memory "
!         y(1:natom,1:nbead) = 0.d0
!    Allocate (  z(natom, nbead), Stat=th)
!      If (th /= 0) Stop " Not Enough Memory "
!         z(1:natom,1:nbead) = 0.d0
!    Allocate ( fx(natom, nbead), Stat=th)
!      If (th /= 0) Stop " Not Enough Memory "
!        fx(1:natom,1:nbead) = 0.d0
!    Allocate ( fy(natom, nbead), Stat=th)
!      If (th /= 0) Stop " Not Enough Memory "
!        fy(1:natom,1:nbead) = 0.d0
!    Allocate ( fz(natom, nbead), Stat=th)
!      If (th /= 0) Stop " Not Enough Memory "
!        fz(1:natom,1:nbead) = 0.d0
!    Allocate  (physmass(natom), Stat =th)
!      If (th /= 0) Stop " Not Enough Memory "
!      physmass (1:natom) = 0.d0
!    Allocate ( ux(natom, nbead), Stat=th)
!      If (th /= 0) Stop " Not Enough Memory "
!        ux(1:natom,1:nbead) = 0.d0
!    Allocate ( uy(natom, nbead), Stat=th)
!      If (th /= 0) Stop " Not Enough Memory "
!        uy(1:natom,1:nbead) = 0.d0
!    Allocate ( uz(natom, nbead), Stat=th)
!      If (th /= 0) Stop " Not Enough Memory "
!        uz(1:natom,1:nbead) = 0.d0
!    Allocate (vux(natom, nbead), Stat=th)
!      If (th /= 0) Stop " Not Enough Memory "
!       vux(1:natom,1:nbead) = 0.d0
!    Allocate (vuy(natom, nbead), Stat=th)
!      If (th /= 0) Stop " Not Enough Memory "
!       vuy(1:natom,1:nbead) = 0.d0
!    Allocate (vuz(natom, nbead), Stat=th)
!      If (th /= 0) Stop " Not Enough Memory "
!       vuz(1:natom,1:nbead) = 0.d0
!    Allocate (fux(natom, nbead), Stat=th)
!      If (th /= 0) Stop " Not Enough Memory "
!       fux(1:natom,1:nbead) = 0.d0
!    Allocate (fuy(natom, nbead), Stat=th)
!      If (th /= 0) Stop " Not Enough Memory "
!       fuy(1:natom,1:nbead) = 0.d0
!    Allocate (fuz(natom, nbead), Stat=th)
!      If (th /= 0) Stop " Not Enough Memory "
!       fuz(1:natom,1:nbead) = 0.d0
!    Allocate (fux_ref(natom, nbead), Stat=th)
!      If (th /= 0) Stop " Not Enough Memory "
!       fux_ref(1:natom,1:nbead) = 0.d0
!    Allocate (fuy_ref(natom, nbead), Stat=th)
!      If (th /= 0) Stop " Not Enough Memory "
!       fuy_ref(1:natom,1:nbead) = 0.d0
!    Allocate (fuz_ref(natom, nbead), Stat=th)
!      If (th /= 0) Stop " Not Enough Memory "
!       fuz_ref(1:natom,1:nbead) = 0.d0
!    Allocate (tnm(nbead, nbead), Stat=th)
!      If (th /= 0) Stop " Not Enough Memory "
!       tnm(1:nbead,1:nbead) = 0.d0
!    Allocate (tnminv(nbead, nbead), Stat=th)
!      If (th /= 0) Stop " Not Enough Memory "
!       tnminv(1:nbead,1:nbead) = 0.d0
!    Allocate (u(nbead, nbead), Stat=th)
!      If (th /= 0) Stop " Not Enough Memory "
!       u(1:nbead,1:nbead) = 0.d0
!    Allocate (uinv(nbead, nbead), Stat=th)
!      If (th /= 0) Stop " Not Enough Memory "
!       uinv(1:nbead,1:nbead) = 0.d0
!    Allocate (xbath(natom, nnhc, nbead), Stat=th)
!      If (th /= 0) Stop " Not Enough Memory "
!       xbath(1:natom,1:nnhc,1:nbead) = 0.d0
!    Allocate (ybath(natom, nnhc, nbead), Stat=th)
!      If (th /= 0) Stop " Not Enough Memory "
!       ybath(1:natom,1:nnhc,1:nbead) = 0.d0
!    Allocate (zbath(natom, nnhc, nbead), Stat=th)
!      If (th /= 0) Stop " Not Enough Memory "
!       zbath(1:natom,1:nnhc,1:nbead) = 0.d0
!    Allocate (vxbath(natom, nnhc, nbead), Stat=th)
!      If (th /= 0) Stop " Not Enough Memory "
!       vxbath(1:natom,1:nnhc,1:nbead) = 0.d0
!    Allocate (vybath(natom, nnhc, nbead), Stat=th)
!      If (th /= 0) Stop " Not Enough Memory "
!       vybath(1:natom,1:nnhc,1:nbead) = 0.d0
!    Allocate (vzbath(natom, nnhc, nbead), Stat=th)
!      If (th /= 0) Stop " Not Enough Memory "
!       vzbath(1:natom,1:nnhc,1:nbead) = 0.d0
!    Allocate (fxbath(natom, nnhc, nbead), Stat=th)
!      If (th /= 0) Stop " Not Enough Memory "
!       fxbath(1:natom,1:nnhc,1:nbead) = 0.d0
!    Allocate (fybath(natom, nnhc, nbead), Stat=th)
!      If (th /= 0) Stop " Not Enough Memory "
!       fybath(1:natom,1:nnhc,1:nbead) = 0.d0
!    Allocate (fzbath(natom, nnhc, nbead), Stat=th)
!      If (th /= 0) Stop " Not Enough Memory "
!       fzbath(1:natom,1:nnhc,1:nbead) = 0.d0
!!YK Changed for reading labels
!    Allocate (alabel(natom), Stat=th)
!      If (th /= 0) Stop " Not Enough Memory "
!       alabel(1:natom) = '  '
!!YK Changed for collecting data from electronic structure calculation
!    Allocate (dipole(nbead), Stat=th)
!      If (th /= 0) Stop " Not Enough Memory "
!       dipole(1:nbead) = 0.d0
!    Allocate (dipolex(nbead), Stat=th)
!      If (th /= 0) Stop " Not Enough Memory "
!       dipolex(1:nbead) = 0.d0
!    Allocate (dipoley(nbead), Stat=th)
!      If (th /= 0) Stop " Not Enough Memory "
!       dipoley(1:nbead) = 0.d0
!    Allocate (dipolez(nbead), Stat=th)
!      If (th /= 0) Stop " Not Enough Memory "
!       dipolez(1:nbead) = 0.d0
!    Allocate (charge(natom,nbead), Stat=th)
!      If (th /= 0) Stop " Not Enough Memory "
!       charge(1:natom,1:nbead) = 0.d0
!!tkawatsu
!    Allocate (hfcc(natom,nbead), Stat=th)
!      If (th /= 0) Stop " Not Enough Memory "
!       hfcc(1:natom,1:nbead) = 0.d0
!!
!    Allocate (nbo(natom,nbead), Stat=th)
!      If (th /= 0) Stop " Not Enough Memory "
!       nbo(1:natom,1:nbead) = 0.d0
!    Allocate (Eenergy(nbead), Stat=th)
!      If (th /= 0) Stop " Not Enough Memory "
!       Eenergy(1:nbead) = 0.d0
!    Allocate (homo(nbead), Stat=th)
!      If (th /= 0) Stop " Not Enough Memory "
!       homo(1:nbead) = 0.d0
!    Allocate (lumo(nbead), Stat=th)
!      If (th /= 0) Stop " Not Enough Memory "
!       lumo(1:nbead) = 0.d0
!
!    If(NCent==1) Then
!
!       If(NColor==1) Then
!
!          Allocate  (rbc11(nnhc), Stat=th)
!             If (th /= 0) Stop " Not Enough Memory "
!                rbc11(1:nnhc) = 0.d0
!          Allocate  (vbc11(nnhc), Stat=th)
!             If (th /= 0) Stop " Not Enough Memory "
!                vbc11(1:nnhc) = 0.d0
!          Allocate  (fbc11(nnhc), Stat=th)
!             If (th /= 0) Stop " Not Enough Memory "
!                 fbc11(1:nnhc) = 0.d0
!
!       Else
!
!          Allocate  (rbc1(nnhc,ncolor), Stat=th)
!             If (th /= 0) Stop " Not Enough Memory "
!                rbc1(1:nnhc,1:ncolor) = 0.d0
!          Allocate  (vbc1(nnhc,ncolor), Stat=th)
!             If (th /= 0) Stop " Not Enough Memory "
!                vbc1(1:nnhc,1:ncolor) = 0.d0
!          Allocate  (fbc1(nnhc,ncolor), Stat=th)
!             If (th /= 0) Stop " Not Enough Memory "
!                fbc1(1:nnhc,1:ncolor) = 0.d0
!       EndIf
!
!    EndIf
!
!    If(NCent==3) Then
!
!       If(NColor==1) Then
!
!          Allocate  (xbc31(natom,nnhc), Stat=th)
!             If (th /= 0) Stop " Not Enough Memory "
!                xbc31(1:natom,1:nnhc) = 0.d0
!          Allocate  (ybc31(natom,nnhc), Stat=th)
!             If (th /= 0) Stop " Not Enough Memory "
!                ybc31(1:natom,1:nnhc) = 0.d0
!          Allocate  (zbc31(natom,nnhc), Stat=th)
!             If (th /= 0) Stop " Not Enough Memory "
!                zbc31(1:natom,1:nnhc) = 0.d0
!          Allocate  (vxbc31(natom,nnhc), Stat=th)
!             If (th /= 0) Stop " Not Enough Memory "
!                vxbc31(1:natom,1:nnhc) = 0.d0
!          Allocate  (vybc31(natom,nnhc), Stat=th)
!             If (th /= 0) Stop " Not Enough Memory "
!                vybc31(1:natom,1:nnhc) = 0.d0
!          Allocate  (vzbc31(natom,nnhc), Stat=th)
!             If (th /= 0) Stop " Not Enough Memory "
!                vzbc31(1:natom,1:nnhc) = 0.d0
!          Allocate  (fxbc31(natom,nnhc), Stat=th)
!             If (th /= 0) Stop " Not Enough Memory "
!                fxbc31(1:natom,1:nnhc) = 0.d0
!          Allocate  (fybc31(natom,nnhc), Stat=th)
!             If (th /= 0) Stop " Not Enough Memory "
!                fybc31(1:natom,1:nnhc) = 0.d0
!          Allocate  (fzbc31(natom,nnhc), Stat=th)
!             If (th /= 0) Stop " Not Enough Memory "
!                fzbc31(1:natom,1:nnhc) = 0.d0
!       Else
!
!          Allocate  (xbc3(natom,nnhc,ncolor), Stat=th)
!             If (th /= 0) Stop " Not Enough Memory "
!                xbc3(1:natom,1:nnhc,1:ncolor) = 0.d0
!          Allocate  (ybc3(natom,nnhc,ncolor), Stat=th)
!             If (th /= 0) Stop " Not Enough Memory "
!                ybc3(1:natom,1:nnhc,1:ncolor) = 0.d0
!          Allocate  (zbc3(natom,nnhc,ncolor), Stat=th)
!             If (th /= 0) Stop " Not Enough Memory "
!                zbc3(1:natom,1:nnhc,1:ncolor) = 0.d0
!          Allocate  (vxbc3(natom,nnhc,ncolor), Stat=th)
!             If (th /= 0) Stop " Not Enough Memory "
!                vxbc3(1:natom,1:nnhc,1:ncolor) = 0.d0
!          Allocate  (vybc3(natom,nnhc,ncolor), Stat=th)
!             If (th /= 0) Stop " Not Enough Memory "
!                vybc3(1:natom,1:nnhc,1:ncolor) = 0.d0
!          Allocate  (vzbc3(natom,nnhc,ncolor), Stat=th)
!             If (th /= 0) Stop " Not Enough Memory "
!                vzbc3(1:natom,1:nnhc,1:ncolor) = 0.d0
!          Allocate  (fxbc3(natom,nnhc,ncolor), Stat=th)
!             If (th /= 0) Stop " Not Enough Memory "
!                fxbc3(1:natom,1:nnhc,1:ncolor) = 0.d0
!          Allocate  (fybc3(natom,nnhc,ncolor), Stat=th)
!             If (th /= 0) Stop " Not Enough Memory "
!                fybc3(1:natom,1:nnhc,1:ncolor) = 0.d0
!          Allocate  (fzbc3(natom,nnhc,ncolor), Stat=th)
!             If (th /= 0) Stop " Not Enough Memory "
!                fzbc3(1:natom,1:nnhc,1:ncolor) = 0.d0
!       EndIf
!
!    EndIf
!!YK
!    Allocate  (dnmmass(natom, nbead), Stat=th)
!       If (th /= 0) Stop " Not Enough Memory "
!          dnmmass(1:natom,1:nbead) = 0.d0
!    Allocate  (fictmass(natom, nbead), Stat=th)
!       If (th /= 0) Stop " Not Enough Memory "
!          fictmass(1:natom,1:nbead) = 0.d0
!    Allocate  (qmass(nbead), Stat=th)
!       If (th /= 0) Stop " Not Enough Memory "
!          qmass(1:nbead) = 0.d0
!!YK added for allowing colors
!    If(NCent==1) Then
!
!       If(NColor==1) Then
!
!          Allocate  (qmcent11(nnhc), Stat=th)
!             If (th /= 0) Stop " Not Enough Memory "
!                qmcent11(1:nnhc) = 0.d0
!
!       Else
!
!          Allocate  (qmcent1(nnhc,ncolor), Stat=th)
!             If (th /= 0) Stop " Not Enough Memory "
!                qmcent1(1:nnhc,1:ncolor) = 0.d0
!
!       EndIf
!
!    EndIf
!
!    If(NCent==3) Then
!
!       If(NColor==1) Then
!
!          Allocate  (qmcent31(nnhc), Stat=th)
!             If (th /= 0) Stop " Not Enough Memory "
!                qmcent31(1:nnhc) = 0.d0
!
!       Else
!
!          Allocate  (qmcent3(nnhc,ncolor), Stat=th)
!             If (th /= 0) Stop " Not Enough Memory "
!                qmcent3(1:nnhc,1:ncolor) = 0.d0
!
!       EndIf
!
!    EndIf
!!YK
!    Allocate   (ysweight(nys), Stat = th)
!      If (th /= 0 ) Stop " Not Enough Memory "
!      ysweight (1:nys) = 0.d0
!
!    If (Order == 4) Then
!      Allocate   (hess(3*natom, 3*natom, nbead), Stat = th)
!        If (th /= 0 ) Stop " Not Enough Memory "
!        hess (1:3*natom,1:3*natom,1:nbead) = 0.d0
!      Allocate (pot_ti(nbead), Stat = th)
!        If (th /= 0 ) Stop " Not Enough Memory "
!        pot_ti(1:nbead)  = 0.d0
!      Allocate ( fx_ti(natom, nbead), Stat=th)
!        If (th /= 0) Stop " Not Enough Memory "
!          fx_ti(1:natom,1:nbead) = 0.d0
!      Allocate ( fy_ti(natom, nbead), Stat=th)
!        If (th /= 0) Stop " Not Enough Memory "
!          fy_ti(1:natom,1:nbead) = 0.d0
!      Allocate ( fz_ti(natom, nbead), Stat=th)
!        If (th /= 0) Stop " Not Enough Memory "
!          fz_ti(1:natom,1:nbead) = 0.d0
!      Allocate ( fx_org(natom, nbead), Stat=th)
!        If (th /= 0) Stop " Not Enough Memory "
!          fx_org(1:natom,1:nbead) = 0.d0
!      Allocate ( fy_org(natom, nbead), Stat=th)
!        If (th /= 0) Stop " Not Enough Memory "
!          fy_org(1:natom,1:nbead) = 0.d0
!      Allocate ( fz_org(natom, nbead), Stat=th)
!        If (th /= 0) Stop " Not Enough Memory "
!          fz_org(1:natom,1:nbead) = 0.d0
!    EndIf
!    If(NForce==3) Then
!      Allocate(no_atom(natom))
!      no_atom(1:natom)=0
!    EndIf
!
!    Return
!
!  End Subroutine
