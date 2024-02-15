Subroutine Print_Ham(ii,ij)

  Use Parameters
  Implicit None

  Integer :: ii,ij
  integer :: Uout

  Open(newunit=Uout,file=trim(address)//'/ham.dat',status='unknown',form='formatted',position='append')
  !Open(Uout,file=trim(address)//'/ham.dat',status='unknown',form='formatted',position='append')
    write(Uout,9998) ij,hamiltonian,potential,dkinetic,qkinetic,ebath,ebath_cent,temp,E_Virial
  close(Uout)

  !If(mod(ij,NSaveVel)==0) Then
  !   Open(igetv,file=trim(address)//'/vel.dat',status='unknown',form='formatted',position='append')
  !   Do imode=1,nbead
  !      Do iatom=1,natom
  !         write(igetv,'(3D23.15)') vux(iatom,imode),vuy(iatom,imode),vuz(iatom,imode)
  !      EndDo
  !   EndDo
  !   close(igetv)

  !   Open(igetvb,file=trim(address)//'/velb.dat',status='unknown',form='formatted',position='append')
  !   do imode = 1, NBEAD
  !      do inhc = 1, NNHC
  !         do iatom = 1, NATOM
  !            write(igetvb,'(3d23.15)') xbath(iatom,inhc,imode),ybath(iatom,inhc,imode),zbath(iatom,inhc,imode)
  !         enddo
  !      enddo
  !   enddo

  !   do imode = 1, NBEAD
  !      do inhc  = 1, NNHC
  !         do iatom = 1, NATOM
  !            write(igetvb,'(3d23.15)') vxbath(iatom,inhc,imode),vybath(iatom,inhc,imode),vzbath(iatom,inhc,imode)
  !         enddo
  !      enddo
  !   enddo

  !   If(NCent==1) Then
  !      If(NColor==1) Then
  !         do inhc=1,nnhc
  !            write(igetvb,'(2d23.15)') rbc11(inhc),vbc11(inhc)
  !         enddo
  !      Else
  !         do icolor=1,ncolor
  !            do inhc=1,nnhc
  !               write(igetvb,'(2d23.15)') rbc1(inhc,icolor),vbc1(inhc,icolor)
  !            enddo
  !         enddo
  !      EndIf
  !   EndIf

  !   If(NCent==3) Then
  !      If(NColor==1) Then

  !         do inhc=1,nnhc
  !            do iatom=1,natom
  !               write(igetvb,'(3d23.15)') xbc31(iatom,inhc),ybc31(iatom,inhc),zbc31(iatom,inhc)
  !            enddo
  !         enddo

  !         do inhc=1,nnhc
  !            do iatom=1,natom
  !               write(igetvb,'(3d23.15)') vxbc31(iatom,inhc),vybc31(iatom,inhc),vzbc31(iatom,inhc)
  !            enddo
  !         enddo

  !      Else

  !         do icolor=1,ncolor
  !            do inhc=1,nnhc
  !               do iatom=1,natom
  !                  write(igetvb,'(3d23.15)') xbc3(iatom,inhc,icolor),ybc3(iatom,inhc,icolor),zbc3(iatom,inhc,icolor)
  !               enddo
  !            enddo
  !         enddo

  !         do icolor=1,ncolor
  !            do inhc=1,nnhc
  !               do iatom=1,natom
  !                  write(igetvb,'(3d23.15)') vxbc3(iatom,inhc,icolor),vybc3(iatom,inhc,icolor),vzbc3(iatom,inhc,icolor)
  !               enddo
  !            enddo
  !         enddo

  !      EndIf
  !   EndIf
  !   close(igetvb)
  !EndIf

  If(ii==0) Then
     Write(*,*) '*********************************&
     &********************************************&
     &********************************************&
     &*******************'
     Write(*,*) '   Step      Hamiltonian        &
     Potential         DKinetic         QKinetic &
                 EBath       EBath_Cent      Temperature    E_Virial '
     Write(*,*) '*********************************&
     &*********************************************&
     &*********************************************&
     &*****************'
     Write(*,9999)  ij,hamiltonian,potential,dkinetic,qkinetic,ebath,ebath_cent,temp,E_Virial
  EndIf
  If(ii==1) Then
     If(mod(ij,100)==0) Then
        Write(*,9999)  ij,hamiltonian,potential,dkinetic,qkinetic,ebath,ebath_cent,temp,E_Virial
     EndIf
  EndIf
  If(ij==nstep) Then
     Write(*,*) '*********************************&
     &********************************************&
     &********************************************&
     &*******************'
  EndIf

 9999 format(i7,8e17.9)
 9998 format(i7,8e23.15) 
 9997 format(a137)
 Return
End Subroutine
