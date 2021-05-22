Subroutine Print_Ham_Classical(istep)

Use Parameters
Implicit None

Integer                            :: istep

If(mod(istep,NSaveVel)==0) Then
   Open(igetv,file=trim(address)//'/vel.dat',status='unknown',form='formatted',position='append')
   Do imode=1,nbead
      Do iatom=1,natom
         write(igetv,'(3D23.15)') vux(iatom,1),vuy(iatom,1),vuz(iatom,1)
      EndDo
   EndDo
   close(igetv)
   Open(igetvb,file=trim(address)//'/velb.dat',status='unknown',form='formatted',position='append')

   If(NCent==1) Then
      If(NColor==1) Then
         do inhc=1,nnhc
            write(igetvb,'(2d23.15)') rbc11(inhc),vbc11(inhc)
         enddo
      Else
         do icolor=1,ncolor
            do inhc=1,nnhc
               write(igetvb,'(2d23.15)') rbc1(inhc,icolor),vbc1(inhc,icolor)
            enddo
         enddo
      EndIf
   EndIf

   If(NCent==3) Then
      If(NColor==1) Then

         do inhc=1,nnhc
            do iatom=1,natom
               write(igetvb,'(3d23.15)') xbc31(iatom,inhc),ybc31(iatom,inhc),zbc31(iatom,inhc)
            enddo
         enddo

         do inhc=1,nnhc
            do iatom=1,natom
               write(igetvb,'(3d23.15)') vxbc31(iatom,inhc),vybc31(iatom,inhc),vzbc31(iatom,inhc)
            enddo
         enddo

      Else

         do icolor=1,ncolor
            do inhc=1,nnhc
               do iatom=1,natom
                  write(igetvb,'(3d23.15)') xbc3(iatom,inhc,icolor),ybc3(iatom,inhc,icolor),zbc3(iatom,inhc,icolor)
               enddo
            enddo
         enddo

         do icolor=1,ncolor
            do inhc=1,nnhc
               do iatom=1,natom
                  write(igetvb,'(3d23.15)') vxbc3(iatom,inhc,icolor),vybc3(iatom,inhc,icolor),vzbc3(iatom,inhc,icolor)
               enddo
            enddo
         enddo

      EndIf
   EndIf
   close(igetvb)
EndIf

If(mod(istep,100)==0) Then
   Write(*,9999)  istep,hamiltonian,potential,dkinetic,ebath_cent,temp
EndIf

Open(iham,file=trim(address)//'/ham.dat',status='unknown',form='formatted',position='append')
  write(iham,9998) istep,hamiltonian,potential,dkinetic,ebath_cent,temp
close(iham)

If(istep==nstep) Then
  Write(*,'(a)') ' ************************************************************************************************'
EndIf

9999 format(i7,5e17.9)
9998 format(i7,5e23.15) 
9997 format(a95)
Return
End Subroutine
