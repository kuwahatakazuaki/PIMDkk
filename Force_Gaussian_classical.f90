Subroutine Force_Gaussian_classical
!Subroutine Force_Gaussian

!$  use omp_lib
  Use Parameters
  use Parameter_tk, only: laddress, addresstmp, address0
  Implicit None

  Character    (Len=90)              :: key1, key2, key3, key4, key5, line
  Integer                            :: iline, numeach, ista, iend, itemp,id,imode2,iatom2
  Double Precision                   :: enetemp
  integer                            :: chr, chl
!    Logical                            :: lopen

  if (trim(version_gaussian) == 'g16') then
    if (theory == 0) then
      key1  = ('SCF Done')
    else if (theory == 1) then
      key1  = ('EUMP2')
    end if
!    key1  = ('SCF Done')
    key2  = ('Number     Number              X              Y              Z')
    key3  = ('Dipole moment')
    key4  = ('Mulliken charges')
    key5  = ('Isotropic Fermi Contact Couplings')
  else
    stop 'ERROR!!! wrog keyword of "version_gaussian"'
  end if

  id=0

!$omp parallel private(line,iatom2,id,iline,enetemp)
!$omp do
    do imode2=1,nbead
!$     id=omp_get_thread_num()

!  Kuwahata 2019/11/05
!    laddress=len_trim(address)+1
    write(addresstmp(laddress+1:laddress+6),'(i5.5,A1)') imode2,'/'
!  End Kuwahata 2019/11/05

     call system('mkdir -p '//trim(addresstmp))
     call system('cp gauss.tmp '//trim(addresstmp))

     If(NGenGau==1) Then
        call system('cp gauss.bss '//trim(addresstmp))
     EndIf

     open(igauss+id,file=trim(addresstmp)//'/gauss.tmp1',status='unknown')
       write(igauss+id,'(a)') '%Chk='//trim(addresstmp)//'/gauss.chk'
       write(igauss+id,'(a)') '%RWF='//trim(addresstmp)
       write(igauss+id,'(a)') '%Int='//trim(addresstmp)
       write(igauss+id,'(a)') '%D2E='//trim(addresstmp)
     close(igauss+id)

     open(igauss+id,file=trim(addresstmp)//'gauss.xyz',status='unknown')
       do iatom2=1,natom
!          write(igauss+id,9999) alabel(iatom2),x(iatom2,imode2)/bohr, &
!                             y(iatom2,imode2)/bohr,z(iatom2,imode2)/bohr
         write(igauss+id,9999) alabel(iatom2),x(iatom2,imode2)*bohr_inv, &
                             y(iatom2,imode2)*bohr_inv,z(iatom2,imode2)*bohr_inv
       enddo
       write(igauss+id,*)
     close(igauss+id)

     call system('cat '//trim(addresstmp)//'gauss.tmp1 > '//trim(addresstmp)//'gauss.com') ! location of Check file
     call system('cat '//trim(addresstmp)//'gauss.tmp >> '//trim(addresstmp)//'gauss.com') ! Setting of Gausian
     call system('cat '//trim(addresstmp)//'gauss.xyz >> '//trim(addresstmp)//'gauss.com')

     If(NGenGau==1) Then
        call system('cat '//trim(addresstmp)//'gauss.bss >> '//trim(addresstmp)//'gauss.com')
     EndIf

     call system(trim(address0)//'g0xrun_p '//trim(addresstmp)//'gauss.com '//trim(addresstmp)//'gauss.log '//trim(addresstmp))
     open(igauss+id,file=trim(addresstmp)//'gauss.log')

!  +++ Reading "SCF Done" +++
     do
       read(igauss+id,'(a)',end=401) line
       iline=index(line,trim(key1))
!       iline=index(trim(line),trim(key1))
       if(iline > 0) exit
     end do

     if (theory == 0) then
       chl = index(line, '=')
       chr = index(line(chl+1:len(line)-chl), '.') + chl
       chr = index(line(chr+1:len(line)-chr), ' ') + chr
       read(line(chl+1:chr-1),*) enetemp  ! For DFT
     else if (theory == 1) then
       read(line(38:60),*) enetemp  ! For MP2
     end if

     Eenergy(imode2)=enetemp
!  +++ End Reading "SCF Done" +++


!  +++ Reading "Atomic Force" +++
     Rewind(igauss+id)
     do
       Read(igauss+id,'(a)',end=402) line
       iline=index(trim(line),trim(key2))
       if(iline > 0) exit
     end do
     read(igauss+id,'()')
     do iatom2=1,natom
        read(igauss+id,'(a)') line
        read(line(24:38),*) fx(iatom2,imode2)
        read(line(39:53),*) fy(iatom2,imode2)
        read(line(54:68),*) fz(iatom2,imode2)
     enddo
!  +++ End Reading "Atomic Force" +++

!  +++ Reading "Dipole moment" +++
     if(nodipole==0) then
       rewind(igauss+id)
       do
         Read(igauss+id,'(a)',end=403) line
         iline=index(trim(line),trim(key3))
         if(iline > 0) exit
       end do
       Read(igauss+id,'(a)') line
       Read(line(20:26),*) dipolex(imode2)
       Read(line(46:52),*) dipoley(imode2)
       Read(line(72:78),*) dipolez(imode2)
     endif
!  +++ End Reading "Dipole moment" +++

!  +++ Reading "Mulliken charge" +++
     if(nocharge==0) then
       rewind(igauss+id)
       do
         Read(igauss+id,'(a)',end=404) line
         iline=index(trim(line),trim(key4))
         if(iline > 0) exit
       end do
       read(igauss+id,'()')
       do iatom2=1,natom
          Read(igauss+id,'(a)') line
          Read(line(11:21),*) charge(iatom2,imode2)
       enddo

       if(nohfcc==0) then
         do
           Read(igauss+id,'(a)',end=405) line
           iline=index(trim(line),trim(key5))
           if(iline > 0) exit
         end do
         read(igauss+id,'()')
         do iatom2=1,natom
           Read(igauss+id,'(a)') line 
           Read(line(36:49),*) hfcc(iatom2,imode2) 
         enddo
       endif
     endif
!  +++ End Reading "Mulliken charge" +++

     close(igauss+id)

     call system('rm -rf '//trim(addresstmp)//'Gau*')
    enddo
!$omp end do
!$omp end parallel

  if(nodipole==0) then
     Open(igetd,file=trim(address)//'/dipole.dat',status='unknown',form='formatted',position='append')
    write(igetd,'(I10)') istepsv
     DO imode=1,nbead
!        write(igetd,9995) dipole(imode),dipolex(imode),dipoley(imode),dipolez(imode)
       write(igetd,8008) dipolex(imode),dipoley(imode),dipolez(imode)
     ENDDO
     Close(igetd)
  endif

  if(nocharge==0) then
     open(igetc,file=trim(address)//'/charge.dat',status='unknown',form='formatted',position='append')

      write(igetc,'(I10)') istepsv
      do imode=1,nbead
        write(igetc,8007) charge(:,imode)
      end do
     Close(igetc)
  endif

!  Open(igetx,file=trim(address)//'/coor.dat',status='unknown',form='formatted',position='append')
  open(igetx,file=trim(address)//'/coor.xyz',status='unknown',form='formatted',position='append')
    write(igetx,'(I5)') natom*nbead
    write(igetx,'(I10)') istepsv
    DO imode=1,nbead
      DO iatom=1,natom
        write(igetx,9999) alabel(iatom),x(iatom,imode)/bohr,y(iatom,imode)/bohr,z(iatom,imode)/bohr
      ENDDO
    ENDDO
  Close(igetx)

! Kuwahata 2019/11/21 skipping writing cent.xyz
!  Open(igetxyz,file=trim(address)//'/cent.xyz',status='unknown',form='formatted',position='append')
!    write(igetxyz,'(I5)') natom
!    write(igetxyz,'(I10)') istepsv
!    DO iatom=1,natom
!       write(igetxyz,9999) alabel(iatom),ux(iatom,1)/bohr,uy(iatom,1)/bohr,uz(iatom,1)/bohr
!    ENDDO
!  Close(igetxyz)
! End Kuwahata 2019/11/21 skipping writing cent.xyz

  if (Save_force .eqv. .True.) then
    open(igetf,file=trim(address)//'/force.dat',status='unknown',form='formatted',position='append')
      do imode=1,nbead
        do iatom=1,natom
          write(igetf,9998) fx(iatom,imode),fy(iatom,imode),fz(iatom,imode)
        end do
      end do
    close(igetf)
  end if


  Open(igete,file=trim(address)//'/ene.dat',status='unknown',form='formatted',position='append')
    write(igete,'(I10)') istepsv
    DO imode=1,nbead
       write(igete,8006) Eenergy(imode)
    ENDDO
  Close(igete)

!  Open(igethl,file='homolumo.dat',status='unknown',form='formatted',position='append')
!  DO imode=1,nbead
!  write(igethl,9997) homo(imode),lumo(imode)
!  ENDDO
!  Close(igethl)

! Kuwahata 2020/01/11
!  do imode=1,nbead
!     do iatom=1,natom
!        fx(iatom,imode)=fx(iatom,imode)/dp
!        fy(iatom,imode)=fy(iatom,imode)/dp
!        fz(iatom,imode)=fz(iatom,imode)/dp
!     enddo
!  enddo
  fx(:,:)=fx(:,:)*dp_inv
  fy(:,:)=fy(:,:)*dp_inv
  fz(:,:)=fz(:,:)*dp_inv
! End Kuwahata 2020/01/11

  potential=0.0D+00
  do imode=1,nbead
     potential=potential+Eenergy(imode)
  enddo
  potential=potential/dp

9999 format(a2,1x,E15.9,1x,E15.9,1x,E15.9)
9998 format(3E23.15)
9997 format(2E23.15)
9996 format(E23.15)
8006 format(F0.10)     ! Potential
8007 format(100F10.6)  ! Charge
8008 format(3F10.5)    ! Dipole
9995 format(4E23.15)

  Return
401 print *, 'ERROR!!: We can not find "SCF Done" in Gaussian output'; stop
402 print *, 'ERROR!!: We can not find "Force" in Gaussian output'; stop
403 print *, 'ERROR!!: We can not find "Dipole moment" in Gaussian output'; stop
404 print *, 'ERROR!!: We can not find "Mulliken charges" in Gaussian output'; stop
405 print *, 'ERROR!!: We can not find "Fermi Contact" in Gaussian output'; stop
End Subroutine

!    key1  = (' SCF Done')
!    key2  = (' Number     Number              X              Y              Z')
!    key3  = (' Dipole        =')
!    key4  = ('Mulliken charges')

