Subroutine Force_Gaussian_classical
  Use Parameters
  use utility, only: program_abort
  Implicit None

  Character (Len=90) :: key1, key2, key3, key4, key5, key6
  character(len=120) :: line
  Integer            :: iline, id,imode2,iatom2
  Double Precision   :: enetemp
  !integer            :: chr, chl
  Integer            :: index1
  integer :: igauss = 20

  key1  = ('SCF Done')
  key6  = ('EUMP2')

  key2  = ('Number     Number              X              Y              Z')
  key3  = ('Dipole moment')
  key5  = ('Isotropic Fermi Contact Couplings')
  if (trim(version) == 'g16') then
    key4  = ('Mulliken charges')
  else if ( trim(version) == 'g09' ) then
    key4  = ('Mulliken atomic charges')
  else
    print *, '"version" is ', version
    call program_abort('ERROR!!! wrong keyword of "version"')
    !stop 'ERROR!!! wrong keyword of "version"'
  end if

  id=0

    do imode2=1,nbead

!  Kuwahata 2019/11/05
    write(addresstmp(laddress+1:laddress+6),'(i5.5,A1)') imode2,'/'
!  End Kuwahata 2019/11/05

     call system('cp gauss.tmp '//trim(addresstmp))

     !If(NGenGau==1) Then
     !   call system('cp gauss.bss '//trim(addresstmp))
     !EndIf

     call system('cat '//trim(addresstmp)//'gauss.tmp1 > '//trim(addresstmp)//'gauss.com') ! location of Check file
     open(igauss+id,file=trim(addresstmp)//'gauss.com',status='old',position='append')
       do iatom2=1,natom
         write(igauss+id,9999) alabel(iatom2), r(:,iatom2,imode2)*bohr_inv
         !write(igauss+id,9999) alabel(iatom2),x(iatom2,imode2)*bohr_inv, y(iatom2,imode2)*bohr_inv,z(iatom2,imode2)*bohr_inv
       enddo
       write(igauss+id,*)
     close(igauss+id)
!     call system('cat '//trim(addresstmp)//'gauss.xyz >> '//trim(addresstmp)//'gauss.com')

     !If(NGenGau==1) Then
     !   call system('cat '//trim(addresstmp)//'gauss.bss >> '//trim(addresstmp)//'gauss.com')
     !EndIf

! Udagawa Start 2021.05.24 --->
     !If(istepsv == 0 .OR. (nRestart ==1 .AND. istepsv == nrstep+1)) then
     If(istepsv == 0 .OR. ( ( Lrestart .eqv. .True. ) .AND. istepsv == nrstep+1)) then
       call system ('sed -e "s/[Gg][Uu][Ee][Ss][Ss]=[Rr][Ee][Aa][Dd]//g" '//trim(addresstmp)//'gauss.com &
                     & > '//trim(addresstmp)//'gauss.com1')
       call system ('mv '//trim(addresstmp)//'gauss.com1 '//trim(addresstmp)//'gauss.com')
     Endif
! <--- Udagawa End 2021.05.24

     call system(trim(address0)//'g0xrun_p '//trim(addresstmp)//'gauss.com '//trim(addresstmp)//'gauss.log '//trim(addresstmp))
     open(igauss+id,file=trim(addresstmp)//'gauss.log')

!  +++ Reading "SCF Done" +++
     do
       read(igauss+id,'(a)',end=401) line
       iline=index(line,trim(key1))
       if(iline > 0) exit
     end do
     index1 = index(line,'=')
     read(line(index1+2:index1+17),*) enetemp ! For DFT

     do
       read(igauss+id,'(a)',end=101) line
       iline=index(line,trim(key6))  ! Reading "SCE Done"
       if(iline > 0) exit
     end do
     read(line(38:60),*) enetemp ! For MP2
101 continue


     !if (theory == 0) then
     !  chl = index(line, '=')
     !  chr = index(line(chl+1:len(line)-chl), '.') + chl
     !  chr = index(line(chr+1:len(line)-chr), ' ') + chr
     !  read(line(chl+1:chr-1),*) enetemp  ! For DFT
     !else if (theory == 1) then
     !  read(line(38:60),*) enetemp  ! For MP2
     !end if

     Eenergy(imode2)=enetemp
     rewind(igauss+id)
!  +++ End Reading "SCF Done" +++

!  +++ Reading "Mulliken charge" +++
     if( Lsave_charge .eqv. .True.) then
     !if(nocharge==0) then
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
     endif
!  +++ End Reading "Mulliken charge" +++

!  +++ Reading "Dipole moment" +++
     if( Lsave_dipole .eqv. .True. ) then
     !if(nodipole==0) then
       do
         Read(igauss+id,'(a)',end=403) line
         iline=index(trim(line),trim(key3))
         if(iline > 0) exit
       end do
       Read(igauss+id,'(a)') line
       Read(line(20:26),*) dipoler(1,imode2)
       Read(line(46:52),*) dipoler(2,imode2)
       Read(line(72:78),*) dipoler(3,imode2)
       !Read(line(98:104),*) dipole(imode2)
     endif
!  +++ End Reading "Dipole moment" +++

!  +++ Reading "Isotropic Fermi" +++
     if( Lsave_hfcc .eqv. .True. ) then
     !if(nohfcc==0) then
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
!  +++ End Reading "Isotropic Fermi" +++

!  +++ Reading "Atomic Force" +++
!     Rewind(igauss+id)
     do
       Read(igauss+id,'(a)',end=402) line
       iline=index(trim(line),trim(key2))
       if(iline > 0) exit
     end do
     read(igauss+id,'()')

     do iatom2=1,natom
        read(igauss+id,'(a)') line
        !read(line(24:38),*) fx(iatom2,imode2)
        !read(line(39:53),*) fy(iatom2,imode2)
        !read(line(54:68),*) fz(iatom2,imode2)
        read(line(24:38),*) fr(1,iatom2,imode2)
        read(line(39:53),*) fr(2,iatom2,imode2)
        read(line(54:68),*) fr(3,iatom2,imode2)
     enddo
!  +++ End Reading "Atomic Force" +++

     close(igauss+id)

     call system('rm -rf '//trim(addresstmp)//'Gau*')
    enddo

! call print_result_cl

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


