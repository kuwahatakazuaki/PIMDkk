subroutine Force_Gaussian_MPI_tk
  use Parameters
  use utility, only: program_abort
  implicit none

  Character (Len=90) :: key1, key2, key3, key4, key5, key6
  character(len=120) :: line
  Integer            :: iline, id,imode2,iatom2
  Integer            :: index1
  Double Precision   :: enetemp
  Integer :: i,j,k
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
  end if

  id=0
!  Call Start_Recv_Send_MPI_tk  ! Gather r (Coordinate)

  do imode2=ista,iend
    write(addresstmp(laddress+1:laddress+6),'(i5.5,A1)') imode2,'/'

    call system('cat '//trim(addresstmp)//'gauss.tmp1 > '//trim(addresstmp)//'gauss.com')


!    open(igauss+id,file=trim(addresstmp)//'gauss.xyz',status='unknown')
    open(igauss+id,file=trim(addresstmp)//'gauss.com',status='old',position='append')
      do iatom2=1,natom
        write(igauss+id,*) alabel(iatom2),r(:,iatom2,imode2)*AUtoAng
      enddo
      write(igauss+id,*)
    close(igauss+id)


    !If(NGenGau==1) Then
    !   call system('cat '//trim(addresstmp)//'gauss.bss >> '//trim(addresstmp)//'gauss.com')
    !EndIf

! Udagawa Start 2021.05.24 --->
    !If(istepsv == 0 .OR. (nRestart ==1 .AND. istepsv == nrstep+1)) then
    If(istepsv == 0 .OR. ( (Lrestart .eqv. .True. ) .AND. istepsv == nrstep+1)) then
      call system ('sed -e "s/[Gg][Uu][Ee][Ss][Ss]=[Rr][Ee][Aa][Dd]//g" '//trim(addresstmp)//'gauss.com  &
        & > '//trim(addresstmp)//'gauss.com1')
      call system ('mv '//trim(addresstmp)//'gauss.com1 '//trim(addresstmp)//'gauss.com')
    Endif
! <--- Udagawa End 2021.05.24

! kuwahata 2021/06/06 for ITO
    call system(trim(address0)//'g0xrun_p '//trim(addresstmp)//'gauss.com '//trim(addresstmp)//'gauss.log '//trim(addresstmp))
    !call system(trim(addresstmp)//'g0xrun_p '//trim(addresstmp)//'gauss.com '//trim(addresstmp)//'gauss.log '//trim(addresstmp))
! End kuwahata 2021/06/06 for ITO

    open(igauss+id,file=trim(addresstmp)//'gauss.log')

!  +++ Reading "SCF Done" or "EUMP2" +++
!  +++ Automatically read EUMP2 with MP2 method +++
!  +++ We don't need to change 'theory' option  +++
     do
       read(igauss+id,'(a)',end=401) line
       iline=index(line,trim(key1))  ! Reading "SCE Done"
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

     Eenergy(imode2)=enetemp
     rewind(igauss+id)
!  +++ End Reading "SCF Done" +++

!  +++ Reading "Mulliken charge" +++
     if( Lsave_charge .eqv. .True.) then
       do
         read(igauss+id,'(a)',end=404) line
         iline=index(trim(line),trim(key4))  ! Reading "Mulliken charge"
         if(iline > 0) exit
       end do
       read(igauss+id,*)
       do iatom2=1,natom
          Read(igauss+id,'(a)') line
          Read(line(11:21),*) charge(iatom2,imode2)
       enddo
     endif
!  +++ End Reading "Mulliken charge" +++

!  +++ Reading "Dipole moment" +++
     if( Lsave_dipole .eqv. .True. ) then
       do
         read(igauss+id,'(a)',end=403) line
         iline=index(trim(line),trim(key3))  ! Reading "Dipole moment"
         if(iline > 0) exit
       end do
       Read(igauss+id,'(a)') line
       Read(line(20:26),*) dipoler(1,imode2)
       Read(line(46:52),*) dipoler(2,imode2)
       Read(line(72:78),*) dipoler(3,imode2)
     endif
!  +++ End Reading "Dipole moment" +++

!  +++ Reading "Isotropic Fermi" +++
     if( Lsave_hfcc .eqv. .True. ) then
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
!     rewind(igauss+id)
     do
       read(igauss+id,'(a)',end=402) line  ! Reading "Force"
       iline=index(line,trim(key2))
       if(iline > 0) exit
     end do
     read(igauss+id,*)
     do iatom2=1,natom
       read(igauss+id,'(a)') line
       read(line(24:38),*) fr(1,iatom2,imode2)
       read(line(39:53),*) fr(2,iatom2,imode2)
       read(line(54:68),*) fr(3,iatom2,imode2)
     enddo
!  +++ End Reading "Atomic Force" +++

     close(igauss+id)

     call system('rm -rf '//trim(addresstmp)//'Gau*')

     fr(:,:,imode2) = fr(:,:,imode2) * dp_inv
   enddo
  fr(:,:,Ista:Iend) = fr(:,:,Ista:Iend) * dp_inv

! Call Start_Send_Recv_MPI_tk

return
401 print *, 'ERROR!!: We can not find "SCF Done" or "EUMP2" in Gaussian output'; stop
402 print *, 'ERROR!!: We can not find "Force" in Gaussian output'; stop
403 print *, 'ERROR!!: We can not find "Dipole moment" in Gaussian output'; stop
404 print *, 'ERROR!!: We can not find "Mulliken charges" in Gaussian output'; stop
405 print *, 'ERROR!!: We can not find "Fermi Contact" in Gaussian output'; stop
end subroutine

