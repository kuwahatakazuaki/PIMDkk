subroutine Force_Gaussian_MPI_tk
  use Parameters
  use utility, only: program_abort, read_val_next, search_line
  implicit none

  character(len=:), allocatable :: key1, key2, key3, key4, key5, key6
  character(len=120) :: line
  integer            :: iline, imode, iatom
  Double Precision   :: enetemp
  integer :: i,j,k, Uinp
  character :: dummyC(10)

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

!  Call Start_Recv_Send_MPI_tk  ! Gather r (Coordinate)

  do imode=ista,iend
    write(addresstmp(laddress+1:laddress+6),'(i5.5,A1)') imode,'/'

    call system('cat '//trim(addresstmp)//'gauss.tmp1 > '//trim(addresstmp)//'gauss.com')

!    open(Uinp,file=trim(addresstmp)//'gauss.xyz',status='unknown')
    open(newunit=Uinp,file=trim(addresstmp)//'gauss.com',status='old',position='append')
      do iatom = 1, Natom
        write(Uinp,*) alabel(iatom),r(:,iatom,imode)*AU2Ang
      end do
      write(Uinp,*)
    close(Uinp)


    !If(NGenGau==1) Then
    !   call system('cat '//trim(addresstmp)//'gauss.bss >> '//trim(addresstmp)//'gauss.com')
    !EndIf

! Udagawa Start 2021.05.24 --->
    if(istepsv == 0 .OR. ( (Lrestart .eqv. .True. ) .AND. istepsv == Irestep+1)) then
      call system ('sed -e "s/[Gg][Uu][Ee][Ss][Ss]=[Rr][Ee][Aa][Dd]//g" '//trim(addresstmp)//'gauss.com  &
        & > '//trim(addresstmp)//'gauss.com1')
      call system ('mv '//trim(addresstmp)//'gauss.com1 '//trim(addresstmp)//'gauss.com')
    end if
! <--- Udagawa End 2021.05.24

! kuwahata 2021/06/06 for ITO
    call system(trim(address0)//'g0xrun_p '//trim(addresstmp)//'gauss.com '//trim(addresstmp)//'gauss.log '//trim(addresstmp))
    !call system(trim(addresstmp)//'g0xrun_p '//trim(addresstmp)//'gauss.com '//trim(addresstmp)//'gauss.log '//trim(addresstmp))
! End kuwahata 2021/06/06 for ITO

    open(newunit=Uinp,file=trim(addresstmp)//'gauss.log')

!  +++ Reading "SCF Done" or "EUMP2" +++
!  +++ Automatically read EUMP2 with MP2 method +++
    call search_line(Uinp,key1,line)
    call read_val_next(line,'=',enetemp)

    do
      read(Uinp,'(a)',end=101) line
      iline=index(line,trim(key6))  ! Reading "SCE Done"
      if(iline > 0) exit
    end do
    read(line(38:60),*) enetemp ! For MP2
101 continue

    pot_bead(imode) = enetemp
    rewind(Uinp)
!  +++ End Reading "SCF Done" +++

!  +++ Reading "Mulliken charge" +++
    if ( Lsave_charge .eqv. .True.) then
      call search_line(Uinp,key4,line)
      read(Uinp,'()')
      do i = 1, Natom
        read(Uinp,*) dummyC(1:2), charge(i,imode)
      end do
    end if
!  +++ End Reading "Mulliken charge" +++

!  +++ Reading "Dipole moment" +++
    if ( Lsave_dipole .eqv. .True. ) then
      call search_line(Uinp,key3,line)
      read(Uinp,*) dummyC(1), dipoler(1,imode), &
                     dummyC(2), dipoler(2,imode), &
                     dummyC(3), dipoler(3,imode)
    end if
!  +++ End Reading "Dipole moment" +++

!  +++ Reading "Isotropic Fermi" +++
    if( Lsave_hfcc .eqv. .True. ) then
      call search_line(Uinp,key5,line)
      read(Uinp,'()')
      do iatom=1,natom
        Read(Uinp,'(a)') line
        Read(line(36:49),*) hfcc(iatom,imode)
      enddo
    endif
!  +++ End Reading "Isotropic Fermi" +++

!  +++ Reading "Atomic Force" +++
    call search_line(Uinp,key2,line)
    read(Uinp,'()')
    do i = 1, Natom
      read(Uinp,*) dummyC(1:2), fr(:,i,imode)
    end do
!  +++ End Reading "Atomic Force" +++

    close(Uinp)
    call system('rm -rf '//trim(addresstmp)//'Gau*')

  enddo
  fr(:,:,Ista:Iend) = fr(:,:,Ista:Iend) * dp_inv


return
401 print *, 'ERROR!!: We can not find "SCF Done" or "EUMP2" in Gaussian output'; stop
402 print *, 'ERROR!!: We can not find "Force" in Gaussian output'; stop
403 print *, 'ERROR!!: We can not find "Dipole moment" in Gaussian output'; stop
404 print *, 'ERROR!!: We can not find "Mulliken charges" in Gaussian output'; stop
405 print *, 'ERROR!!: We can not find "Fermi Contact" in Gaussian output'; stop

!contains

end subroutine

