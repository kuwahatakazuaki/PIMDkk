subroutine Force_Gaussian_MPI_tk
  use Parameters
  use utility, only: program_abort
  implicit none

  character(len=:), allocatable :: key1, key2, key3, key4, key5, key6
  character(len=120) :: line
  integer            :: iline, imode, iatom
  Double Precision   :: enetemp
  integer :: i,j,k, igauss
  !integer :: igauss = 20
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


!    open(igauss,file=trim(addresstmp)//'gauss.xyz',status='unknown')
    open(newunit=igauss,file=trim(addresstmp)//'gauss.com',status='old',position='append')
      do iatom = 1, Natom
        write(igauss,*) alabel(iatom),r(:,iatom,imode)*AUtoAng
      end do
      write(igauss,*)
    close(igauss)


    !If(NGenGau==1) Then
    !   call system('cat '//trim(addresstmp)//'gauss.bss >> '//trim(addresstmp)//'gauss.com')
    !EndIf

! Udagawa Start 2021.05.24 --->
    !If(istepsv == 0 .OR. (nRestart ==1 .AND. istepsv == nrstep+1)) then
    if(istepsv == 0 .OR. ( (Lrestart .eqv. .True. ) .AND. istepsv == nrstep+1)) then
      call system ('sed -e "s/[Gg][Uu][Ee][Ss][Ss]=[Rr][Ee][Aa][Dd]//g" '//trim(addresstmp)//'gauss.com  &
        & > '//trim(addresstmp)//'gauss.com1')
      call system ('mv '//trim(addresstmp)//'gauss.com1 '//trim(addresstmp)//'gauss.com')
    end if
! <--- Udagawa End 2021.05.24

! kuwahata 2021/06/06 for ITO
    call system(trim(address0)//'g0xrun_p '//trim(addresstmp)//'gauss.com '//trim(addresstmp)//'gauss.log '//trim(addresstmp))
    !call system(trim(addresstmp)//'g0xrun_p '//trim(addresstmp)//'gauss.com '//trim(addresstmp)//'gauss.log '//trim(addresstmp))
! End kuwahata 2021/06/06 for ITO

    open(newunit=igauss,file=trim(addresstmp)//'gauss.log')

!  +++ Reading "SCF Done" or "EUMP2" +++
!  +++ Automatically read EUMP2 with MP2 method +++
    call search_line(igauss,key1,line)
    call read_val_next(line,'=',enetemp)

    do
      read(igauss,'(a)',end=101) line
      iline=index(line,trim(key6))  ! Reading "SCE Done"
      if(iline > 0) exit
    end do
    read(line(38:60),*) enetemp ! For MP2
101 continue

    Eenergy(imode) = enetemp
    rewind(igauss)
!  +++ End Reading "SCF Done" +++

!print *, imode, enetemp
!  +++ Reading "Mulliken charge" +++
    if ( Lsave_charge .eqv. .True.) then
      call search_line(igauss,key4,line)
      read(igauss,'()')
      do i = 1, Natom
        read(igauss,*) dummyC(1:2), charge(i,imode)
      end do
     ! read(igauss,*)
     ! do iatom=1,natom
     !    Read(igauss,'(a)') line
     !    Read(line(11:21),*) charge(iatom,imode)
     ! enddo
    end if
!  +++ End Reading "Mulliken charge" +++

!print *, imode, charge(1,1)
!  +++ Reading "Dipole moment" +++
    if ( Lsave_dipole .eqv. .True. ) then
      call search_line(igauss,key3,line)
      read(igauss,*) dummyC(1), dipoler(1,imode), &
                     dummyC(2), dipoler(2,imode), &
                     dummyC(3), dipoler(3,imode)
    end if
!  +++ End Reading "Dipole moment" +++

!  +++ Reading "Isotropic Fermi" +++
    if( Lsave_hfcc .eqv. .True. ) then
      call search_line(igauss,key5,line)
      !do
      !  Read(igauss,'(a)',end=405) line
      !  iline=index(trim(line),trim(key5))
      !  if(iline > 0) exit
      !end do
      read(igauss,'()')
      do iatom=1,natom
        Read(igauss,'(a)') line
        Read(line(36:49),*) hfcc(iatom,imode)
      enddo
    endif
!  +++ End Reading "Isotropic Fermi" +++

!  +++ Reading "Atomic Force" +++
    call search_line(igauss,key2,line)
    read(igauss,'()')
    do i = 1, Natom
      read(igauss,*) dummyC(1:2), fr(:,i,imode)
    end do
    !read(igauss,*)
    !do iatom=1,natom
    !  read(igauss,'(a)') line
    !  read(line(24:38),*) fr(1,iatom,imode)
    !  read(line(39:53),*) fr(2,iatom,imode)
    !  read(line(54:68),*) fr(3,iatom,imode)
    !enddo
!  +++ End Reading "Atomic Force" +++

    close(igauss)
    call system('rm -rf '//trim(addresstmp)//'Gau*')

    !fr(:,:,imode) = fr(:,:,imode) * dp_inv
  enddo
  fr(:,:,Ista:Iend) = fr(:,:,Ista:Iend) * dp_inv


return
401 print *, 'ERROR!!: We can not find "SCF Done" or "EUMP2" in Gaussian output'; stop
402 print *, 'ERROR!!: We can not find "Force" in Gaussian output'; stop
403 print *, 'ERROR!!: We can not find "Dipole moment" in Gaussian output'; stop
404 print *, 'ERROR!!: We can not find "Mulliken charges" in Gaussian output'; stop
405 print *, 'ERROR!!: We can not find "Fermi Contact" in Gaussian output'; stop

contains

  subroutine read_val_next(line,key,val)
    character(len=*), intent(in) :: key, line
    real(8), intent(out) :: val
    integer :: pos
    pos = index(line,key)
    read( line(pos+1:),* ) val
  end subroutine read_val_next

  subroutine search_line(Iunit,key,line)
    integer, intent(in) :: Iunit
    character(len=*), intent(in) :: key
    character(len=*), intent(out) :: line

    do
      read(Iunit,'(a)',end=401) line
      if (index(line,key) > 0 ) exit
    end do
    return
    401 print *, 'ERROR!!: There is no key: "'//key//'"' ; stop
  end subroutine search_line

end subroutine

