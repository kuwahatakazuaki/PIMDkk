subroutine Force_VASP_MPI
  use Parameters
  use utility, only: program_abort, read_val_next, search_line
  implicit none

  character(Len=130) :: line
  character(Len=32)  :: key1, key2, key3
  character :: Cdummy(3)
  integer :: i,j,k
  integer :: iline, Imode, Iatom !, igauss = 20
  integer :: Uout, Uinp
  real(8) :: enetemp, dummy
  !real(8), parameter :: ev_to_hartree  = 1.0 / 27.21138505
  !real(8), parameter :: eVAng_HartBohr = 0.5291772108 / 27.21138505
  key1  = ('energy  without entropy')
  key2  = ('TOTAL-FORCE')
  key3  = ('external pressure')

  !Call Start_Recv_Send_MPI_tk

  do Imode = Ista, Iend

    write(addresstmp(laddress+1:laddress+6),'(i5.5,A1)') Imode,'/' 

!need file -> INCAR, POTCAR, KPOINTS, LATTICE-------------------
! +++ Calculation VASP +++
    call system('cp LATTICE '//trim(addresstmp)//'POSCAR')
    open(newunit=Uout,file=trim(addresstmp)//'POSCAR',status='old',position='append')
      do Iatom=1,natom
        write(Uout,*) r(:,Iatom,Imode)*AUtoAng
        !write(Uout,9999) r(:,Iatom,Imode)*AUtoAng
      enddo
      write(Uout,*)
    close(Uout)
    call system('./vasp_run.sh '//trim(addresstmp))
! +++ End Calculation VASP +++

! +++ Reading VASP output +++
    open(newunit=Uinp,file=trim(addresstmp)//'OUTCAR')

! +++ Reading "external pressure" +++
      do
        read(Uinp,'(a)',end=403) line
        iline=index(line,trim(key3))
        if(iline > 0) exit
      enddo
      read(line,*) Cdummy(1:3), pressure(Imode)
! +++ End  "external pressure" +++

! +++ Reading "TOTAL-FORCE" +++
!      rewind(Uinp)
      do
        Read(Uinp,'(a)',end=402) line
        iline=index(line(49:60),trim(key2))
        if(iline > 0) exit
      end do
      read(Uinp,'()')
      do Iatom = 1, natom
        read(Uinp,*) dummy, dummy, dummy, fr(:,Iatom,Imode)
      enddo
! +++ End Reading "TOTAL-FORCE" +++

! +++ Reading "energy  without entropy" +++
      do
        read(Uinp,'(a)',end=401) line  !read OUTCAR
        iline=index(line,trim(key1))
        if(iline > 0) exit
      enddo
      read(line(32:45),*) enetemp
      !enetemp = enetemp * ev_to_hartree
      Eenergy(Imode) = enetemp * eVtoAU
! +++ End Reading "energy  without entropy" +++

    close(Uinp)

    !fr(:,:,Imode)=fr(:,:,Imode)*eVAng_HartBohr*dp_inv
    fr(:,:,Imode)=fr(:,:,Imode)*eVAng2AU*dp_inv
  enddo


!9999 format(3F24.16)
!9998 format(3E23.15)
!9997 format(2E23.15)
!9996 format(E23.15)
!9995 format(4E23.15)

return
401 print *, 'ERROR!!: We can not find "energy  without entropy" in VASP output'; stop
402 print *, 'ERROR!!: We can not find "TOTAL-FORCE" in VASP output'; stop
403 print *, 'ERROR!!: We can not find "external pressure" in VASP output'; stop
end subroutine



