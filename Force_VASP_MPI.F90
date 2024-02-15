subroutine Force_VASP_MPI
  use Parameters
  !use MPI
  implicit none

  character(Len=130) :: line
  character(Len=32)  :: key1, key2, key3
  character :: Cdummy(3)
  Integer :: i,j,k
  Integer :: iline, id,imode2,iatom2, igauss = 20
  real(8) :: enetemp, dummy
  real(8), parameter :: ev_to_hartree  = 1.0 / 27.21138505
  real(8), parameter :: eVAng_HartBohr = 0.5291772108 / 27.21138505
  key1  = ('energy  without entropy')
  key2  = ('TOTAL-FORCE')
  key3  = ('external pressure')

  id=0
  Call Start_Recv_Send_MPI_tk

  do imode2=ista,iend

    write(addresstmp(laddress+1:laddress+6),'(i5.5,A1)') imode2,'/' 

!need file -> INCAR, POTCAR, KPOINTS, LATTICE-------------------
! +++ Calculation VASP +++

    call system('cp LATTICE '//trim(addresstmp)//'POSCAR')
! For Wisteria
!call system('cp  vasp_std '//trim(addresstmp)//'vasp_std')
!call system('cp  vasp.sh '//trim(addresstmp)//'vasp.sh')
! End Wisteria

    open(igauss+id,file=trim(addresstmp)//'POSCAR',status='old',position='append')
      do iatom2=1,natom
        write(igauss+id,9999) r(:,iatom2,imode2)*bohr_inv
        !write(igauss+id,9999) x(iatom2,imode2)*bohr_inv, y(iatom2,imode2)*bohr_inv, z(iatom2,imode2)*bohr_inv
      enddo
      write(igauss+id,*)
    close(igauss+id)

    call system('./vasp_run.sh '//trim(addresstmp))
! For Wisteria
!call system( 'cd  '//trim(addresstmp)//'; pjsub vasp.sh' )
!call system( 'cd  '//trim(addresstmp)//'; pwd;  mpiexec  -n 1 ./vasp_std' )
! End Wisteria

! +++ End Calculation VASP +++

! +++ Reading VASP output +++
    open(igauss+id,file=trim(addresstmp)//'OUTCAR')

! +++ Reading "external pressure" +++
      do
        read(igauss+id,'(a)',end=403) line
        iline=index(line,trim(key3))
        if(iline > 0) exit
      enddo
      read(line,*) Cdummy(1:3), pressure(imode2)
! +++ End  "external pressure" +++

! +++ Reading "TOTAL-FORCE" +++
!      rewind(igauss+id)
      do
        Read(igauss+id,'(a)',end=402) line
        iline=index(line(49:60),trim(key2))
        if(iline > 0) exit
      end do
      read(igauss+id,'()')
      do iatom2=1,natom
        read(igauss+id,*) dummy, dummy, dummy, fr(:,iatom2,imode2)
        !read(igauss+id,*) dummy, dummy, dummy, fx(iatom2,imode2), fy(iatom2,imode2), fz(iatom2,imode2)
      enddo
! +++ End Reading "TOTAL-FORCE" +++

! +++ Reading "energy  without entropy" +++
      do
        read(igauss+id,'(a)',end=401) line  !read OUTCAR
        iline=index(line,trim(key1))
        if(iline > 0) exit
      enddo
      read(line(32:45),*) enetemp
      enetemp = enetemp * ev_to_hartree
      Eenergy(imode2)=enetemp
! +++ End Reading "energy  without entropy" +++

    close(igauss+id)

    !fx(:,imode2)=fx(:,imode2)*eVAng_HartBohr*dp_inv
    !fy(:,imode2)=fy(:,imode2)*eVAng_HartBohr*dp_inv
    !fz(:,imode2)=fz(:,imode2)*eVAng_HartBohr*dp_inv
    fr(:,:,imode2)=fr(:,:,imode2)*eVAng_HartBohr*dp_inv
!print *, imode2, pressure(imode2)
  enddo

  call Start_Send_Recv_MPI_tk

9999 format(3F24.16)
9998 format(3E23.15)
9997 format(2E23.15)
9996 format(E23.15)
9995 format(4E23.15)

return
401 print *, 'ERROR!!: We can not find "energy  without entropy" in VASP output'; stop
402 print *, 'ERROR!!: We can not find "TOTAL-FORCE" in VASP output'; stop
403 print *, 'ERROR!!: We can not find "external pressure" in VASP output'; stop
end subroutine



