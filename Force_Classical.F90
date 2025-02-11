Subroutine Force_Classical
  use Parameters
  use utility, only: program_abort
  implicit none

  r(:,:,1) = ur(:,:,1)

  select case(Iforce)
    case(1)
      call Force_MOPAC
    case(6)
      call Force_Gaussian_classical
    case(8)
      call Force_VASP_cl
    case(15)
      call Force_Harmonic
    case(16)
      call Force_Double_Morse

    case(21)
      call force_nnp_araidai
    case(22)
      call force_nnp_matlantis
#ifdef _LAMMPS_
    case(24)
      call force_LAMMPS
#endif

    case(31)
      call force_spcf
    case default
      call program_abort('ERROR!!! Wrong "Iforce" option')
  end select

  if (Iumb > 0) call calc_umbrella

  fur(:,:,1)=fr(:,:,1)

  !potential = sum(Eenergy(:)) * dp_inv
  potential = Eenergy(1)

  Return
contains

! +++++++++++++++++++++++
! +++++ Force_MOPAC +++++
! +++++++++++++++++++++++
  subroutine Force_MOPAC
    Use Parameters
    use utility, only: program_abort
    implicit none
    integer :: imode

    Character    (Len=90)              :: key1, key2, key3, key4, key5, key6, key7, line
    Character    (Len=5)               :: nobeadtmp
    Character    (Len=3)               :: atemp
    character(len=10)  :: Cdummy(10)
    integer :: Idummy
    real(8) :: Ddummy
    integer :: i, imopac = 20
    Integer                            :: itemp,iatom2,icount,id,iline
    Double Precision                   :: enetemp
    key1  = ('TOTAL ENERGY')
! Need to change the "Key2" depend on version of Mopac
!    key2  = ('PARAMETER     ATOM    TYPE            VALUE       GRADIENT')
    key2  = ('PARAMETER     ATOM    TYPE            VALUE        GRADIENT') !For  MOPAC2016 version
    key3  = ('DIPOLE           X         Y         Z       TOTAL')
    key4  = ('MULLIKEN POPULATIONS AND CHARGES')
    key5  = ('HOMO LUMO ENE')

    id=0
    do imode=ista,iend
      write(addresstmp(laddress+1:laddress+6),'(i5.5,A1)') imode,'/'
       call system('cp mopac.mop '//trim(addresstmp))


      open(imopac+id,file=trim(addresstmp)//'mopac.mop',status='old',position='append')
        do i = 1, natom
         write(imopac+ID,*) alabel(i),r(:,i,imode)*AUtoAng
        end do
      close(imopac)

     call system( './g0xrun_p '//trim(addresstmp))
     open(imopac+ID,file=trim(addresstmp)//'mopac.out')
!  +++ Reading "ENERGY" +++
     do
       read(imopac+id,'(a)',end=401) line
       if ( index(line,trim(key1)) > 0 ) exit  ! Reading "SCE Done"
     end do
     read(line,*) Cdummy(1:3), enetemp
     Eenergy(imode) = enetemp * 0.0367493238d0
!  +++ End Reading "ENERGY" +++

!  +++ Reading "Atomic Force" +++
    do
       read(imopac+id,'(a)',end=402) line
       if ( index(line,trim(key2)) > 0 ) exit
    end do
    do i = 1, natom
      read(imopac+id,*) Idummy, Idummy, Cdummy(1:3), Ddummy, fr(1,i,imode)
      read(imopac+id,*) Idummy, Idummy, Cdummy(1:3), Ddummy, fr(2,i,imode)
      read(imopac+id,*) Idummy, Idummy, Cdummy(1:3), Ddummy, fr(3,i,imode)
    end do
    fr(:,:,imode) = -fr(:,:,imode)*0.00159362D0 * AUtoAng * dp_inv
!  +++ End Reading "Atomic Force" +++


!  +++ Reading "Dipole moment" +++
    if( Lsave_dipole .eqv. .True. ) then
      do
         read(imopac+id,'(a)',end=403) line
         if ( index(line,trim(key3)) > 0 ) exit
      end do

! Udagawa Start 2021.05.24 --->
      do i = 1, 2
        read(imopac+id,*) line
      end do
      read(imopac+id,*) line, dipoler(:,imode)
! <--- Udagawa End 2021.05.24
    end if
!  +++ End Reading "Dipole moment" +++
     close(imopac+ID)
       call system('rm -rf '//trim(addresstmp)//'mopac.*')
    enddo

    return
401 call program_abort('ERROR!!: We can not find "TOTAL ENERGY" in mopac.out')
402 call program_abort('ERROR!!: We can not find "Atomic Force" in mopac.out')
403 call program_abort('ERROR!!: We can not find "DIPOLE" in mopac.out')
  end subroutine Force_MOPAC
! +++++++++++++++++++++++++++
! +++++ End Force_MOPAC +++++
! +++++++++++++++++++++++++++

! +++++++++++++++++++++++++
! +++++ Force_VASP_cl +++++
! +++++++++++++++++++++++++
  subroutine Force_VASP_cl
  use Parameters
  implicit none

  character(Len=130) :: line
  character(Len=32)  :: key1, key2, key3
  character :: Cdummy(3)
  Integer :: i,j,k
  Integer :: iline, id,imode2,iatom2
  integer :: igauss = 20
  real(8) :: enetemp, dummy
  real(8), parameter :: ev_to_hartree  = 1.0 / 27.21138505
  real(8), parameter :: eVAng_HartBohr = 0.5291772108 / 27.21138505
  key1  = ('energy  without entropy')
  key2  = ('TOTAL-FORCE')
  key3  = ('external pressure')

  id=0
!  Call Start_Recv_Send_MPI_tk

  do imode2=ista,iend

    write(addresstmp(laddress+1:laddress+6),'(i5.5,A1)') imode2,'/'

!need file -> INCAR, POTCAR, KPOINTS, LATTICE-------------------
! +++ Calculation VASP +++
    call system('cp LATTICE '//trim(addresstmp)//'POSCAR')

    open(igauss+id,file=trim(addresstmp)//'POSCAR',status='old',position='append')
      do iatom2=1,natom
        write(igauss+id,9999) r(:,iatom2,imode2)*AUtoAng
      enddo
      write(igauss+id,*)
    close(igauss+id)

    call system('./vasp_run.sh '//trim(addresstmp))
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
      read(igauss+id,'()') !skip 1line
      do iatom2=1,natom
        read(igauss+id,*) dummy, dummy, dummy, fr(:,iatom2,imode2)
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
    fr(:,:,imode2)=fr(:,:,imode2)*eVAng_HartBohr*dp_inv

  enddo


  9999 format(3F24.16)
  9998 format(3E23.15)
  9997 format(2E23.15)
  9996 format(E23.15)
  9995 format(4E23.15)

  return
  401 print *, 'ERROR!!: We can not find "energy  without entropy" in VASP output'; stop
  402 print *, 'ERROR!!: We can not find "TOTAL-FORCE" in VASP output'; stop
  403 print *, 'ERROR!!: We can not find "external pressure" in VASP output'; stop
  end subroutine Force_VASP_cl

! +++++++++++++++++++++++++
! +++++ Force_VASP_cl +++++
! +++++++++++++++++++++++++
end subroutine Force_Classical
