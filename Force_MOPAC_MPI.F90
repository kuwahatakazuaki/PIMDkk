Subroutine Force_MOPAC_MPI
    Use Parameters
    use utility, only: program_abort
    Implicit None

    Character    (Len=90)              :: key1, key2, key3, key4, key5, key6, key7, line
    Character    (Len=3)               :: atemp
    character(len=10)  :: Cdummy(10)
    integer :: Idummy
    real(8) :: Ddummy
    Integer                            :: itemp,iatom2,icount,id,iline
    integer :: i, j, k, imode
    integer :: imopac = 20
    Double Precision                   :: enetemp
    key1  = ('TOTAL ENERGY')
!    key2  = ('PARAMETER     ATOM    TYPE            VALUE       GRADIENT')
    key2  = ('PARAMETER     ATOM    TYPE            VALUE        GRADIENT') !For MOPAC2016 version
    key3  = ('DIPOLE           X         Y         Z       TOTAL')
    key4  = ('MULLIKEN POPULATIONS AND CHARGES')
    key5  = ('HOMO LUMO ENE')


    id=0
  !Call Start_Recv_Send_MPI_tk
    do imode=ista,iend
      write(addresstmp(laddress+1:laddress+6),'(i5.5,A1)') imode,'/'
       call system('cp mopac.mop '//trim(addresstmp))

      open(imopac+id,file=trim(addresstmp)//'mopac.mop',status='old',position='append')
        do i = 1, natom
         write(imopac+ID,*) alabel(i),r(:,i,imode)*AUtoAng
        end do
      close(imopac)

!       call system('/opt/mopac/MOPAC2009.exe '//trim(address)//'/'//nobeadtmp//'/mopac.mop')
     call system( './g0xrun_p '//trim(addresstmp))
     open(imopac+ID,file=trim(addresstmp)//'mopac.out')

!  +++ Reading "ENERGY" +++
     do
       read(imopac+id,'(a)',end=401) line
       if ( index(line,trim(key1)) > 0 ) exit  ! Reading "SCE Done"
     end do
     read(line,*) Cdummy(1:3), enetemp
     pot_bead(imode) = enetemp * 0.0367493238d0
!  +++ End Reading "ENERGY" +++

!  +++ Reading "Atomic Force" +++
    do
       read(imopac+id,'(a)',end=402) line
       if ( index(line,trim(key2)) > 0 ) exit
    end do
    do i = 1, natom
      !read(imopac+id,*) Idummy, Idummy, Cdummy(1:3), Ddummy, fx(i,imode)
      !read(imopac+id,*) Idummy, Idummy, Cdummy(1:3), Ddummy, fy(i,imode)
      !read(imopac+id,*) Idummy, Idummy, Cdummy(1:3), Ddummy, fz(i,imode)
      read(imopac+id,*) Idummy, Idummy, Cdummy(1:3), Ddummy, fr(1,i,imode)
      read(imopac+id,*) Idummy, Idummy, Cdummy(1:3), Ddummy, fr(2,i,imode)
      read(imopac+id,*) Idummy, Idummy, Cdummy(1:3), Ddummy, fr(3,i,imode)
    end do
    fr(:,:,imode) = -fr(:,:,imode)*0.00159362D0 * AUtoAng * dp_inv
!  +++ End Reading "Atomic Force" +++


!  +++ Reading "Dipole moment" +++
     if( Lsave_dipole .eqv. .True. ) then
    !if ( nodipole == 0 ) then
      do
         read(imopac+id,'(a)',end=403) line
         if ( index(line,trim(key3)) > 0 ) exit
      end do
! Udagawa Start 2021.05.24 --->
      do i = 1, 2
        read(imopac+id,*) line
      end do
      read(imopac+id,*) line, dipoler(:,imode)
      !read(imopac+id,*) line, dipolex(imode), dipoley(imode), dipolez(imode),dipole(imode)
! <--- Udagawa End 2021.05.24
!      read(imopac+id,*) dipolex(imode), dipoley(imode), dipolez(imode),dipole(imode)
    end if
!  +++ End Reading "Dipole moment" +++


       close(imopac+ID)
       call system('rm -rf '//trim(addresstmp)//'mopac.*')
    enddo

!Call Start_Send_Recv_MPI_tk


!    potential=0.0D+00
!    do imode=1,nbead
!       potential=potential+pot_bead(imode)
!    enddo
!    potential=potential/dp

return
401 call program_abort('ERROR!!: We can not find "TOTAL ENERGY" in mopac.out')
402 call program_abort('ERROR!!: We can not find "Atomic Force" in mopac.out')
403 call program_abort('ERROR!!: We can not find "DIPOLE" in mopac.out')
end subroutine Force_MOPAC_MPI


