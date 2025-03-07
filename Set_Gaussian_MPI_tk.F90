subroutine Set_Gaussian_MPI_tk
  use Parameters
  use utility, only: program_abort
  implicit none
  integer :: i,j,k,imode, Uout
  integer :: access
  logical :: Ltemp, Lg0xrun

  if ( MyRank == 0 ) then
    inquire(FILE='./gauss.tmp',EXIST=Ltemp)
    inquire(FILE='./g0xrun_p', EXIST=Lg0xrun)
    if ( Ltemp   .eqv. .false. ) call err_set_tmp
    if ( Lg0xrun .eqv. .false. ) call err_set_g0xrun
    if ( (Ltemp .and. Lg0xrun) .eqv. .false. ) &
        call program_abort('== Stop Input files are missing ==')
    !if ( access("./gauss.tmp", " ") /= 0) call err_set_tmp
    !if ( access("./g0xrun_p", " ")  /= 0) call err_set_g0xrun
  end if

  do imode=ista,iend
     write(addresstmp(laddress+1:laddress+6),'(i5.5,A1)') imode,'/'
     call system('mkdir -p '//trim(addresstmp))
     call system('cp gauss.tmp '//trim(addresstmp))
! kuwahata 2021/06/06 for ITO
     call system('cp -f g0xrun_p '//trim(address0))
! End kuwahata 2021/06/06 for ITO
     !If(NGenGau==1) Then
     !   call system('cp gauss.bss '//trim(addresstmp))
     !EndIf
     open(newunit=Uout,file=trim(addresstmp)//'gauss.tmp1',status='unknown')
       write(Uout,'(a)') '%Chk='//trim(addresstmp)//'gauss.chk'
       write(Uout,'(a)') '%RWF='//trim(addresstmp)
       write(Uout,'(a)') '%Int='//trim(addresstmp)
       write(Uout,'(a)') '%D2E='//trim(addresstmp)
     close(Uout)
     call system('cat '//trim(addresstmp)//'gauss.tmp >> '//trim(addresstmp)//'gauss.tmp1')

  !! Udagawa Start 2021.05.24 --->
  !     If(istepsv == 0 .OR. (nRestart ==1 .AND. istepsv == nrstep+1)) then
  !       call system ('sed -e "s/[Gg][Uu][Ee][Ss][Ss]=[Rr][Ee][Aa][Dd]//g" '//trim(addresstmp)//'gauss.tmp1  &
  !         & > '//trim(addresstmp)//'gauss.tmp2')
  !       call system ('mv '//trim(addresstmp)//'gauss.tmp2 '//trim(addresstmp)//'gauss.tmp1')
  !     Endif
  !! <--- Udagawa End 2021.05.24
  enddo

  return
contains

  subroutine err_set_tmp
    integer :: Nunit
    open(newunit=Nunit,file='gauss.tmp',status='new')
      write(Nunit,'(a)') "# B3LYP/6-31G force SCF=XQC"
      write(Nunit,'(a)') ""
      write(Nunit,'(a)') "title"
      write(Nunit,'(a)') ""
      write(Nunit,'(a)') "0 1"
    close(Nunit)
    print *, 'ERROR!! "gauss.tmp" is NOT exist!'
    print *, 'Please use "gauss.tmp" template'
  end subroutine err_set_tmp

  subroutine err_set_g0xrun
    integer :: Nunit
    open(newunit=Nunit,file='g0xrun_p',status='new')
      write(Nunit,'(a)') "#!/bin/bash"
      write(Nunit,'(a)') "export GAUSS_SCRDIR=$3"
      write(Nunit,'(a)') "/usr/local/g16b01/pgi/g16/g16 < $1 >& $2"
    close(Nunit)
    print *, 'ERROR!! "g0xrun_p" is NOT exist!'
    print *, 'Please use "g0xrun_p" template'
  end subroutine err_set_g0xrun

end subroutine Set_Gaussian_MPI_tk

