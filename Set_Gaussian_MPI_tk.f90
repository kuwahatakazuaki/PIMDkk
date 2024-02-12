Subroutine Set_Gaussian_MPI_tk
  Use MPI
  Use Parameters
  use utility, only: program_abort
  Implicit None
  Integer   :: i,j,k,id, imode
  integer :: access
  integer :: igauss = 20

if ( MyRank == 0 ) then
  if      ( access("./gauss.tmp", " ") .ne. 0) then
    call err_set_tmp
  else if ( access("./g0xrun_p", " ") .ne. 0) then
    call err_set_g0xrun
  end if
end if


  id=0
do imode=ista,iend
   write(addresstmp(laddress+1:laddress+6),'(i5.5,A1)') imode,'/'
   call system('mkdir -p '//trim(addresstmp))
   call system('cp gauss.tmp '//trim(addresstmp))
! kuwahata 2021/06/06 for ITO
    call system('cp -f g0xrun_p '//trim(address0))
!    call system('cp g0xrun_p '//trim(addresstmp))
! End kuwahata 2021/06/06 for ITO
   !If(NGenGau==1) Then
   !   call system('cp gauss.bss '//trim(addresstmp))
   !EndIf
   open(igauss+id,file=trim(addresstmp)//'gauss.tmp1',status='unknown')
     write(igauss+id,'(a)') '%Chk='//trim(addresstmp)//'gauss.chk'
     write(igauss+id,'(a)') '%RWF='//trim(addresstmp)
     write(igauss+id,'(a)') '%Int='//trim(addresstmp)
     write(igauss+id,'(a)') '%D2E='//trim(addresstmp)
   close(igauss+id)
   call system('cat '//trim(addresstmp)//'gauss.tmp >> '//trim(addresstmp)//'gauss.tmp1')

!! Udagawa Start 2021.05.24 --->
!     If(istepsv == 0 .OR. (nRestart ==1 .AND. istepsv == nrstep+1)) then
!       call system ('sed -e "s/[Gg][Uu][Ee][Ss][Ss]=[Rr][Ee][Aa][Dd]//g" '//trim(addresstmp)//'gauss.tmp1  &
!         & > '//trim(addresstmp)//'gauss.tmp2')
!       call system ('mv '//trim(addresstmp)//'gauss.tmp2 '//trim(addresstmp)//'gauss.tmp1')
!     Endif
!! <--- Udagawa End 2021.05.24

enddo

Return
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
    call program_abort('')
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
    call program_abort('')
  end subroutine err_set_g0xrun

End Subroutine Set_Gaussian_MPI_tk

