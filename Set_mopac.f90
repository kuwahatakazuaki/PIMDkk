subroutine Set_mopac
  use MPI
  use Parameters
  use utility, only: program_abort
  Implicit None
  Integer   :: i,j,k, imode
  integer :: access

if (MyRank == 0) then
  if      ( access("./mopac.mop", " ") .ne. 0 ) then
    call err_set_mopac_mop
  else if ( access("./g0xrun_p",  " ") .ne. 0 ) then
    call err_set_mopac_g0x
  else
    print *, '"mopac.mop" and "g0xrun_p" are exist -> OK!!'
  end if
end if


do imode=ista,iend
   write(addresstmp(laddress+1:laddress+6),'(i5.5,A1)') imode,'/'
   call system('mkdir -p '//trim(addresstmp))
   call system('cp mopac.mop '//trim(addresstmp))
!   call system('cp g0xrun_p  '//trim(addresstmp))
enddo
!  call program_abort("HELLO")

return

contains

  subroutine err_set_mopac_mop
    integer :: Nunit
    open(newunit=Nunit,file='mopac.mop',status='new')
      write(Nunit,'(a)') "PM7 1SCF XYZ PRECISE GRADIENTS NOSYM NOREOR GEO-OK"
      write(Nunit,'(a)') "comment1"
      write(Nunit,'(a)') "comment2"
    close(Nunit)
    print *, 'ERROR!! "mopac.tmp" is NOT exist!'
    print *, 'Please use "mopac.tmp" template'
    call program_abort('')
  end subroutine err_set_mopac_mop

  subroutine err_set_mopac_g0x
    integer :: Nunit
    open(newunit=Nunit,file='g0xrun_p',status='new')
      write(Nunit,'(a)') "#!/bin/bash"
      write(Nunit,'(a)') "~/mopacCent5/MOPAC2016.exe $1 2> /dev/null"
    close(Nunit)
    print *, 'ERROR!! "g0xrun_p" is NOT exist!'
    print *, 'Please use "g0xrun_p" template'
    call program_abort('')
  end subroutine err_set_mopac_g0x

end subroutine Set_mopac




