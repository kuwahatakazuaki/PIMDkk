subroutine Set_siesta
  use Parameters
  Implicit None
  Integer   :: i,j,k,id,imode

  id=0
do imode=ista,iend
   write(addresstmp(laddress+1:laddress+6),'(i5.5,A1)') imode,'/'
   call system('mkdir -p '//trim(addresstmp))
   call system('cp InputFile/* '//trim(addresstmp))
enddo

return
end subroutine

