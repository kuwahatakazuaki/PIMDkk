Subroutine Set_Gaussian_MPI_tk

!$  use omp_lib
  Use MPI
  Use Parameters
  Use Parameter_tk
  Implicit None
  Integer   :: i,j,k,id

  id=0
!$omp parallel private(addresstmp,line,iatom2,id,iline,enetemp)
!$omp do
do imode=ista,iend
!$omp  id=omp_get_thread_num()
   write(addresstmp(laddress+1:laddress+6),'(i5.5,A1)') imode,'/'
   call system('mkdir -p '//trim(addresstmp))
!   wait(i)
   call system('cp gauss.tmp '//trim(addresstmp))
   If(NGenGau==1) Then
      call system('cp gauss.bss '//trim(addresstmp))
   EndIf
   open(igauss+id,file=trim(addresstmp)//'gauss.tmp1',status='unknown')
     write(igauss+id,'(a)') '%Chk='//trim(addresstmp)//'gauss.chk'
     write(igauss+id,'(a)') '%RWF='//trim(addresstmp)
     write(igauss+id,'(a)') '%Int='//trim(addresstmp)
     write(igauss+id,'(a)') '%D2E='//trim(addresstmp)
   close(igauss+id)
   call system('cat '//trim(addresstmp)//'gauss.tmp >> '//trim(addresstmp)//'gauss.tmp1')
enddo
!$omp end do
!$omp end parallel

Return
End Subroutine

