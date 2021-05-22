subroutine setup_indivi
  use global_variable
  implicit none
  select case(Nforce)
    case(6)
      call set_gaussian
    case(9)
      call set_siesta
  end select

contains

  subroutine set_gaussian
    integer :: imode, Ugauss
    character(len=:), allocatable :: path_scr_sub
    character(len=5) :: temp_chr

    do imode = Ista, Iend
      write(temp_chr,'(I5.5)') imode
      path_scr_sub = path_scr//'/'//temp_chr//'/'
      call system('mkdir -p '//path_scr_sub)
      call system('cp gauss.tmp '//path_scr_sub)
      call system('cp g0xrun_p  '//path_scr)
      if (Lgengau .eqv. .True.) call system('cp gauss.bss '//path_scr_sub)
      open(newunit=Ugauss, file=path_scr_sub//'gauss.tmp1')
        write(Ugauss,'(a)') '%Chk='//path_scr_sub
        write(Ugauss,'(a)') '%RWF='//path_scr_sub
        write(Ugauss,'(a)') '%Int='//path_scr_sub
        write(Ugauss,'(a)') '%D2E='//path_scr_sub
      close(Ugauss)
    end do
    call system('cat '//path_scr_sub//'gauss.tmp >> '//path_scr_sub//'gauss.tmp1')
  end subroutine set_gaussian

  subroutine set_siesta
  end subroutine set_siesta

end subroutine setup_indivi

!Subroutine Set_Gaussian_MPI_tk
!
!!$  use omp_lib
!  Use MPI
!  Use Parameters
!  Use Parameter_tk
!  Implicit None
!  Integer   :: i,j,k,id
!
!  id=0
!!$omp parallel private(addresstmp,line,iatom2,id,iline,enetemp)
!!$omp do
!do imode=ista,iend
!!$omp  id=omp_get_thread_num()
!   write(addresstmp(laddress+1:laddress+6),'(i5.5,A1)') imode,'/'
!   call system('mkdir -p '//trim(addresstmp))
!!   wait(i)
!   call system('cp gauss.tmp '//trim(addresstmp))
!   If(NGenGau==1) Then
!      call system('cp gauss.bss '//trim(addresstmp))
!   EndIf
!   open(igauss+id,file=trim(addresstmp)//'gauss.tmp1',status='unknown')
!     write(igauss+id,'(a)') '%Chk='//trim(addresstmp)//'gauss.chk'
!     write(igauss+id,'(a)') '%RWF='//trim(addresstmp)
!     write(igauss+id,'(a)') '%Int='//trim(addresstmp)
!     write(igauss+id,'(a)') '%D2E='//trim(addresstmp)
!   close(igauss+id)
!   call system('cat '//trim(addresstmp)//'gauss.tmp >> '//trim(addresstmp)//'gauss.tmp1')
!enddo
!!$omp end do
!!$omp end parallel
!
!Return
!End Subroutine

