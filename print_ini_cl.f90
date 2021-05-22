subroutine print_ini_cl
use Parameters
implicit none
integer :: i
character(len=5) :: cfmt


! Kuwahata 2019/11/05
!   Followings are moved from 'Print_Ham_Classical.f90'
Open(iham,file=trim(address)//'/ham.dat',status='unknown',form='formatted',position='append')
  Write(iham,'(a)') ' ****************************************************&
                    &**********************************************************************'
  Write(iham,'(a)') '   Step      Hamiltonian             Potential             &
                    &  DKinetic              EBath_Cent             Temperature '
  Write(iham,'(a)') ' ****************************************************&
                    &**********************************************************************'
close(iham)

  Write(*,'(a)') ' ************************************************************************************************'
  Write(*,'(a)') '    Step     Hamiltonian      Potential        DKinetic        EBath_Cent      Temperature '
  Write(*,'(a)') ' ************************************************************************************************'
! End Kuwahata 2019/11/05

  if(nocharge==0) then
    open(igetc,file=trim(address)//'/charge.dat',status='unknown',form='formatted',position='append')
!      write(cfmt,'(I5)') Natom
!      write(igetc,'(a,'//trim(cfmt)//'a10)') "# ", alabel(:)
      write(igetc,'(a)',advance='no') "#"
      do i = 1, Natom
        write(igetc,'(a7,I2,1x)',advance='no') alabel(i), i
      end do
    close(igetc)
  endif
return
end subroutine print_ini_cl

