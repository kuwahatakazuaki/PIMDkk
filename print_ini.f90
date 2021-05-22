subroutine print_ini
use global_variable
implicit none
integer :: iham, igetc
integer :: i
character(len=5) :: cfmt


if ( Nbead == 1) then
  call print_ham_cl
else if ( Nbead >= 2) then
  call print_ham_qm
end if


if( Lcharge .eqv. .True. ) then
  open(newunit=igetc,file=trim(path_result)//'/charge.dat',status='unknown',form='formatted',position='append')
!      write(cfmt,'(I5)') Natom
!      write(igetc,'(a,'//trim(cfmt)//'a10)') "# ", alabel(:)
    write(igetc,'(a)',advance='no') "#"
    do i = 1, Natom
      write(igetc,'(a7,I2,1x)',advance='no') alabel(i), i
    end do
  close(igetc)
endif


return

contains

  subroutine print_ham_qm
    open(newunit=iham,file=trim(path_result)//'/ham.dat',status='unknown',form='formatted',position='append')
    write(iham,'("#  ",a)')  repeat('*',188)
    write(iham,'(a)') &
      &'#   Step      Hamiltonian            Potential               DKinetic                &
       QKinetic               EBath               EBath_Cent              Temperature            E_Virial   '
    write(iham,'("#  ",a)')  repeat('*',188)
!    write(iham,'(a)') &
!      &'#  **************************************************************************************&
!      &*****************************************************************************************************'
    close(iham)

    write(*,'(a)')  repeat('*',143)
    write(*,'(a)') &
      &'   Step    Hamiltonian       Potential         DKinetic         &
      QKinetic           EBath          EBath_Cent       Temperature      E_Virial '
    write(*,'(a)')  repeat('*',143)
!    write(*,'(a)') &
!      &'******************************************************************************&
!      &*****************************************************************'
   close(iham)
  end subroutine print_ham_qm


  subroutine print_ham_cl
    open(newunit=iham,file=trim(path_result)//'/ham.dat',status='unknown',form='formatted',position='append')
      write(iham,'("# ",a)') repeat('*',121)
!      write(iham,'(a)') ' ****************************************************&
!                        &**********************************************************************'
      write(iham,'(a)') '#   Step     Hamiltonian             Potential             &
                        &  DKinetic              EBath_Cent             Temperature '
      write(iham,'("# ",a)') repeat('*',121)
    close(iham)

    write(*,'(a)') ' ************************************************************************************************'
    write(*,'(a)') '    Step     Hamiltonian      Potential        DKinetic        EBath_Cent      Temperature '
    write(*,'(a)') ' ************************************************************************************************'
  end subroutine print_ham_cl

end subroutine print_ini


!subroutine print_ini_cl
!use global_variable
!implicit none
!integer :: iham, igetc
!integer :: i
!character(len=5) :: cfmt
!
!  if( Lcharge .eqv. .True. ) then
!    open(newunit=igetc,file=trim(path_result)//'/charge.dat',status='unknown',form='formatted',position='append')
!!      write(cfmt,'(I5)') Natom
!!      write(igetc,'(a,'//trim(cfmt)//'a10)') "# ", alabel(:)
!      write(igetc,'(a)',advance='no') "#"
!      do i = 1, Natom
!        write(igetc,'(a7,I2,1x)',advance='no') alabel(i), i
!      end do
!    close(igetc)
!  endif
!return
!end subroutine print_ini_cl

