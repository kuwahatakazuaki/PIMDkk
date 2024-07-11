subroutine Set_VASP
!$  use omp_lib
  use Parameters
  use utility, only: program_abort
  implicit none
  integer   :: i,j,k,imode
  integer :: access

if ( MyRank == 0 ) then
  if ( access("./LATTICE", " ")     /= 0)  call err_LATTICE
  if ( access("./vasp_run.sh", " ") /= 0)  call err_exefile
end if

do imode=ista,iend
   write(addresstmp(laddress+1:laddress+6),'(i5.5,A1)') imode,'/'

   call system('mkdir -p '//trim(addresstmp))  ! need ??
   call system('cp INCAR POTCAR KPOINTS '//trim(addresstmp))
   call system('cp LATTICE '//trim(addresstmp)//'POSCAR')

! Machine learning version
   call system('ls ML_AB ML_FF 2>/dev/null && cp ML_AB ML_FF '//trim(addresstmp))
! Machine learning version
!   call system('cp vasp_run.sh '//trim(addresstmp))
enddo

return
contains

  subroutine err_LATTICE
    integer :: Nunit
    open(newunit=Nunit,file='LATTICE',status='new')
      write(Nunit,'(a)') "H2O molecule                                                  "
      write(Nunit,'(a)') "1.0                                                           "
      write(Nunit,'(a)') "        7.9376997948         0.0000000000         0.0000000000"
      write(Nunit,'(a)') "        0.0000000000         7.9376997948         0.0000000000"
      write(Nunit,'(a)') "        0.0000000000         0.0000000000         7.9376997948"
      write(Nunit,'(a)') "    O    H                                                    "
      write(Nunit,'(a)') "    1    2                                                    "
    close(Nunit)
    print *, 'ERROR!! "LATTICE" is NOT exist!'
    print *, 'Please use "LATTICE" template'
    call program_abort('')
  end subroutine err_LATTICE

  subroutine err_exefile
    integer :: Nunit
    open(newunit=Nunit,file='vasp_run.sh',status='new')
      write(Nunit,'(a)') '#!/bin/bash'
      write(Nunit,'(a)') 'cd $1                                                         '
      write(Nunit,'(a)') 'exe_vasp="/home/vasp_guest/.local/vasp.5.4.4.intel_serial/bin"'
      write(Nunit,'(a)') '$exe_vasp/vasp_std 1> log.out 2> err.out                      '
    close(Nunit)
    print *, 'ERROR!! "vasp_run.sh" is NOT exist!'
    print *, 'Please use "vasp_run.sh" template'
    call program_abort('')
  end subroutine err_exefile

end subroutine Set_VASP
