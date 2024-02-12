Subroutine Set_VASP

!$  use omp_lib
  Use MPI
  Use Parameters
  use utility, only: program_abort
  Implicit None
  Integer   :: i,j,k,id, imode
  integer :: access

if ( MyRank == 0 ) then
  if      ( access("./LATTICE", " ") .ne. 0) then
    call err_LATTICE
  else if ( access("./vasp_run.sh", " ") .ne. 0) then
    call err_set_g0xrun
  end if
end if


  id=0
!$omp parallel private(addresstmp,line,iatom2,id,iline,enetemp)
!$omp do
do imode=ista,iend
!$  id=omp_get_thread_num()
   write(addresstmp(laddress+1:laddress+6),'(i5.5,A1)') imode,'/'

   call system('mkdir -p '//trim(addresstmp))  ! need ??
   call system('cp INCAR POTCAR KPOINTS '//trim(addresstmp))
   call system('cp LATTICE '//trim(addresstmp)//'POSCAR')

! Machine learning version
   call system('ls ML_AB ML_FF 2>/dev/null && cp ML_AB ML_FF '//trim(addresstmp))
! Machine learning version
!   call system('cp vasp_run.sh '//trim(addresstmp))

enddo
!$omp end do
!$omp end parallel

Return
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

  subroutine err_set_g0xrun
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
  end subroutine err_set_g0xrun

End Subroutine Set_VASP

