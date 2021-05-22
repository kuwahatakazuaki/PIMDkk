subroutine force_siesta

  use Parameters
  use Parameter_tk
  use MPI
  implicit none

!  Character (Len=90) :: key1, key2, key3, key4, key5 !, line
  character(len=120) :: line
  character(len=:), allocatable :: out_FA, out_log
  Integer            :: iline, id,imode2,iatom2, index1
!  Integer            :: index1
  Double Precision   :: enetemp
  Integer :: i,j,k, dummy
  integer :: atom_num(3) = [1,2,2]
  real(8), parameter :: eVtoAU    = 1.0d0/27.21162
  real(8), parameter :: AngtoAU = bohr
  real(8), parameter :: AUtoAng = bohr_inv

  out_FA = 'H2O'//'.FA'
  out_log= 'log.out'


  id=0
  Call Start_Recv_Send_MPI_tk

   do imode2=ista,iend
     write(addresstmp(laddress+1:laddress+6),'(i5.5,A1)') imode2,'/'

     open(igauss+id,file=trim(addresstmp)//'temp.xyz',status='unknown')
!       do iatom2=1,natom
!         write(igauss+id,9999) alabel(iatom2),x(iatom2,imode2)*bohr_inv, &
!                            y(iatom2,imode2)*bohr_inv,z(iatom2,imode2)*bohr_inv
!       enddo
       do i = 1, Natom
         write(igauss+id,*) x(i,imode2)*bohr_inv, y(i,imode2)*bohr_inv, z(i,imode2)*bohr_inv, atom_num(i)
       end do
     close(igauss+id)

     call system('cat '//trim(addresstmp)//'input1   >  '//trim(addresstmp)//'input.fdf')
     call system('cat '//trim(addresstmp)//'temp.xyz >> '//trim(addresstmp)//'input.fdf')
     call system('cat '//trim(addresstmp)//'input2   >> '//trim(addresstmp)//'input.fdf')

     if (NGenGau == 1) call system('cat '//trim(addresstmp)//'gauss.bss >> '//trim(addresstmp)//'input.fdf')

!     call system(trim(address0)//'g0xrun_p '//trim(addresstmp)//'input.fdf '//trim(addresstmp)//'gauss.log '//trim(addresstmp))
     call system(trim(addresstmp)//'/run_siesta.sh '//trim(addresstmp))

     open(igauss+id,file=trim(addresstmp)//out_log)

!  +++ Reading "SCF Done" +++
       do
         read(igauss+id,'(a)',end=401) line
         if ( index(line,'siesta: E_KS(eV)') > 0 ) exit
       end do
       read(line(31:41),*) Eenergy(imode2)
     close(igauss+id)
!  +++ End Reading "SCF Done" +++

!  +++ Reading "Atomic Force" +++
!     rewind(igauss+id)
     open(igauss+id,file=trim(addresstmp)//out_FA)
       read(igauss+id,'()')
       do i = 1, Natom
         read(igauss+id,*) dummy, fx(i,imode2), fy(i,imode2), fz(i,imode2)
       end do
     close(igauss+id)
!  +++ End Reading "Atomic Force" +++

     call system('rm -rf '//trim(addresstmp)//'input.fdf')
     call system('rm -rf '//trim(addresstmp)//'INPUT_TMP*')
     call system('rm -rf '//trim(addresstmp)//'*log')

     fx(:,imode2) = fx(:,imode2) * dp_inv * eVtoAU * AUtoAng
     fy(:,imode2) = fy(:,imode2) * dp_inv * eVtoAU * AUtoAng
     fz(:,imode2) = fz(:,imode2) * dp_inv * eVtoAU * AUtoAng
     Eenergy(imode2) = Eenergy(imode2) * eVtoAU
   enddo


Call Start_Send_Recv_MPI_tk

! call print_gaussian
call print_siesta


9999 format(a2,1x,E15.9,1x,E15.9,1x,E15.9)
9998 format(3E23.15)
9997 format(2E23.15)
9996 format(E23.15)
9995 format(4E23.15)

return
401 print *, 'ERROR!!: We can not find "SCF Done" or "EUMP2" in Gaussian output'; stop
402 print *, 'ERROR!!: We can not find "Force" in Gaussian output'; stop
403 print *, 'ERROR!!: We can not find "Dipole moment" in Gaussian output'; stop
404 print *, 'ERROR!!: We can not find "Mulliken charges" in Gaussian output'; stop
405 print *, 'ERROR!!: We can not find "Fermi Contact" in Gaussian output'; stop
end subroutine force_siesta


!     if (theory == 0) then
!       index1 = index(line,'=')
!       read(line(index1+2:index1+17),*) enetemp ! For DFT
!     else if (theory == 1) then
!       read(line(38:60),*) enetemp ! For MP2
!     end if

!     Eenergy(imode2)=enetemp
!       do
!         read(igauss+id,'(a)',end=402) line  ! Reading "Force"
!         iline=index(line,trim(key2))
!         if(iline > 0) exit
!       end do
!       read(igauss+id,*)
!       do iatom2=1,natom
!          read(igauss+id,'(a)') line
!          read(line(24:38),*) fx(iatom2,imode2)
!          read(line(39:53),*) fy(iatom2,imode2)
!          read(line(54:68),*) fz(iatom2,imode2)
!       enddo
