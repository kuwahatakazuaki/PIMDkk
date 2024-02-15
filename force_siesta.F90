subroutine force_siesta
  use Parameters
  !use MPI

  implicit none
  character(len=120) :: line
  character(len=:), allocatable :: out_FA, out_log
  Integer            :: iline, id,imode2,iatom2, index1
  Double Precision   :: enetemp
  Integer :: i,j,k, dummy, igauss = 20

  out_FA = 'siesta.FA'
  out_log= 'log.out'


  id=0
  Call Start_Recv_Send_MPI_tk

   do imode2=ista,iend
     write(addresstmp(laddress+1:laddress+6),'(i5.5,A1)') imode2,'/'

!     open(igauss+id,file=trim(addresstmp)//'temp.xyz',status='unknown')
     call system('cp '//trim(addresstmp)//'input '//trim(addresstmp)//'input.fdf')
     open(igauss+id,file=trim(addresstmp)//'input.fdf',status='old', position='append')
       write(igauss+id,'("AtomicCoordinatesFormat  Ang")')
       write(igauss+id,'("%block AtomicCoordinatesAndAtomicSpecies")')
       do i = 1, Natom
         write(igauss+id,*) r(:,i,imode2)*bohr_inv
         !write(igauss+id,*) x(i,imode2)*bohr_inv, y(i,imode2)*bohr_inv, z(i,imode2)*bohr_inv, atom_num(i)
       end do
       write(igauss+id,'("%endblock AtomicCoordinatesAndAtomicSpecies")')
     close(igauss+id)

!     call system('cat '//trim(addresstmp)//'input1   >  '//trim(addresstmp)//'input.fdf')
!     call system('cat '//trim(addresstmp)//'temp.xyz >> '//trim(addresstmp)//'input.fdf')
!     call system('cat '//trim(addresstmp)//'input2   >> '//trim(addresstmp)//'input.fdf')

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
         read(igauss+id,*) dummy, fr(:,i,imode2)
         !read(igauss+id,*) dummy, fx(i,imode2), fy(i,imode2), fz(i,imode2)
       end do
     close(igauss+id)
!  +++ End Reading "Atomic Force" +++

     call system('rm -rf '//trim(addresstmp)//'input.fdf')
     call system('rm -rf '//trim(addresstmp)//'INPUT_TMP*')
     call system('rm -rf '//trim(addresstmp)//'*log')

     !fx(:,imode2) = fx(:,imode2) * dp_inv * eVtoAU * AUtoAng
     !fy(:,imode2) = fy(:,imode2) * dp_inv * eVtoAU * AUtoAng
     !fz(:,imode2) = fz(:,imode2) * dp_inv * eVtoAU * AUtoAng
     fr(:,:,imode2) = fr(:,:,imode2) * dp_inv * eVtoAU * AUtoAng
     Eenergy(imode2) = Eenergy(imode2) * eVtoAU
   enddo


Call Start_Send_Recv_MPI_tk
call print_result_qm


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


