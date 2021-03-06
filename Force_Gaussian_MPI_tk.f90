subroutine Force_Gaussian_MPI_tk

!$  use omp_lib
  use Parameters
  use Parameter_tk
  use MPI
  implicit none

  Character (Len=90) :: key1, key2, key3, key4, key5 !, line
  character(len=120) :: line
  Integer            :: iline, id,imode2,iatom2
  Integer            :: index1
  Double Precision   :: enetemp
  Integer :: i,j,k

  if (trim(version_gaussian) == 'g16') then
    select case (theory)
      case(0)
        key1  = ('SCF Done')
      case(1)
        key1  = ('EUMP2')
      case default
        stop 'ERROR!!! Wrong "theory" option'
    end select

    key2  = ('Number     Number              X              Y              Z')
    key3  = ('Dipole moment')
    key4  = ('Mulliken charges')
    key5  = ('Isotropic Fermi Contact Couplings')
  else
    print *, '"version_gaussian" is ', version_gaussian
    stop 'ERROR!!! wrog keyword of "version_gaussian"'
  end if

  id=0
  Call Start_Recv_Send_MPI_tk
!  Call Recv_Send_MPI_tk2

!$omp parallel private(line,iatom2,id,iline,enetemp)
!$omp do
   do imode2=ista,iend
     write(addresstmp(laddress+1:laddress+6),'(i5.5,A1)') imode2,'/'


     call system('cat '//trim(addresstmp)//'gauss.tmp1 > '//trim(addresstmp)//'gauss.com')


!     open(igauss+id,file=trim(addresstmp)//'gauss.xyz',status='unknown')
     open(igauss+id,file=trim(addresstmp)//'gauss.com',status='old',position='append')
       do iatom2=1,natom
!         write(igauss+id,9999) alabel(iatom2),x(iatom2,imode2)*bohr_inv, &
!                            y(iatom2,imode2)*bohr_inv,z(iatom2,imode2)*bohr_inv
         write(igauss+id,*) alabel(iatom2),x(iatom2,imode2)*bohr_inv, &
                            y(iatom2,imode2)*bohr_inv,z(iatom2,imode2)*bohr_inv
       enddo
       write(igauss+id,*)
     close(igauss+id)

!     call system('cat '//trim(addresstmp)//'gauss.xyz >> '//trim(addresstmp)//'gauss.com')

     If(NGenGau==1) Then
        call system('cat '//trim(addresstmp)//'gauss.bss >> '//trim(addresstmp)//'gauss.com')
     EndIf

     call system(trim(address0)//'g0xrun_p '//trim(addresstmp)//'gauss.com '//trim(addresstmp)//'gauss.log '//trim(addresstmp))

     open(igauss+id,file=trim(addresstmp)//'gauss.log')

!  +++ Reading "SCF Done" +++
     do
       read(igauss+id,'(a)',end=401) line
       iline=index(line,trim(key1))  ! Reading "SCE Done"
       if(iline > 0) exit
     end do

     if (theory == 0) then
       index1 = index(line,'=')
       read(line(index1+2:index1+17),*) enetemp ! For DFT
     else if (theory == 1) then
       read(line(38:60),*) enetemp ! For MP2
     end if

     Eenergy(imode2)=enetemp
!  +++ End Reading "SCF Done" +++

!  +++ Reading "Mulliken charge" +++
     if(nocharge==0) then
!       rewind(igauss+id)
       do
         read(igauss+id,'(a)',end=404) line
         iline=index(trim(line),trim(key4))  ! Reading "Mulliken charge"
         if(iline > 0) exit
       end do
       read(igauss+id,*)
       do iatom2=1,natom
          Read(igauss+id,'(a)') line
          Read(line(11:21),*) charge(iatom2,imode2)
       enddo
       if(nohfcc==0) then
         do
           Read(igauss+id,'(a)',end=405) line
           iline=index(trim(line),trim(key5))
           if(iline > 0) exit
         end do
         read(igauss+id,'()')
         do iatom2=1,natom
           Read(igauss+id,'(a)') line
           Read(line(36:49),*) hfcc(iatom2,imode2)
         enddo
       endif
     endif
!  +++ End Reading "Mulliken charge" +++

!  +++ Reading "Dipole moment" +++
!     Rewind(igauss+id)
     if(nodipole==0) then
       do
         read(igauss+id,'(a)',end=403) line
         iline=index(trim(line),trim(key3))  ! Reading "Dipole moment"
         if(iline > 0) exit
       end do
       Read(igauss+id,'(a)') line
       Read(line(20:26),*) dipolex(imode2)
       Read(line(46:52),*) dipoley(imode2)
       Read(line(72:78),*) dipolez(imode2)
       Read(line(98:104),*) dipole(imode2)
     endif
!  +++ End Reading "Dipole moment" +++

!  +++ Reading "Atomic Force" +++
!     rewind(igauss+id)
     do
       read(igauss+id,'(a)',end=402) line  ! Reading "Force"
!       iline=index(trim(line),trim(key2))
       iline=index(line,trim(key2))
       if(iline > 0) exit
     end do
     read(igauss+id,*)
     do iatom2=1,natom
        read(igauss+id,'(a)') line
        read(line(24:38),*) fx(iatom2,imode2)
        read(line(39:53),*) fy(iatom2,imode2)
        read(line(54:68),*) fz(iatom2,imode2)
     enddo
!  +++ End Reading "Atomic Force" +++


     close(igauss+id)

     call system('rm -rf '//trim(addresstmp)//'Gau*')

     fx(:,imode2) = fx(:,imode2) * dp_inv
     fy(:,imode2) = fy(:,imode2) * dp_inv
     fz(:,imode2) = fz(:,imode2) * dp_inv
   enddo

!$omp end do
!$omp end parallel

Call Start_Send_Recv_MPI_tk
! Call Send_Recv_MPI_tk2

! Kuwahata 2020/01/12
! Call Write_MPI_tk
! if (MyRank == 0)
call print_gaussian
! End Kuwahata 2020/01/12


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
end subroutine

