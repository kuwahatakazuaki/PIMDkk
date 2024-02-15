subroutine set_nnp_araidai
  use Parameters
  implicit none
  Integer   :: i,j,k! ,id
  integer :: imode

do imode=ista,iend
   write(addresstmp(laddress+1:laddress+6),'(i5.5,A1)') imode,'/'
   call system('mkdir -p '//trim(addresstmp))
   call system('cp -r nnp_files '//trim(addresstmp))
   call system('cp training.data_1 '//trim(addresstmp))
   call system('cp n2training '//trim(addresstmp))
enddo

return
end subroutine set_nnp_araidai


subroutine force_nnp_araidai
  use Parameters, &
    !only: Eenergy, x, y, z, fx, fy, fz, Natom, bohr_inv, dp_inv, alabel, &
    only: Eenergy, r, fr, Natom, bohr_inv, dp_inv, alabel, &
          addresstmp, ista, iend, laddress
  implicit none

  character(Len=130) :: line
  integer :: imode, i, Uout, Uinp
  character :: dummyC
  real(8), parameter :: ev_to_hartree  = 1.0 / 27.21138505
  real(8), parameter :: eVAng_HartBohr = 0.5291772108 / 27.21138505

  Call Start_Recv_Send_MPI_tk
  do imode=ista,iend
    write(addresstmp(laddress+1:laddress+6),'(i5.5,A1)') imode,'/' 

    open(newunit=Uout,file=trim(addresstmp)//'training.data_1',status='replace',position='append')
      write(Uout,*) "begin"
      write(Uout,*) "comment   LiNiO2_633"
      write(Uout,*) "lattice   5.7560070000E+00  0.0000000000E+00  0.0000000000E+00"
      write(Uout,*) "lattice   0.0000000000E+00  4.9848490000E+00  0.0000000000E+00"
      write(Uout,*) "lattice   0.0000000000E+00  0.0000000000E+00  1.4190003000E+01"
      do i=1,Natom
        write(Uout,9998) &
         "atom", r(:,i,imode)*bohr_inv, alabel(i), 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0
         !"atom", x(i,imode)*bohr_inv, y(i,imode)*bohr_inv, z(i,imode)*bohr_inv, alabel(i), 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0
      enddo
      write(Uout,*) "energy   -2.0683566596E+02"
      write(Uout,*) "charge    0.0000000000E+00"
      write(Uout,*) "end"
    close(Uout)

    call system('cd '//trim(addresstmp)//' ; mpiexec.hydra -n 1 ./n2training ; cd ../.. ')
    !call execute_command_line( './do_nnp.sh '//trim(addresstmp)  )

    open(newunit=Uinp,file=trim(addresstmp)//'foce_npp.out',status='old')
      read(Uinp,*) Eenergy(imode)
      do i = 1, Natom
        read(Uinp,*) dummyC, fr(:,i,imode)
        !read(Uinp,*) dummyC, fx(i,imode), fy(i,imode), fz(i,imode)
      end do
    close(Uinp)

    !fx(:,imode)=fx(:,imode)*eVAng_HartBohr*dp_inv
    !fy(:,imode)=fy(:,imode)*eVAng_HartBohr*dp_inv
    !fz(:,imode)=fz(:,imode)*eVAng_HartBohr*dp_inv
    fr(:,:,imode)=fr(:,:,imode)*eVAng_HartBohr*dp_inv
    Eenergy(imode) = Eenergy(imode) * ev_to_hartree
  end do

  call Start_Send_Recv_MPI_tk

return
9999 format(3F24.16)
9998 format(A, 3E17.9, x, A, 6E12.4)

end subroutine force_nnp_araidai

