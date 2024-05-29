subroutine set_nnp_araidai
  use Parameters
  implicit none
  Integer   :: i,j,k
  integer :: imode

do imode = Ista, Iend
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
    only: Eenergy, r, fr, Natom, AUtoAng, eVtoAU, dp_inv, alabel, &
          addresstmp, Ista, Iend, laddress
  implicit none

  character(Len=130) :: line
  integer :: imode, i, Uout, Uinp
  character :: dummyC
  !real(8), parameter :: eVtoAU  = 1.0d0 / 27.21138505
  real(8), parameter :: eVAng_HartBohr = 0.5291772108d0 / 27.21138505d0

  do imode = Ista, Iend
    write(addresstmp(laddress+1:laddress+6),'(i5.5,A1)') imode,'/'

    open(newunit=Uout,file=trim(addresstmp)//'training.data_1',status='replace',position='append')
      write(Uout,*) "begin"
      write(Uout,*) "comment   LiNiO2_633"
      write(Uout,*) "lattice   5.7560070000E+00  0.0000000000E+00  0.0000000000E+00"
      write(Uout,*) "lattice   0.0000000000E+00  4.9848490000E+00  0.0000000000E+00"
      write(Uout,*) "lattice   0.0000000000E+00  0.0000000000E+00  1.4190003000E+01"
      do i=1,Natom
        write(Uout,9998) &
         "atom", r(:,i,imode)*AUtoAng, alabel(i), 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0
      enddo
      write(Uout,*) "energy   -2.0683566596E+02"
      write(Uout,*) "charge    0.0000000000E+00"
      write(Uout,*) "end"
    close(Uout)

    call system('cd '//trim(addresstmp)//' ; mpiexec.hydra -n 1 ./n2training ; cd ../.. ')

    open(newunit=Uinp,file=trim(addresstmp)//'foce_npp.out',status='old')
      read(Uinp,*) Eenergy(imode)
      do i = 1, Natom
        read(Uinp,*) dummyC, fr(:,i,imode)
        !read(Uinp,*) dummyC, fx(i,imode), fy(i,imode), fz(i,imode)
      end do
    close(Uinp)

    fr(:,:,imode)=fr(:,:,imode)*eVAng_HartBohr*dp_inv
    Eenergy(imode) = Eenergy(imode) * eVtoAU
  end do


return
9999 format(3F24.16)
9998 format(A, 3E17.9, x, A, 6E12.4)
end subroutine force_nnp_araidai

subroutine force_nnp_aenet
end subroutine force_nnp_aenet


