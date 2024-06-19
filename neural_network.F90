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
  use utility, only: program_abort

  implicit none
  character(Len=130) :: line, inp_train(8)
  integer :: imode, i, Uout, Uinp, ios
  character :: dummyC
  real(8), parameter :: eVAng_HartBohr = 0.5291772108d0 / 27.21138505d0

  open(newunit=Uinp,file='training.data_1',status='old',iostat=ios)
    if ( ios /= 0) then
      call program_abort('ERROR!!: There is no file of "training.data_1"')
    end if
    do i = 1, 5
      read(Uinp,'(a)') inp_train(i)
    end do
    do i = 1, Natom
      read(Uinp,'()')
    end do
    do i = 6, 8
      read(Uinp,'(a)') inp_train(i)
    end do
  close(Uinp)

!do i = 1, 8
!  print *, inp_train(i)
!end do

  do imode = Ista, Iend
    write(addresstmp(laddress+1:laddress+6),'(i5.5,A1)') imode,'/'

    open(newunit=Uout,file=trim(addresstmp)//'training.data_1',status='replace')
      do i = 1, 5
        write(Uout,'(a)') trim(inp_train(i))
      end do
      do i=1,Natom
        write(Uout,9998) &
         "atom", r(:,i,imode)*AUtoAng, alabel(i), 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0
      enddo
      do i = 6, 8
        write(Uout,'(a)') trim(inp_train(i))
      end do
    close(Uout)

    call system('cd '//trim(addresstmp)//' ; ./n2training > /dev/null ; cd ../.. ')
    !call system('cd '//trim(addresstmp)//' ; mpiexec.hydra -n 1 ./n2training ; cd ../.. ')

    open(newunit=Uinp,file=trim(addresstmp)//'foce_npp.out',status='old')
      read(Uinp,*) Eenergy(imode)
      do i = 1, Natom
        read(Uinp,*) dummyC, fr(:,i,imode)
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


