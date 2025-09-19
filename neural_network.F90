subroutine force_nnp_matlantis
!!! === Nproc much be 1 at current verion === !!!
  use Parameters, &
    only: pot_bead, r, fr, Natom, AU2Ang, eV2AU, dp_inv, alabel, &
          addresstmp, Ista, Iend, laddress, MyRank, Nbead, eVAng2AU, &
          Lperiodic
  use utility, only: program_abort
  implicit none
  integer :: Nfile, Imode, i, j
  integer :: Uinp, Uout, ios
  character(len=13), allocatable :: Fout(:)
  character(len=12) :: char_num
  character(len=:), allocatable :: command, Fforce, Fenergy
  integer :: access, Iaccess

  Fforce  = 'forces.out'
  Fenergy = 'energy.out'

  Nfile = Iend - Ista + 1
  allocate(Fout(Nfile))
  if ( Lperiodic .eqv. .True. ) then
    call input_periodic
  else
    call input_nonperio
  end if

  if (access(Fforce, ' ') == 0) call system('rm -f '//Fforce)
  if (access(Fenergy, ' ') == 0) call system('rm -f '//Fenergy)

! +++ Execution Matlantis +++
  write(char_num,'(" ",I0, " ",I0, L2)') Ista, Iend, Lperiodic
  command = ' ./run_matlantis.py '//trim(addresstmp)//trim(char_num)
  call system(command)
! +++ Execution Matlantis +++

  !if ( MyRank == 0 ) then
  open(newunit=Uinp,file=Fforce,status='old',iostat=ios)
    if ( ios /= 0 ) then
      call program_abort('ERROR!!! There is no "'//Fforce//'"')
    end if
    do j = 1, Nbead
      do i = 1, Natom
        read(Uinp,*) fr(:,i,j)
      end do
    end do
  close(Uinp)
  open(newunit=Uinp,file=Fenergy,status='old',iostat=ios)
    if ( ios /= 0 ) then
      call program_abort('There is no "'//Fenergy//'"')
    end if
    do j = 1, Nbead
      read(Uinp,*) pot_bead(j)
    end do
  close(Uinp)

  fr(:,:,:)  = fr(:,:,:) * eVAng2AU * dp_inv
  pot_bead(:) = pot_bead(:) * eV2AU
  !end if

contains

  subroutine input_periodic
    do Imode = Ista, Iend
      write(Fout(Imode),'("str",I5.5,".vasp")') Imode
    end do
    do Imode = Ista, Iend
      call system( 'cp LATTICE '//trim(addresstmp)//trim(Fout(Imode)) )
      open(newunit=Uout,file=trim(addresstmp)//trim(Fout(Imode)),status='old',position='append')
        do i = 1, Natom
          write(Uout,*) r(:,i,Imode) * AU2Ang
        end do
        write(Uout,*)
      close(Uout)
    end do
  end subroutine input_periodic

  subroutine input_nonperio
    do Imode = Ista, Iend
      write(Fout(Imode),'("str",I5.5,".xyz")') Imode
    end do
    do Imode = Ista, Iend
      open( newunit=Uout,file=trim(addresstmp)//trim(Fout(Imode)) )
        write(Uout,*) Natom
        write(Uout,*) 'Properties=species:S:1:pos:R:3 pbc="F F F"'
        do i = 1, Natom
          write(Uout,*) alabel(i), r(:,i,Imode)*AU2Ang
        end do
      close(Uout)
    end do
  end subroutine input_nonperio

end subroutine force_nnp_matlantis

subroutine set_nnp_matlantis
  use Parameters, only: Lperiodic, Nproc
  use utility, only: program_abort
  implicit none
  character(len=:), allocatable :: name_file
  integer :: access

  if ( Nproc /= 1 ) then
    call program_abort('ERROR!!! current version only for Nproc == 1')
  end if

  name_file = 'run_matlantis.py'
  if ( access(name_file,' ') /= 0 ) then
    call program_abort('ERROR!!! There is no "'//name_file//'"')
  end if

  if ( Lperiodic .eqv. .True. ) then
    if ( access("./LATTICE", " ")     /= 0)  call err_LATTICE
  end if

contains

  subroutine err_LATTICE
    integer :: Nunit
    open(newunit=Nunit,file='LATTICE',status='new')
      write(Nunit,'(a)') trim("H2O molecule                                                  ")
      write(Nunit,'(a)') trim("1.0                                                           ")
      write(Nunit,'(a)') trim("        7.9376997948         0.0000000000         0.0000000000")
      write(Nunit,'(a)') trim("        0.0000000000         7.9376997948         0.0000000000")
      write(Nunit,'(a)') trim("        0.0000000000         0.0000000000         7.9376997948")
      write(Nunit,'(a)') trim("    O    H                                                    ")
      write(Nunit,'(a)') trim("    1    2                                                    ")
      write(Nunit,'(a)') trim("Cartesian                                                     ")
    close(Nunit)
    print *, 'ERROR!! "LATTICE" is NOT exist!'
    print *, 'Please use "LATTICE" template'
    call program_abort('')
  end subroutine err_LATTICE
end subroutine set_nnp_matlantis


subroutine set_nnp_araidai
  use Parameters
  implicit none
  Integer   :: i,j,k
  integer :: Imode

  do Imode = Ista, Iend
    write(addresstmp(laddress+1:laddress+6),'(i5.5,A1)') Imode,'/'
    call system('mkdir -p '//trim(addresstmp))
    call system('cp -r nnp_files '//trim(addresstmp))
    call system('cp training.data_1 '//trim(addresstmp))
    call system('cp n2training '//trim(addresstmp))
  enddo
return
end subroutine set_nnp_araidai


subroutine force_nnp_araidai
  use Parameters, &
    only: pot_bead, r, fr, Natom, AU2Ang, eV2AU, dp_inv, alabel, &
          addresstmp, Ista, Iend, laddress, eVAng2AU
  use utility, only: program_abort

  implicit none
  character(Len=130) :: line, inp_train(8)
  integer :: Imode, i, Uout, Uinp, ios
  character :: dummyC
  !real(8), parameter :: eVAng_HartBohr = 0.5291772108d0 / 27.21138505d0

  open(newunit=Uinp,file='training.data_1',status='old',iostat=ios)
    if ( ios /= 0) then
      call program_abort('ERROR!!: There is no "training.data_1"')
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

  do Imode = Ista, Iend
    write(addresstmp(laddress+1:laddress+6),'(i5.5,A1)') Imode,'/'

    open(newunit=Uout,file=trim(addresstmp)//'training.data_1',status='replace')
      do i = 1, 5
        write(Uout,'(a)') trim(inp_train(i))
      end do
      do i=1,Natom
        write(Uout,9998) &
         "atom", r(:,i,Imode)*AU2Ang, alabel(i), 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0
      enddo
      do i = 6, 8
        write(Uout,'(a)') trim(inp_train(i))
      end do
    close(Uout)

    call system('cd '//trim(addresstmp)//' ; ./n2training > /dev/null ; cd ../.. ')
    !call system('cd '//trim(addresstmp)//' ; mpiexec.hydra -n 1 ./n2training ; cd ../.. ')

    open(newunit=Uinp,file=trim(addresstmp)//'foce_npp.out',status='old')
      read(Uinp,*) pot_bead(Imode)
      do i = 1, Natom
        read(Uinp,*) dummyC, fr(:,i,Imode)
      end do
    close(Uinp)

    fr(:,:,Imode)=fr(:,:,Imode)*eVAng2AU*dp_inv
    pot_bead(Imode) = pot_bead(Imode) * eV2AU
  end do

return
9999 format(3F24.16)
9998 format(A, 3E17.9, x, A, 6E12.4)
end subroutine force_nnp_araidai

subroutine force_nnp_aenet
end subroutine force_nnp_aenet


