subroutine set_Iforce
  use Parameters
  use utility, only: program_abort
  implicit none

  select case(Iforce)
    case(1)
      call Set_mopac
    case(6)
      call Set_Gaussian_MPI_tk  ! Set for chk, rwf, etc.
    case(8)
      call Set_VASP
    case(9)
      call Set_siesta
    !case(10)
    !  call Set_Molcas_MPI
    case(21)
      call set_nnp_araidai
    case(22)
      call set_nnp_matlantis
  end select

  return

contains


  subroutine Set_siesta
    use Parameters
    implicit none
    integer   :: i,j,k,imode

    do imode=ista,iend
      write(addresstmp(laddress+1:laddress+6),'(i5.5,A1)') imode,'/'
      call system('mkdir -p '//trim(addresstmp))
      call system('cp InputFile/* '//trim(addresstmp))
    enddo

  return
  end subroutine Set_siesta

end subroutine set_Iforce
