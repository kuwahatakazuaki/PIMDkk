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
end subroutine set_Iforce
