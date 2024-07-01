subroutine Force_New_MPI_tk
  use Parameters
  use utility, only: program_abort
#ifdef _mpi_
  use mpi_module, only: MyMPI_gather_fr, MyMPI_gather_others, MyMPI_scatter_r
#endif
  implicit none
  integer :: i

#ifdef _mpi_
  call MyMPI_scatter_r
#endif

  select case(Iforce)

! === Start Ab initio calculation ===
    case(1)
      call Force_MOPAC_MPI
    case(6)
      call Force_Gaussian_MPI_tk
    case(8)
      call Force_VASP_MPI
    case(9)
      call force_siesta
    !case(10)
    !  call Force_Molcas_MPI
! === End Ab initio calculation ===


! === Start model calculation ===
    case(11)
      call Force_model_Morse
    !case(12)
    !  call Force_model_DoubleWell1D
    !case(13)
    !  call Force_model_DoubleWell3D
    case(15)
      call Force_Harmonic
    case(16)
      call Force_Double_Morse
! === End model calculation ===

    case(21)
      call force_nnp_araidai
    !case(22)
    !  call force_nnp_aenet
    case(31)
      call force_spcf
    case default
      stop 'ERROR!!! Wrong "Iforce" option'
  end select

  if (Iumb > 0) call calc_umbrella

#ifdef _mpi_
  call MyMPI_gather_fr
  call MyMPI_gather_others
#endif
! +++ Start Writting output +++
  call print_result_qm
! +++ End Writting output +++
  potential = sum(Eenergy(:)) * dp_inv

  return
end subroutine Force_New_MPI_tk
