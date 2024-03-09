Subroutine Force_New_MPI_tk
  use Parameters
  use utility, only: program_abort
  use mpi_module
  implicit none

#ifdef _mpi_
  !call MyMPI_bcast_r
  call MyMPI_scatter_r
#endif
  select case(Iforce)
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
    case(11)
      call Force_model_Morse
    case(12)
      call Force_model_DoubleWell1D
    case(13)
      call Force_model_DoubleWell3D
    !case(14)
    !  call Force_DoubleHarmonic
    case(15)
      call Force_Harmonic
    !case(16)
    !  call Force_Harmonic_single
    case(21)
      call force_nnp_araidai
    case(31)
      call force_spcf
    case default
      stop 'ERROR!!! Wrong "Iforce" option'
  end select

#ifdef _mpi_
  call MyMPI_gather_fr
#endif
! +++ Start Writting output +++
  call print_result_qm
! +++ End Writting output +++


  if (umbrella_sampling > 0) call calc_umbrella

  Return
End Subroutine
