Subroutine Force_New_MPI_tk
  Use Parameters

  select case(NForce)
    case(1)
      call Force_MOPAC_MPI
    case(2)
      goto 900 !Call Force_Gaussian_MPI
    case(3)
      call Force_DFTB_MPI
!    case(4)
!      call Force_Turbomole_MPI
!    case(5)
!      call Force_Gamess_MPI_tk
    case(6)
      call Force_Gaussian_MPI_tk
!    case(7)
!      call Force_Turbomole_MPI_tk
    case(8)
      call Force_VASP_MPI
    case(9)
      call force_siesta
    case(11)
      call Force_model_Morse
    case(12)
      call Force_model_DoubleWell
    case default
      stop 'ERROR!!! Wrong "NForce" option'
  end select


  Return
  900 continue
  print *, 'We no longer use "Nforce = 2"'
  print *, 'Please use "Nforce = 6"'
  stop
End Subroutine
