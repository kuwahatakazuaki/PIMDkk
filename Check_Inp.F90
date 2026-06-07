subroutine Check_Inp
  use Parameters
#ifdef _MACE_
  use mace_force_config, only: get_mace_model_path
#endif
  implicit none
  integer :: Iatom, Uout
#ifdef _MACE_
  character(len=256) :: mace_model_path
#endif

  if ( MyRank == 0 ) then
    open(newunit=Uout,file=Fout,status='old',position='append')
      write(Uout,'(a)')' +++++ Input Check +++++   '
      write(Uout,9999) ' +++++ Isimulation         ', Isimulation
      write(Uout,9997) ' +++++ Simulation type     ', trim(name_simulation)
      write(Uout,9999) ' +++++ Number of Atoms     ', Natom
      write(Uout,9999) ' +++++ Number of Dim       ', Ndim
      write(Uout,9999) ' +++++ Number of Beads     ', Nbead
      write(Uout,9999) ' +++++ Number of Steps     ', Nstep
      if ( Isimulation == 3 ) then
        write(Uout,9999) ' +++++ Number of HMC dyn   ', Ndyn
      end if
#ifdef _mpi_
      write(Uout,9999) ' +++++ Number of Nproc     ', Nproc
#endif
      write(Uout,9998) ' +++++ Given Temperature   ', Temperature
      write(Uout,9998) ' +++++ Given time step     ', dt / fs2AU
      write(Uout,9998) ' +++++ Output step (fs)    ', dt / fs2AU * out_step
      write(Uout,9999) ' +++++ Method of Centr NHC ', Ncent
      write(Uout,9999) ' +++++ Length of Centr NHC ', Nnhc
      write(Uout,9999) ' +++++ Flag for Force Calc ', Iforce
      if ( Isimulation == 3 ) then
        write(Uout,9995) ' +++++ Flag for Dual HMC   ', Ldual
        if ( Ldual ) then
          write(Uout,9999) ' +++++ Target force calc   ', dual_Iforce
        end if
      end if
      write(Uout,9999) ' +++++ Seed for Random    ', Iseed
      write(Uout,9997) ' +++++ Address of Result   ', trim(dir_result)
      write(Uout,9997) ' +++++ Address of Scratch  ', trim(dir_scr)
      write(Uout,9995) ' +++++ Flag for Restart    ', Lrestart
      write(Uout,9995) ' +++++ Flag for Periodic   ', Lperiodic
      write(Uout,'(a)') ' +++++ lattice (Angstrom)'
      write(Uout,'(3F16.8)') lattice(1,:)
      write(Uout,'(3F16.8)') lattice(2,:)
      write(Uout,'(3F16.8)') lattice(3,:)
      write(Uout,*)

#ifdef _MACE_
      if ( Iforce == 25 ) then
        call get_mace_model_path(mace_model_path)
        write(Uout,'(a)') ' +++++ MACE force field +++++'
        write(Uout,'(a,a)') ' +++++ MACE model      ', trim(mace_model_path)
        write(Uout,'(a,a)') ' +++++ MACE device     ', trim(device)
        write(Uout,*)
      end if
#endif

      if ( Icons > 0 ) then
        write(Uout,'(a)')' +++++ Harmonic constraint +++++'
        write(Uout,9999) ' +++++ Icons                     ', Icons
        write(Uout,9999) ' +++++ cons_atom1                ', cons_atom1
        write(Uout,9999) ' +++++ cons_atom2                ', cons_atom2
        write(Uout,9999) ' +++++ cons_atom3                ', cons_atom3
        write(Uout,9998) ' +++++ cons_val (Angstrom)       ', cons_val
        write(Uout,9994) ' +++++ cons_strength (Ha/bohr^2) ', cons_strength
        write(Uout,*)
      end if

      !if ( Lrestart .eqv. .False. ) then
        write(Uout,*) ' +++++ Atom, Mass, Coordnate(x,y,z) +++++++'
        if ( Langstrom .eqv. .True. ) then
          do Iatom = 1, Natom
            write(Uout,9996) alabel(iatom), PhysMass(iatom)/amu2AU, ur(:,Iatom,1)*AU2Ang
          end do
        else
          do Iatom = 1, Natom
            write(Uout,9996) alabel(iatom), PhysMass(iatom)/amu2AU, ur(:,Iatom,1)
          end do
        end if
        write(Uout,'(a,a)')  '  ',repeat('+',42)
        write(Uout,*)
      !end if
    close(Uout)
  end if

  return
9999 Format(A26,I12)
9998 Format(A26,F12.4)
9997 Format(A26,8X,A)
!9997 Format(A26,A80)
9996 Format(A4,1X,4F12.6)
9995 Format(A26,L12)
9994 Format(A26,ES10.2)
end subroutine Check_Inp
