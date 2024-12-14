subroutine Check_Inp
  use Parameters
  implicit none
  integer :: Iatom, Uout

  if ( MyRank == 0 ) then
    open(newunit=Uout,file=Fout,status='old',position='append')
      write(Uout,'(a)')' +++++ Input Check +++++   '
      write(Uout,9999) ' +++++ Isimulation         ', Isimulation
      write(Uout,9997) ' +++++ Simulation type     ', name_simulation
      write(Uout,9999) ' +++++ Number of Atoms     ', Natom
      write(Uout,9999) ' +++++ Number of Beads     ', Nbead
      write(Uout,9999) ' +++++ Number of Steps     ', Nstep
#ifdef _mpi_
      write(Uout,9999) ' +++++ Number of Nproc     ', Nproc
#endif
      write(Uout,9998) ' +++++ Given Temperature   ', Temperature
      write(Uout,9998) ' +++++ Given time step     ', dt / fs2AU
      write(Uout,9998) ' +++++ Output step (fs)    ', dt / fs2AU * out_step
      !write(Uout,9998) ' +++++ Adiabaticity param. ', gamma1
      write(Uout,9999) ' +++++ Method of Centr NHC ', ncent
      write(Uout,9999) ' +++++ Length of Centr NHC ', nnhc
      write(Uout,9995) ' +++++ Flag for Restart    ', Lrestart
      write(Uout,9999) ' +++++ Flag for Force Calc ', Iforce
      write(Uout,9999) ' +++++ Seed for Random No.1', Iseeds(1)
      write(Uout,9999) ' +++++ Seed for Random No.2', Iseeds(2)
      write(Uout,9999) ' +++++ Seed for Random No.3', Iseeds(3)
      write(Uout,9999) ' +++++ Seed for Random No.4', Iseeds(4)
      write(Uout,9997) ' +++++ Address of Result   ', trim(dir_result)
      write(Uout,9997) ' +++++ Address of Scratch  ', trim(dir_scr)
      write(Uout,*)

      if ( Iumb > 0 ) then
        write(Uout,'(a)')' +++++ Umbrella sampling +++++   '
        write(Uout,9999) ' +++++ Iumb                      ', Iumb
        write(Uout,9994) ' +++++ umb_cons                  ', umb_cons
        write(Uout,9999) ' +++++ umb_atom1                 ', umb_atom1
        write(Uout,9999) ' +++++ umb_atom2                 ', umb_atom2
        write(Uout,9999) ' +++++ umb_atom3                 ', umb_atom3
        write(Uout,*)
      end if

      !if ( Lrestart .eqv. .False. ) then
        write(Uout,*) ' +++++ Atomic Label, Mass, and Coords +++++'
        if ( Langstrom .eqv. .True. ) then
          do Iatom = 1, Natom
            write(Uout,9996) alabel(iatom), PhysMass(iatom)/factmass, ur(:,Iatom,1)*AUtoAng
          end do
        else
          do Iatom = 1, Natom
            write(Uout,9996) alabel(iatom), PhysMass(iatom)/factmass, ur(:,Iatom,1)
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
9997 Format(A26,6X,A)
!9997 Format(A26,A80)
9996 Format(A4,1X,4F12.6)
9995 Format(A26,L12)
9994 Format(A26,ES12.4)
end subroutine Check_Inp


