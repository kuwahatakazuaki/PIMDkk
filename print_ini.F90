subroutine print_ini
  use Parameters
  use utility, only: program_abort
  implicit none
  integer :: i, Uham, Uchar, Uhfc, Uout, Udip, Ucon, Upre

if (MyRank == 0) then

! +++ Start Nbead == 1 +++
  if (Nbead == 1) then

    open(newunit=Uham,file=trim(dir_result)//'/ham.dat',status='replace',form='formatted',position='append')
      write(Uham,'("# ",a)') repeat('*',87)
      write(Uham,'(a)') '#  1Step 2 Hamiltonian   3 Temperature   4 Potential     5 DKinetic      6 EBath_Cent '
      write(Uham,'("# ",a)') repeat('*',87)
    close(Uham)

    open(newunit=Uout,file=Fout,status='old',position='append')
      write(Uout,'(" ",a)') repeat('*',95)
      write(Uout,'(a)') '   Step  Hamiltonian  Temperature  Potential    DKinetic     EBath_Cent    Date & Time'
      write(Uout,'(" ",a)') repeat('*',95)
    close(Uout)

! +++ End Nbead == 1 +++
  else if (Nbead > 1) then
! +++ Start Nbead > 1 +++

    open(newunit=Uham,file=trim(dir_result)//'/ham.dat',status='replace',form='formatted',position='append')
      write(Uham,'("# ",a)')  repeat('*',132)
      write(Uham,'(a)') &
  &'# 1Step  2 Hamiltonian   3 Temperature   4 Potential     5 DKinetic      6 QKinetic      7 EBath         &
  8 EBath_Cent    9 E_Virial   '
      write(Uham,'("#  ",a)')  repeat('*',132)
    close(Uham)

    open(newunit=Uout,file=Fout,status='old',position='append')
      write(Uout,'(" ",a)')  repeat('*',121)
      write(Uout,'(a)') &
        '   Step  Hamiltonian  Temperature  Potential    DKinetic     QKinetic     EBath        EBath_Cent     Date & Time'
      write(Uout,'(" ",a)')  repeat('*',121)
    close(Uout)

  end if
! +++ End Nbead > 1 +++

  if( Lsave_charge .eqv. .True.) then
    open(newunit=Uchar,file=trim(dir_result)//'/charge.dat',status='replace',form='formatted',position='append')
      write(Uchar,'(a)',advance='no') "#"
      do i = 1, Natom
        write(Uchar,'(a7,I2,1x)',advance='no') alabel(i), i
      end do
    close(Uchar)
  endif

  if( Lsave_hfcc .eqv. .True. ) then
    open(newunit=Uhfc,file=trim(dir_result)//'/hfcc.dat',status='replace',form='formatted',position='append')
      write(Uhfc,'(a)',advance='no') "#"
      do i = 1, Natom
        write(Uhfc,'(a9,I2,1x)',advance='no') alabel(i), i
      end do
    close(Uhfc)
  end if

  if( Lsave_dipole .eqv. .True. ) then
    open(newunit=Udip,file=trim(dir_result)//'/dipole.dat',status='replace',form='formatted',position='append')
      write(Udip,'("# Nbead = ",I10)') Nbead
    close(Udip)
  endif

  if ( Icons > 0 ) then
    open(newunit=Ucon,file=trim(dir_result)//'/constraint.dat',status='replace',form='formatted',position='append')
      write(Ucon,'(a)') '#  1step 2time_fs       3cons_cv       4dVdcons_cv 5potential_cons'
    close(Ucon)
  end if

  if ( Lperiodic ) then
    open(newunit=Upre,file=trim(dir_result)//'/pressure.dat',status='replace',form='formatted',position='append')
      write(Upre,'(a)') '# 1step 2time_fs 3P_cv_bar 4P_prim_bar 5P_inst_bar 6volume_A3 7a_A 8b_A 9c_A 10W_pot_au'
    close(Upre)
  end if

end if

return
end subroutine print_ini
