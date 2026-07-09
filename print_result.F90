subroutine print_result
  use Parameters
  use utility, only: program_abort
  implicit none
  integer :: i,j,k, imode, iatom
  integer :: Upre, Udip, Uchar, Uhfc, Ucoor, Ufor, Uene, Ucon
  character(len=5) :: pbc_xyz
  real(8) :: cell_len(3)

if (MyRank==0) then

  if ( Lsave_dipole .eqv. .True. ) then
    open(newunit=Udip,file=trim(dir_result)//'/dipole.dat',status='unknown',form='formatted',position='append')
      write(Udip,'("#",I10)') istepsv
      do imode=1,nbead
        write(Udip,8008) dipoler(:,imode)
      end do
    close(Udip)
  end if

  if ( Lsave_charge .eqv. .True.) then
    open(newunit=Uchar,file=trim(dir_result)//'/charge.dat',status='unknown',form='formatted',position='append')
      write(Uchar,'("#", I10)') istepsv
      do imode=1,nbead
        write(Uchar,8007) charge(:,imode)
      end do
    close(Uchar)
  end if

  if ( Lsave_hfcc .eqv. .True. ) then
    open(newunit=Uhfc,file=trim(dir_result)//'/hfcc.dat',status='unknown',form='formatted',position='append')
    write(Uhfc,'("#", I10)') istepsv
    do imode=1,nbead
      write(Uhfc,8004) hfcc(:,imode)
    end do
    close(Uhfc)
  end if

  if ( Lperiodic ) then
    do i = 1, 3
      cell_len(i) = sqrt(sum(lattice(i,:)**2))
    end do
    open(newunit=Upre,file=trim(dir_result)//'/pressure.dat',status='unknown',form='formatted',position='append')
      write(Upre,8010) istepsv, dble(istepsv) * dt / fs2AU, &
                       press_cv / bar2AU, press_prim / bar2AU, press_inst / bar2AU, &
                       vol * AU2Ang**3, cell_len(:), W_pot
    close(Upre)
  end if

  open(Ucoor,file=trim(dir_result)//'/coor.xyz',status='unknown',form='formatted',position='append')
    if (Lperiodic .eqv. .True.) then
      pbc_xyz = 'T T T'
    else
      pbc_xyz = 'F F F'
    end if
    write(Ucoor,'(I5)') natom*nbead
    if (Lperiodic .eqv. .True.) then
      write(Ucoor,8011) lattice(1,:), lattice(2,:), lattice(3,:), pbc_xyz, istepsv
    else
      write(Ucoor,'(a,a,a,i0)') 'Properties=species:S:1:pos:R:3 pbc="', pbc_xyz, '" step=', istepsv
    end if
    do imode=1,nbead
      do iatom=1,natom
        write(Ucoor,9999) alabel(iatom),r(:,iatom,imode)*AU2Ang
      end do
    end do
  close(Ucoor)

  if (Lsave_force .eqv. .True.) then
    open(newunit=Ufor,file=trim(dir_result)//'/force.dat',status='unknown',form='formatted',position='append')
      write(Ufor,'("#",I10)') istepsv
      do imode=1,Nbead
        do iatom=1,natom
          write(Ufor,8005) fr(:,iatom,imode)
        end do
      end do
    close(Ufor)
  end if

  if (Lsave_energy.eqv. .True.) then
    open(newunit=Uene,file=trim(dir_result)//'/ene.dat',status='unknown',form='formatted',position='append')
      write(Uene,'("#",I10)') istepsv
      do imode=1,Nbead
        write(Uene,8006) pot_bead(imode)
      end do
    close(Uene)
  end if

  if ( Icons > 0 ) then
    open(newunit=Ucon,file=trim(dir_result)//'/constraint.dat',status='unknown',form='formatted',position='append')
      write(Ucon,8009) istepsv, dble(istepsv) * dt / fs2AU, cons_cv_ave * AU2Ang, &
                       dVdcons_cv, potential_cons
    close(Ucon)
  end if

end if


!9999 format(a2,1x,E15.9,1x,E15.9,1x,E15.9)
9999 format(a,E15.7,E15.7,E15.7)
!9999 format(a,1x,E15.9,1x,E15.9,1x,E15.9)
9998 format(3E23.15)
9997 format(2E23.15)
9996 format(E23.15)
8004 format(100F12.5)  ! HFCC
8005 format(3E15.7)   ! Force
! 8005 format(3F15.10)   ! Force
8006 format(F0.10)     ! Potential
8007 format(100F10.6)  ! Charge
8008 format(4F10.5)    ! Dipole
8009 format(I7,4E15.7)
8010 format(I7,9E17.9)
8011 format('Lattice="',9(ES23.15,1X),'" Properties=species:S:1:pos:R:3 pbc="',a,'" step=',i0)
9995 format(4E23.15)

return
end subroutine print_result
