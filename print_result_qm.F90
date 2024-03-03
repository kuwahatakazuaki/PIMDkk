subroutine print_result_qm
  use Parameters
  use utility, only: program_abort
  implicit none
  integer :: i,j,k, imode, iatom
  integer :: Upre, Udip, Uchar, Uhfc, Ucoor, Ufor, Uene

if(MyRank==0) Then

  if( Lsave_dipole .eqv. .True. ) then
    open(newunit=Udip,file=trim(address)//'/dipole.dat',status='unknown',form='formatted',position='append')
      write(Udip,'("#",I10)') istepsv
      do imode=1,nbead
        write(Udip,8008) dipoler(:,imode)
      end do
    close(Udip)
  endif

  if( Lsave_charge .eqv. .True.) then
    open(newunit=Uchar,file=trim(address)//'/charge.dat',status='unknown',form='formatted',position='append')
      write(Uchar,'("#", I10)') istepsv
      do imode=1,nbead
        write(Uchar,8007) charge(:,imode)
      end do
    close(Uchar)
  endif

  if( Lsave_hfcc .eqv. .True. ) then
    open(newunit=Uhfc,file=trim(address)//'/hfcc.dat',status='unknown',form='formatted',position='append')
    write(Uhfc,'("#", I10)') istepsv
    do imode=1,nbead
      write(Uhfc,8004) hfcc(:,imode)
    end do
    close(Uhfc)
  end if

  if ( Iforce == 8 ) then
    open(newunit=Upre,file=trim(address)//'/pressure.dat',status='unknown',form='formatted',position='append')
      write(Upre,'("#", I10)') istepsv
      do imode=1,nbead
        write(Upre,'(F7.2)') pressure(imode)
      end do
    close(Upre)

    open(newunit=Upre,file=trim(address)//'/PV.dat',status='unknown',form='formatted',position='append')
      write(Upre,*) istepsv, PV
    close(Upre)
  end if

!  open(igetxyz,file=trim(address)//'/cent.xyz',status='unknown',form='formatted',position='append')
!    write(igetxyz,'(I5)') natom
!    write(igetxyz,'(I10)') istepsv
!    do iatom=1,natom
!      write(igetxyz,9999) alabel(iatom),ux(iatom,1)*bohr_inv,uy(iatom,1)*bohr_inv,uz(iatom,1)*bohr_inv
!    end do
!  close(igetxyz)

  open(Ucoor,file=trim(address)//'/coor.xyz',status='unknown',form='formatted',position='append')
    write(Ucoor,'(I5)') natom*nbead
    write(Ucoor,'(I10)') istepsv
    do imode=1,nbead
      do iatom=1,natom
        write(Ucoor,9999) alabel(iatom),r(:,iatom,imode)*AUtoAng
      end do
    end do
  close(Ucoor)

  if (Lsave_force .eqv. .True.) then
    open(newunit=Ufor,file=trim(address)//'/force.dat',status='unknown',form='formatted',position='append')
      write(Ufor,'("#",I10)') istepsv
      do imode=1,nbead
        do iatom=1,natom
          write(Ufor,8005) fr(:,iatom,imode)
        end do
      end do
    close(Ufor)
  end if

!  open(newunit=Uene,file=trim(address)//'/ene.dat',status='unknown',form='formatted',position='append')
!    write(Uene,'("#",I10)') istepsv
!    do imode=1,nbead
!      write(Uene,8006) Eenergy(imode)
!    end do
!  close(Uene)

  potential = sum(Eenergy(:)) * dp_inv
EndIf


!9999 format(a2,1x,E15.9,1x,E15.9,1x,E15.9)
9999 format(a,1x,E15.9,1x,E15.9,1x,E15.9)
9998 format(3E23.15)
9997 format(2E23.15)
9996 format(E23.15)
8004 format(100F12.5)  ! HFCC
8005 format(3F15.10)   ! Force
8006 format(F0.10)     ! Potential
8007 format(100F10.6)  ! Charge
8008 format(4F10.5)    ! Dipole
9995 format(4E23.15)

return
end subroutine print_result_qm
