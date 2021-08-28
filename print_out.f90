module print_out


contains
<<<<<<< HEAD

subroutine print_result
=======
subroutine print_gaussian
>>>>>>> origin/main
  use global_variable
  implicit none
  integer :: i, j
  integer :: Ucoor, Ucent, Udipo, Uchar, Ufor, Uene

  open(newunit=Ucent, file=path_result//'/cent.xyz', position='append')
    write(Ucent, '(I5)') Natom
    write(Ucent,'(I10)') Istep
    do i = 1, Natom
      write(Ucent,9999) alabel(i), u(:,i,1) * AUtoAng
    end do
  close(Ucent)

  open(newunit=Ucoor, file=path_result//'/coor.xyz', position='append')
    write(Ucoor, '(I5)') Natom
    write(Ucoor,'(I10)') Istep
    do j = 1, Nbead
      do i = 1, Natom
        write(Ucoor,9999) alabel(i), r(:,i,j) * AUtoAng
      end do
    end do
  close(Ucoor)

  open(newunit=Uene, file=path_result//'/ene.dat', position='append')
    write(Uene,'(I10)') Istep
    do j = 1, Nbead
      write(Uene,8006) energy(j)
    end do
  close(Uene)

  if ( Ldipole .eqv. .True. ) then
    open(newunit=Udipo, file=path_result//'/dipole.dat', position='append')
      write(Udipo,'(I10)') Istep
      do j = 1, Nbead
        write(Udipo, 8008) dipole(:,j)
      end do
    close(Udipo)
  end if

  if ( Lcharge .eqv. .True. ) then
    open(newunit=Uchar, file=path_result//'/charge.dat', position='append')
      write(Uchar,'(I10)') Istep
      do j = 1, Nbead
        write(Uchar, 8007) charge(:,j)
      end do
    close(Uchar)
  end if

  if ( Lforce .eqv. .True. ) then
    open(newunit=Ufor, file=path_result//'/force.dat', position='append')
      write(Ufor,'(I10)') Istep
      do j = 1, Nbead
        do i = 1, Natom
          write(Ufor,8005) f(:,i,j)
        end do
      end do
    close(Ufor)
  end if

!  potential=0.D0
!  DO imode=1,nbead
!     potential=potential+Eenergy(imode)
!  ENDDO
!  potential=potential*dp_inv

!9998 format(3E23.15)
!9997 format(2E23.15)
!9996 format(E23.15)
!9995 format(4E23.15)

return
9999 format(a2,1x,E15.9,1x,E15.9,1x,E15.9)
8008 format(4F10.5)    ! Dipole
8007 format(100F10.6)  ! Charge
8006 format(F0.10)     ! Potential
8005 format(3F15.10)   ! Force
<<<<<<< HEAD
end subroutine print_result
=======
end subroutine print_gaussian
>>>>>>> origin/main

end module print_out

!  open(igetxyz,file=trim(address)//'/cent.xyz',status='unknown',form='formatted',position='append')
!    write(igetxyz,'(I5)') natom
!    write(igetxyz,'(I10)') istepsv
!    do iatom=1,natom
!      write(igetxyz,9999) alabel(iatom),ux(iatom,1)*bohr_inv,uy(iatom,1)*bohr_inv,uz(iatom,1)*bohr_inv
!    end do
!  close(igetxyz)

!  open(igetx,file=trim(address)//'/coor.xyz',status='unknown',form='formatted',position='append')
!    write(igetx,'(I5)') natom*nbead
!    write(igetx,'(I10)') istepsv
!    do imode=1,nbead
!      do iatom=1,natom
!        write(igetx,9999) alabel(iatom),x(iatom,imode)*bohr_inv,y(iatom,imode)*bohr_inv,z(iatom,imode)*bohr_inv
!      end do
!    end do
!  close(igetx)

!  if(nodipole==0) then
!    Open(igetd,file=trim(address)//'/dipole.dat',status='unknown',form='formatted',position='append')
!      write(igetd,'(I10)') istepsv
!      do imode=1,nbead
!        write(igetd,8008) dipolex(imode),dipoley(imode),dipolez(imode),dipole(imode)
!      end do
!    close(igetd)
!  endif
!  if(nocharge==0) then
!    open(igetc,file=trim(address)//'/charge.dat',status='unknown',form='formatted',position='append')
!      if(nohfcc==0) then
!        DO imode=1,nbead
!          DO iatom=1,natom
!             write(igetc,9997) charge(iatom,imode), hfcc(iatom,imode)
!          ENDDO
!        ENDDO
!      else
!        write(igetc,'(I10)') istepsv
!        do imode=1,nbead
!          write(igetc,8007) charge(:,imode)
!        end do
!      endif
!    close(igetc)
!  endif
!
!  if (Save_force .eqv. .True.) then
!    open(igetf,file=trim(address)//'/force.dat',status='unknown',form='formatted',position='append')
!      do imode=1,nbead
!        do iatom=1,natom
!          write(igetf,9998) fx(iatom,imode),fy(iatom,imode),fz(iatom,imode)
!        end do
!      end do
!    close(igetf)
!  end if
!
!  open(igete,file=trim(address)//'/ene.dat',status='unknown',form='formatted',position='append')
!    write(igete,'(I10)') istepsv
!    do imode=1,nbead
!      write(igete,8006) Eenergy(imode)
!    end do
!  close(igete)
