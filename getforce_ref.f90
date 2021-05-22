subroutine getforce_ref
  use global_variable
  implicit none
  integer :: j, xyz

  fu_ref(:,:,1) = 0.0d0
  do j = 2, Nbead
    do xyz = 1, 3
      fu_ref(xyz,:,j) = -dnmmass(:,j) * omega_p2 * u(xyz,:,j)
    end do
  end do
end subroutine getforce_ref

!Subroutine Getforce_Ref
!
!Use Parameters
!Implicit None
!
!  do iatom = 1, natom
!     fux_ref(iatom,1) = 0.d0
!     fuy_ref(iatom,1) = 0.d0
!     fuz_ref(iatom,1) = 0.d0
!  enddo
!
!  do imode = 1, nbead
!     do iatom = 1, natom
!        fux_ref(iatom,imode) = -dnmmass(iatom,imode)*omega_p2*ux(iatom,imode)
!        fuy_ref(iatom,imode) = -dnmmass(iatom,imode)*omega_p2*uy(iatom,imode)
!        fuz_ref(iatom,imode) = -dnmmass(iatom,imode)*omega_p2*uz(iatom,imode)
!     enddo
!  enddo
!
!
!Return
!End Subroutine
