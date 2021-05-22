Subroutine Write_MPI_tk

  Use Parameters
  Use Parameter_tk
  Use MPI
  Implicit None
!  Logical :: lopen
  Integer :: i,j,k

! #if 0
!     If(MyRank==kproc) then
!       if(nodipole==0) then
!         Do imode=1,k
!          Call MPI_Wait(ireqdx(imode),mstatus,IERR)
!          Call MPI_Wait(ireqdy(imode),mstatus,IERR)
!          Call MPI_Wait(ireqdz(imode),mstatus,IERR)
!         Enddo
!       endif
!       if(nocharge==0) then
!         Do imode=1,k
!          Call MPI_Wait(ireqc(imode),mstatus,IERR)
!         Enddo
!       endif
!     Else
!        if(nodipole==0) then
!         Call MPI_Wait(ireqdx(1),mstatus,IERR)
!         Call MPI_Wait(ireqdy(1),mstatus,IERR)
!         Call MPI_Wait(ireqdz(1),mstatus,IERR)
!        endif
!        if(nocharge==0) then
!         Call MPI_Wait(ireqc(1),mstatus,IERR)
!        endif
!     Endif
!     If(MyRank==0) then
!       Do imode=1,k
!        Call MPI_Wait(ireqe(imode),mstatus,IERR)
!        Call MPI_Wait(ireqfx(imode),mstatus,IERR)
!        Call MPI_Wait(ireqfy(imode),mstatus,IERR)
!        Call MPI_Wait(ireqfz(imode),mstatus,IERR)
!       Enddo
!     Else
!        Call MPI_Wait(ireqe(1),mstatus,IERR)
!        Call MPI_Wait(ireqfx(1),mstatus,IERR)
!        Call MPI_Wait(ireqfy(1),mstatus,IERR)
!        Call MPI_Wait(ireqfz(1),mstatus,IERR)
!     Endif
! #endif

! kproc = 0
  If(MyRank==kproc) Then
    if(nodipole==0) then
!      Do imode=1,k
!        Call MPI_Wait(ireqdx(imode),mstatus,IERR)
!        Call MPI_Wait(ireqdy(imode),mstatus,IERR)
!        Call MPI_Wait(ireqdz(imode),mstatus,IERR)
!      Enddo
      write(igetd,'(I10)') istepsv
      do imode=1,nbead
!        write(igetd,9998) dipolex(imode),dipoley(imode),dipolez(imode)
        write(igetd,8008) dipolex(imode),dipoley(imode),dipolez(imode)
      end do
    endif
    if(nocharge==0) then
!      Do imode=1,k
!         Call MPI_Wait(ireqc(imode),mstatus,IERR)
!      Enddo
      if(nohfcc==0) then
        DO imode=1,nbead
          DO iatom=1,natom
             write(igetc,9997) charge(iatom,imode), hfcc(iatom,imode)
          ENDDO
        ENDDO
      else
        write(igetc,'(I10)') istepsv
        do imode=1,nbead
          write(igetc,8007) charge(:,imode)
!           DO iatom=1,natom
!              write(igetc,9996) charge(iatom,imode)
!           ENDDO
        end do
      endif
    endif
!  Else
!     if(nodipole==0) then
!      Call MPI_Wait(ireqdx(1),mstatus,IERR)
!      Call MPI_Wait(ireqdy(1),mstatus,IERR)
!      Call MPI_Wait(ireqdz(1),mstatus,IERR)
!     endif
!     if(nocharge==0) then
!      Call MPI_Wait(ireqc(1),mstatus,IERR)
!     endif
  Endif

  If(MyRank==0) Then
     write(igetxyz,'(I5)') natom
     write(igetxyz,'(I10)') istepsv
     DO iatom=1,natom
        write(igetxyz,9999) alabel(iatom),ux(iatom,1)*bohr_inv,uy(iatom,1)*bohr_inv,uz(iatom,1)*bohr_inv
     ENDDO
     write(igetx,'(I5)') natom*nbead
     write(igetx,'(I10)') istepsv
     DO imode=1,nbead
        DO iatom=1,natom
           write(igetx,9999) alabel(iatom),x(iatom,imode)*bohr_inv,y(iatom,imode)*bohr_inv,z(iatom,imode)*bohr_inv
        ENDDO
     ENDDO
!     Do imode=0,nprocs-2
!       Call MPI_Wait(ireqfx(imode),mstatus,IERR)
!       Call MPI_Wait(ireqfy(imode),mstatus,IERR)
!       Call MPI_Wait(ireqfz(imode),mstatus,IERR)
!     Enddo
     if (Save_force .eqv. .True.) then
       DO imode=1,nbead
          DO iatom=1,natom
             write(igetf,9998) fx(iatom,imode),fy(iatom,imode),fz(iatom,imode)
          ENDDO
       ENDDO
     end if
!     Do imode=0,nprocs-2
!       Call MPI_Wait(ireqe(imode),mstatus,IERR)
!     Enddo
     DO imode=1,nbead
        write(igete,9996) Eenergy(imode)
     ENDDO
     potential=0.D0
     DO imode=1,nbead
        potential=potential+Eenergy(imode)
     ENDDO
     potential=potential*dp_inv
!    DO imode=1,nbead
!       write(igethl,9997) homo(imode),lumo(imode)
!    ENDDO
!  Else
!     Call MPI_Wait(ireqe(1),mstatus,IERR)
!     Call MPI_Wait(ireqfx(1),mstatus,IERR)
!     Call MPI_Wait(ireqfy(1),mstatus,IERR)
!     Call MPI_Wait(ireqfz(1),mstatus,IERR)
  EndIf

9999 format(a2,1x,E15.9,1x,E15.9,1x,E15.9) 
9998 format(3E23.15)
9997 format(2E23.15)
9996 format(E23.15)
8007 format(100F10.6)  ! Charge
8008 format(3F10.5)    ! Dipole
9995 format(4E23.15)

Return
End Subroutine
