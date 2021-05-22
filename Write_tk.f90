  Subroutine Write_tk

    Use Parameters
    Use Parameter_tk
!    Use MPI
    Implicit None

    Logical                            :: lopen
    Integer :: i,j,k


#if 1
!    If(MyRank==kproc) Then
          DO imode=1,nbead
            write(igetd,9998) dipolex(imode),dipoley(imode),dipolez(imode)
          ENDDO
          DO imode=1,nbead
             DO iatom=1,natom
                write(igetc,9996) charge(iatom,imode)
             ENDDO
          ENDDO
!    Endif
#endif

!    If(MyRank==0) Then
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
       DO imode=1,nbead
          DO iatom=1,natom
             write(igetf,9998) fx(iatom,imode),fy(iatom,imode),fz(iatom,imode)
          ENDDO
       ENDDO
       DO imode=1,nbead
          write(igete,9996) Eenergy(imode)
       ENDDO
       potential=0.D0
       DO imode=1,nbead
          potential=potential+Eenergy(imode)
       ENDDO
       potential=potential*dp_inv
!    EndIf


 9999 format(a2,1x,E15.9,1x,E15.9,1x,E15.9) 
 9998 format(3E23.15) 
 9997 format(2E23.15) 
 9996 format(E23.15) 
 9995 format(4E23.15) 

    Return
  End Subroutine
