Subroutine Getfnm

  Use Parameters
  Implicit None
  integer :: i, j, k
  Integer   :: jmode


!   /*  initialize array  */
  Do imode = 1, nbead
     Do iatom = 1, natom
        fux(iatom,imode) = 0.d0
        fuy(iatom,imode) = 0.d0
        fuz(iatom,imode) = 0.d0
     Enddo
  Enddo


!  print *, "f"
!  do i =1, Natom
!    do j = 1, Nbead
!      print *, i,j,fx(i,j), fy(i,j), fz(i,j)
!    end do
!  end do
!  stop "HERE"

!   /*  transformation                      *
!    *  fu(i) = fu(i) + sum_j fx(j)*tnm(j,i)  */
  Do iatom = 1, natom
     Do imode = 1, nbead
        Do jmode = 1, nbead
           fux(iatom,imode) = fux(iatom,imode) + fx(iatom,jmode)*tnm(jmode,imode)
           fuy(iatom,imode) = fuy(iatom,imode) + fy(iatom,jmode)*tnm(jmode,imode)
           fuz(iatom,imode) = fuz(iatom,imode) + fz(iatom,jmode)*tnm(jmode,imode)
        Enddo
     Enddo
  Enddo

!  do i = 1, Nbead
!    do j = 1, Nbead
!      print *, tnm(i,j)
!    end do
!  end do

Return
End Subroutine
