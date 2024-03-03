
! +++ gaussian distribution of normal modes +++
! +++ corresponding to free particle        +++

Subroutine NM_Position
  use Parameters
  use utility, only: program_abort
  implicit none
  Double Precision :: gasd, usigma
  integer :: i, j


  If( Lrandom_coor .eqv. .True. ) Then
    do j = 2, nbead
      do i = 1, natom
        usigma = dsqrt(1.d0/beta/omega_p2/dnmmass(i,j))
        call gasdev(gasd); ur(1,i,j) = usigma*gasd
        call gasdev(gasd); ur(2,i,j) = usigma*gasd
        call gasdev(gasd); ur(3,i,j) = usigma*gasd
      enddo
    enddo
  Else
    ur(:,:,2:Nbead) = 0.d0
  EndIf

  Return
End Subroutine
