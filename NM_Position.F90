
! +++ gaussian distribution of normal modes +++
! +++ corresponding to free particle        +++

subroutine NM_Position
  use Parameters
  use utility, only: program_abort, gasdev
  implicit none
  real(8) :: gasd, usigma
  integer :: i, j, Idim


  if ( Lrandom_coor .eqv. .True. ) then
    do j = 2, Nbead
      do i = 1, Natom
        usigma = dsqrt(1.d0/beta/omega_p2/dnmmass(i,j))
        !ur(1,i,j) = usigma*gasdev()
        !ur(2,i,j) = usigma*gasdev()
        !ur(3,i,j) = usigma*gasdev()
        do Idim = 1, Ndim
          ur(Idim,i,j) = usigma*gasdev()
        end do
      enddo
    enddo
  else
    ur(:,:,2:Nbead) = 0.d0
  end if

  return
end subroutine
