
! +++ gaussian distribution of normal modes +++
! +++ corresponding to free particle        +++

subroutine NM_Position
  use Parameters
  use utility, only: program_abort, fill_gaussian_random
  implicit none
  real(8) :: usigma
  integer :: i, j


  if ( Lrandom_coor .eqv. .True. ) then
    do j = 2, Nbead
      do i = 1, Natom
        usigma = dsqrt(1.d0/beta/omega_p2/dnmmass(i,j))
        call fill_gaussian_random(ur(:,i,j))
        ur(:,i,j) = usigma * ur(:,i,j)
      enddo
    enddo
  else
    ur(:,:,2:Nbead) = 0.d0
  end if

  return
end subroutine
