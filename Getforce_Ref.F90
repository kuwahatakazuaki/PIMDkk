Subroutine Getforce_Ref
  use Parameters
  use utility, only: program_abort
  implicit none
  integer :: i, j

  fur_ref(:,:,1) = 0.0d0
  do j = 2, nbead
    do i = 1, natom
      fur_ref(:,i,j) = -dnmmass(i,j)*omega_p2*ur(:,i,j)
    enddo
  enddo

return
End Subroutine
