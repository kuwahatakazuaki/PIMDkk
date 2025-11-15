subroutine getforce_ref_nor
!subroutine Getforce_Ref
  use Parameters
  use utility, only: program_abort
  implicit none
  integer :: i, j

  fur_ref(:,:,1) = 0.0d0
  do j = 2, Nbead
    do i = 1, Natom
      fur_ref(:,i,j) = -dnmmass(i,j)*omega_p2*ur(:,i,j)
    end do
  end do

return
end subroutine getforce_ref_nor
!end subroutine Getforce_Ref
