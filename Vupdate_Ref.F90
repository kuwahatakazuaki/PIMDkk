subroutine Vupdate_Ref
  use Parameters
  use utility, only: program_abort
  implicit none
  integer:: j, i

  do j = 2, Nbead
    do i = 1, Natom
      vur(:,i,j) = vur(:,i,j) + 0.5d0*dt_ref*fur_ref(:,i,j)/fictmass(i,j)
    end do
  end do
end subroutine
