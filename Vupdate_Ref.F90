Subroutine Vupdate_Ref
  use Parameters
  use utility, only: program_abort
  implicit none
  integer:: j, i

  do j = 2, nbead
    do i = 1, natom
      vur(:,i,j) = vur(:,i,j) + 0.5d0*dt_ref*fur_ref(:,i,j)/fictmass(i,j)
    enddo
  enddo

End Subroutine
