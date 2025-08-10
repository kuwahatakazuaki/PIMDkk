function get_kinetic_ene() result(kine)
  use Parameters
  use utility, only: norm_seq
  implicit none
  real(8) :: kine
  integer :: i, j
  kine = 0.0d0
  do i = 1, Natom
    do j = 1, Nbead
      kine = kine + fictmass(i,j) * norm_seq(vur(:,i,j))
    end do
  end do
  kine = kine * 0.5d0
end function get_kinetic_ene

