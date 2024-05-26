module mod_model_force
  implicit none
  real(8), parameter :: fcons = 1.0d-1
  real(8), parameter :: eps = 1.0d-10
! Parameters in Atomic unit
  ! === Hydrogen molecule ===
  !real(8), parameter :: r_e   = 1.41014d0
  !real(8), parameter :: De    = 0.1745d0
  !real(8), parameter :: width = 1.0213d0
  ! === Hydrogen molecule ===
  real(8), parameter :: r_e   = 1.41014d0
  real(8), parameter :: De    = 3.0d-3
  real(8), parameter :: width = 2.0d0

contains
  real(8) function harmonic(r1,r0)
    real(8), intent(in) :: r1(3), r0(3)
    real(8) :: rij(3)
    rij(:) = r1(:)-r0(:)
    harmonic = fcons * dot_product(rij,rij)
  end function harmonic

  function Fharmo(r1,r0) result(Fr)
    real(8), intent(in) :: r1(3), r0(3)
    real(8) :: Fr(3)
    real(8) :: rij(3)
    rij(:) = r1(:)-r0(:)
    Fr(:)  = -2.0d0 * fcons * rij(:)
  end function Fharmo

  real(8) function morse(r1,r0)
    real(8), intent(in) :: r1(3), r0(3)
    real(8) :: dis, temp

    dis = norm2(r1(:)-r0(:))
    temp = 1.0d0 - exp(-width*dis)
    morse = De*temp**2
  end function morse

  function Fmorse(r1,r0) result(Fr)
    real(8), intent(in) :: r1(3), r0(3)
    real(8) :: Fr(3)
    real(8) :: dis, power, e(3)
    dis = norm2(r1(:)-r0(:))
    if (dis < eps) then
      Fr(:) = 0.0d0
    else
      power = exp(-width*dis)
      e(:)  = (r1(:)-r0(:))/dis
      Fr(:) = -2.0d0*De*width*power*(1.0d0-power)*e(:)
    end if
  end function Fmorse


end module mod_model_force


