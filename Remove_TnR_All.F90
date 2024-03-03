Subroutine Remove_TnR_All
  Use Parameters
  use utility, only: outer_product, program_abort, calc_determinant33, calc_inv_mat33
  Implicit None
  real(8), parameter :: delta = 0.00001d0
  Double Precision :: detmat! ,gradx,grady,gradz
  Double Precision, Dimension(1:3,1:3) :: dinvmat
  integer :: j, i
  real(8) :: sumvr(3), comr(3), gradr(3), totmas
  real(8) :: inertia(3,3), ang_vel(3), moment(3)
!
!   Remove Translation and Rotation Velocity of 
!   COM for Centroids and Noncentroids
!
  do j=1, Nbead
    comr(:) = 0.0d0
    sumvr(:) = 0.0d0
    ang_vel(:) = 0.0d0
!
! Calculate Center of Mass and Translational Velocity of COM
!
    do i = 1, natom
      sumvr(:) = sumvr(:)  + vur(:,i,j)*fictmass(i,j)
      comr(:)  = comr(:)   + ur(:,i,j)*fictmass(i,j)
      !totmas = totmas + fictmass(i,j)
    enddo
    totmas = sum(fictmass(:,j))
    sumvr(:) = sumvr(:) / totmas
    comr(:)  = comr(:)  / totmas

!
! Calculate Rotational Velocity of COM
! Construct the Momenta 
!
    inertia(:,:) = 0.0d0
    do i = 1, natom
      !cixx = cixx + fictmass(i,j) * ((uy(i,j)-comy)**2+(uz(i,j)-comz)**2)
      !ciyy = ciyy + fictmass(i,j) * ((ux(i,j)-comx)**2+(uz(i,j)-comz)**2)
      !cizz = cizz + fictmass(i,j) * ((uy(i,j)-comy)**2+(ux(i,j)-comx)**2)
      !cixy = cixy - fictmass(i,j) *  (ux(i,j)-comx) *  (uy(i,j)-comy)
      !ciyz = ciyz - fictmass(i,j) *  (uz(i,j)-comz) *  (uy(i,j)-comy)
      !cizx = cizx - fictmass(i,j) *  (ux(i,j)-comx) *  (uz(i,j)-comz)
      inertia(:,:) = inertia(:,:) + calc_inertia(fictmass(i,j),ur(:,i,j)-comr(:))
    enddo
!
! Obtain the Inverse Matrix of Momenta
!
    detmat = calc_determinant33(inertia)
!if ( detmat <= 0.00001d0 ) goto 100
    dinvmat(:,:) = 0.0d0
    if ( detmat <= delta ) then
      do i = 1, 3
        if ( abs(inertia(i,i)) > delta ) dinvmat(i,i) = 1/inertia(i,i)
      end do
    else
      dinvmat(:,:) = calc_inv_mat33(inertia)
!      do i = 1, natom
!        moment(:) = fictmass(i,j) * outer_product(ur(:,i,j)-comr(:), vur(:,i,j))
!        ang_vel(:) = ang_vel(:) + matmul(dinvmat,moment)
!      enddo
    end if

    moment(:) = 0.0d0
    do i = 1, natom
      moment(:) = moment(:) + fictmass(i,j) * outer_product(ur(:,i,j)-comr(:), vur(:,i,j))
    enddo
    ang_vel(:) = matmul(dinvmat,moment)

!100   continue
!
! Subtract Translation and Rotation Velocities of COM from 
! Velocities of Atoms
!
    do i = 1, natom
      vur(:,i,j) = vur(:,i,j) - sumvr(:)
      vur(:,i,j) &
        = vur(:,i,j) - outer_product(ang_vel,ur(:,i,j)-comr(:))
    enddo

!moment(:) = 0.0d0
!do i = 1, natom
!  moment(:) = moment(:) + fictmass(i,j) * outer_product(ur(:,i,j)-comr(:), vur(:,i,j))
!enddo
!print *, moment(:)

  enddo
  Return
contains

  function calc_inertia(m,r) result(L)
    real(8), intent(in) :: m,r(3)
    real(8) :: L(3,3)
    real(8) :: x, y, z
    x = r(1)
    y = r(2)
    z = r(3)
    L(1,1) = m * (y**2+z**2)
    L(2,2) = m * (z**2+x**2)
    L(3,3) = m * (x**2+y**2)
    L(1,2) = (-1) * m*x*y
    L(1,3) = (-1) * m*z*x
    L(2,3) = (-1) * m*y*z
    L(2,1) = L(1,2)
    L(3,1) = L(1,3)
    L(3,2) = L(2,3)
  end function calc_inertia

End Subroutine Remove_TnR_All
