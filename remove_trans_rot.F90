subroutine remove_trans_rot_cent
  use Parameters
  use utility, only: crossr_product, program_abort, calc_determinant33, calc_inv_mat33, calc_inertia
  implicit none
  real(8), parameter :: delta = 0.00001d0
  real(8) :: dinvmat(3,3), detmat
  integer :: i
  real(8) :: sumvr(3), comr(3), gradr(3), totmas
  real(8) :: inertia(3,3), ang_vel(3), moment(3)
!
! Remove Translation and Rotation Velocity of
! COM for the centroid mode only (j = 1)
!
  comr(:) = 0.0d0
  sumvr(:) = 0.0d0
  ang_vel(:) = 0.0d0
!
! Calculate Center of Mass and Translational Velocity of COM
!
  do i = 1, Natom
    sumvr(:) = sumvr(:) + vur(:,i,1)*fictmass(i,1)
    comr(:)  = comr(:)  + ur(:,i,1)*fictmass(i,1)
  end do
  totmas = sum(fictmass(:,1))
  sumvr(:) = sumvr(:) / totmas
  comr(:)  = comr(:)  / totmas

!
! Calculate Rotational Velocity of COM
! Construct the Momenta of Inertia
!
  inertia(:,:) = 0.0d0
  do i = 1, Natom
    inertia(:,:) = inertia(:,:) + calc_inertia(fictmass(i,1),ur(:,i,1)-comr(:))
  end do
!
! Obtain the Inverse Matrix of Momenta
!
  detmat = calc_determinant33(inertia)
  dinvmat(:,:) = 0.0d0

  if ( detmat <= delta ) then
    do i = 1, 3
      if ( abs(inertia(i,i)) > delta ) dinvmat(i,i) = 1/inertia(i,i)
    end do
  else
    dinvmat(:,:) = calc_inv_mat33(inertia)
  end if

  moment(:) = 0.0d0
  do i = 1, Natom
    moment(:) = moment(:) + fictmass(i,1) * crossr_product(ur(:,i,1)-comr(:), vur(:,i,1))
  end do
  ang_vel(:) = matmul(dinvmat,moment)

!
! Subtract Translation and Rotation Velocities of COM from Velocities of Atoms
!
  do i = 1, Natom
    vur(:,i,1) = vur(:,i,1) - sumvr(:)
    vur(:,i,1) &
      = vur(:,i,1) - crossr_product(ang_vel,ur(:,i,1)-comr(:))
  end do
  return

end subroutine remove_trans_rot_cent


!subroutine Remove_TnR_All
subroutine remove_trans_rot_beads
  use Parameters
  use utility, only: crossr_product, program_abort, calc_determinant33, calc_inv_mat33, calc_inertia
  implicit none
  real(8), parameter :: delta = 0.00001d0
  real(8) :: dinvmat(3,3), detmat
  integer :: i, j
  real(8) :: sumvr(3), comr(3), gradr(3), totmas
  real(8) :: inertia(3,3), ang_vel(3), moment(3)
!
! Remove Translation and Rotation Velocity of 
! COM for Centroids and Noncentroids
!
  do j=1, Nbead
    comr(:) = 0.0d0
    sumvr(:) = 0.0d0
    ang_vel(:) = 0.0d0
!
! Calculate Center of Mass and Translational Velocity of COM
!
    do i = 1, Natom
      sumvr(:) = sumvr(:) + vur(:,i,j)*fictmass(i,j)
      comr(:)  = comr(:)  + ur(:,i,j)*fictmass(i,j)
    end do
    totmas = sum(fictmass(:,j))
    sumvr(:) = sumvr(:) / totmas
    comr(:)  = comr(:)  / totmas

!
! Calculate Rotational Velocity of COM
! Construct the Momenta of Inertia
!
    inertia(:,:) = 0.0d0
    do i = 1, Natom
      inertia(:,:) = inertia(:,:) + calc_inertia(fictmass(i,j),ur(:,i,j)-comr(:))
    end do
!
! Obtain the Inverse Matrix of Momenta
!
    detmat = calc_determinant33(inertia)
    dinvmat(:,:) = 0.0d0

    if ( detmat <= delta ) then
      do i = 1, 3
        if ( abs(inertia(i,i)) > delta ) dinvmat(i,i) = 1/inertia(i,i)
      end do
    else
      dinvmat(:,:) = calc_inv_mat33(inertia)
    end if

    moment(:) = 0.0d0
    do i = 1, Natom
      moment(:) = moment(:) + fictmass(i,j) * crossr_product(ur(:,i,j)-comr(:), vur(:,i,j))
    end do
    ang_vel(:) = matmul(dinvmat,moment)

!
! Subtract Translation and Rotation Velocities of COM from Velocities of Atoms
!
    do i = 1, Natom
      vur(:,i,j) = vur(:,i,j) - sumvr(:)
      vur(:,i,j) &
        = vur(:,i,j) - crossr_product(ang_vel,ur(:,i,j)-comr(:))
    end do

  end do
  return

end subroutine remove_trans_rot_beads
