subroutine init_velocity
  use global_variable
  use utility
  implicit none
  real(8) :: gasd, vsigma
  integer :: imode, iatom
  integer :: i,j

!
!     /*  vsigma: standard deviation of Maxwell distribution  */
!

  do imode = 1, nbead
    do iatom = 1, natom
      vsigma = dsqrt(1.d0/beta/fictmass(iatom,imode))
      do i = 1, 3
        call gasdev(gasd)
        vu(i,iatom,imode) = vsigma*gasd
      end do
    enddo
  enddo

  call remove_rotation
  call tune_tempe

return
contains

subroutine remove_rotation
!  use global_variable ! u(3,Natom,Nbead), vu(3,Natom,Nbead)
!  implicit none
  integer :: i, j, xyz
  integer :: iatom, imode
  real(8) :: sumv(3), com(3), ave(3), grad(3)
  real(8) :: ci11, ci22, ci33, ci12, ci23, ci31
  real(8) :: detmat, totmas
  real(8) :: dinvmat(3,3)
  real(8) :: tempu(3,Natom), tempv(3,Natom)

!
!   Remove Translation and Rotation Velocity of 
!   COM for Centroids and Noncentroids
!
  do j = 1, Nbead
    sumv(:) = 0.0d0
    com(:)  = 0.0d0
    ave(:)  = 0.0d0
    ci11    = 0.0d0
    ci22    = 0.0d0
    ci33    = 0.0d0
    ci12    = 0.0d0
    ci23    = 0.0d0
    ci31    = 0.0d0
    totmas  = 0.0d0

    do i = 1, Natom
      sumv(:) = sumv(:) + vu(:,i,j) * fictmass(i,j)
      com(:)  = com(:)  +  u(:,i,j) * fictmass(i,j)
    end do

    totmas = sum(fictmass(:,j))
    sumv(:) = sumv(:) / totmas
    com(:)  = com(:)  / totmas


    do xyz = 1, 3
      tempu(xyz,:) = u(xyz,:,j) - com(xyz)
    end do

    do i = 1, Natom
      ci11 = ci11 + fictmass(i,j) * ( tempu(2,i)**2 + tempu(3,i)**2 )
      ci22 = ci22 + fictmass(i,j) * ( tempu(3,i)**2 + tempu(1,i)**2 )
      ci33 = ci33 + fictmass(i,j) * ( tempu(1,i)**2 + tempu(2,i)**2 )
      ci12 = ci12 - fictmass(i,j) *   tempu(1,i)    * tempu(2,i)
      ci23 = ci23 - fictmass(i,j) *   tempu(2,i)    * tempu(3,i)
      ci31 = ci31 - fictmass(i,j) *   tempu(3,i)    * tempu(1,i)
    end do



    detmat = ci11*ci22*ci33 + 2.0*ci12*ci23*ci31 - ci11*ci23*ci23 - ci31*ci22*ci31 - ci12*ci12*ci33
    if ( detmat <= 0.00001d0 ) goto 100

!      detmat_inv = 1.0d0 / detmat
      dinvmat(1,:) = [ ci22*ci33-ci23*ci23 , ci23*ci31-ci33*ci12 , ci12*ci23-ci31*ci22 ]
      dinvmat(2,:) = [ dinvmat(1,2)        , ci11*ci33-ci31*ci31 , ci31*ci12-ci11*ci23 ]
      dinvmat(3,:) = [ dinvmat(1,3)        , dinvmat(2,3)        , ci11*ci22-ci12*ci12 ]
      dinvmat(:,:) = dinvmat(:,:) / detmat

    do i = 1, Natom
      grad(1) = tempu(2,i)*vu(3,i,j) - tempu(3,i)*vu(2,i,j)
      grad(2) = tempu(3,i)*vu(1,i,j) - tempu(1,i)*vu(3,i,j)
      grad(3) = tempu(1,i)*vu(2,i,j) - tempu(2,i)*vu(1,i,j)
      do xyz = 1, 3
        ave(xyz) = ave(xyz) + fictmass(i,j) * dot_product(dinvmat(xyz,:),grad(:))
      end do
    end do

100 continue


    do i = 1, Natom
      vu(1,i,j) = vu(1,i,j) - sumv(1) - ave(2)*tempu(3,i) + ave(3)*tempu(2,i)
      vu(2,i,j) = vu(2,i,j) - sumv(2) - ave(3)*tempu(1,i) + ave(1)*tempu(3,i)
      vu(3,i,j) = vu(3,i,j) - sumv(3) - ave(1)*tempu(2,i) + ave(2)*tempu(1,i)
    end do
  end do
end subroutine remove_rotation

subroutine tune_tempe
!  use global_variable
!  use utility
!  implicit none
  real(8) :: tempe, scalefac

!  call kinetic_energy(tempe)
  tempe = kinetic_energy()
  tempe = 2.0d0 * tempe / (3.0d0 * dble(Natom)) / dble(Nbead) / KtoAU
  scalefac = dsqrt( temperature / tempe ) 
  vu(:,:,:) = vu(:,:,:) * scalefac
end subroutine tune_tempe

end subroutine
