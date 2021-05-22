!subroutine tune_tempe
!  use global_variable
!  use utility
!  implicit none
!  real(8) :: tempe, scalefac
!
!!  call kinetic_energy(tempe)
!  tempe = kinetic_energy()
!  tempe = tempe / (3.0d0 * dble(Natom)) / dble(Nbead) / KtoAU
!  scalefac = dsqrt( temperature / tempe ) 
!  print *, "tempe",tempe
!  vu(:,:,:) = vu(:,:,:) * scalefac
!end subroutine tune_tempe

!subroutine kinetic_energy(kine)
!  use global_variable
!  real(8), intent(out) :: kine
!  integer :: i, j
!  kine = 0.0d0
!  do i = 1, Natom
!    do j = 1, Nbead
!      kine = kine + fictmass(i,j) * dot_product(vu(:,i,j),vu(:,i,j))
!    end do
!  end do
!! --- Be care!! Original program doesn't have the factor of 0.5 ---
!  kine = 0.5 * kine
!end subroutine kinetic_energy

!Subroutine Temp_ctr
!
!  Use Parameters
!  Implicit None
!  Double Precision :: temp_scale, tempi
!
!  tempi=0.0D+00
!  call kinetic_energy(tempi)
!  tempi = 0.5d0*tempi
!  tempi = 2.d0*tempi/(3.d0*dble(natom))/boltz
!  tempi = tempi/dble(nbead)
!!    Write(*,*) 'Temperature calculated from Initial Velocities = ',tempi      
!!    Write(*,*) 'Rescale the Velocities to Set to Temperature   = ',temperature
!
!  temp_scale = dsqrt(temperature/tempi)
!!      Write(*,*) 'Rescale Paramter   = ',temp_scale
!  do imode = 1, nbead
!     do iatom = 1, natom
!        vux(iatom,imode) = vux(iatom,imode)*temp_scale
!        vuy(iatom,imode) = vuy(iatom,imode)*temp_scale
!        vuz(iatom,imode) = vuz(iatom,imode)*temp_scale
!     enddo
!  enddo
!
!Return
!End Subroutine
