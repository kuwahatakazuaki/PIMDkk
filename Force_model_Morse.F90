! Obtaining force(fx) from position(x)

! +++++++++++++++++++++++++
! +++ Force_model_Morse +++
! +++++++++++++++++++++++++
subroutine Force_model_Morse
  Use Parameters, & 
    only: r, fr, Natom, Nbead, Eenergy, potential, &
          alabel, dp_inv, address, istepsv, &
          Lsave_force, beta, &
          AUtoAng, AngtoAU
  implicit none
  integer :: i, j
  integer :: Udis, Ucoor, Ufor, Uene, imode
  real(8) :: f_two(3), power, dis !, Epoten
  real(8) :: rij(3), temp

  !! +++ Constants for conversion +++
  !real(8), parameter :: eVtoAU    = 1.0d0/27.21162
  !real(8), parameter :: AngtoAU = 1/0.529177249d0
  !real(8), parameter :: AUtoAng = 0.529177249d0
  !
  !! +++ Parameters for Morse Potential +++
  !real(8), parameter :: De    = 4.519d0  * eVtoAU   ! eV
  !real(8), parameter :: r0    = 0.74d0   * AngtoAU  ! angstrom
  !real(8), parameter :: curvat = 1.981d0 * AUtoAng  ! 1/Angstrom

  ! Parameters in Atomic unit
  ! Hydrogen molecule
    real(8), parameter :: r0    = 1.41014d0
    real(8), parameter :: De    = 0.1745d0
    real(8), parameter :: curvat = 1.0213d0

  ! +++ Calculating Forcde which atom (i) feels from atom (j) +++
  fr(:,:,:) = 0.0d0
  do imode = 1, Nbead
    do i = 1, Natom
      do j = i+1, Natom
        rij(:) = r(:,i,imode)-r(:,j,imode)
        dis = dsqrt( dot_product( rij(:),rij(:) ) )
        power = -curvat * (dis-r0)
        f_two(:) = 2 * curvat * De * exp(power) * (exp(power) - 1) * (rij(:))/dis
        fr(:,i,imode) = fr(:,i,imode) + f_two(:)
        fr(:,j,imode) = fr(:,j,imode) - f_two(:)
      end do
    end do
  end do
  fr(:,:,:) = fr(:,:,:) * dp_inv
  ! +++ End Calculating Forcde which atom (i) feels from atom (j) +++

  ! +++ Calculating enetemp +++
  open(newunit=Udis,file=trim(address)//'/distance.dat',status='unknown',position='append')
    Eenergy(:) = 0.0d0
    write(Udis,*) "# ", istepsv
    do imode = 1, Nbead
      do i = 1, Natom
        do j = i+1, Natom
          rij(:) = r(:,i,imode)-r(:,j,imode)
          dis = dsqrt( dot_product( rij(:),rij(:) ) )
          temp = 1 - exp(-curvat * (dis-r0))
          Eenergy(imode) = Eenergy(imode) + De * temp * temp
          write(Udis,*) dis * AUtoAng
        end do
      end do
    end do
  close(Udis)
  ! +++ End Calculating enetemp +++

  9998 format(3E23.15)
  9999 format(a2,1x,E15.9,1x,E15.9,1x,E15.9)
  return
end subroutine Force_model_Morse
! +++++++++++++++++++++++++++++
! +++ End Force_model_Morse +++
! +++++++++++++++++++++++++++++

! ++++++++++++++++++++++
! +++ Force_Harmonic +++
! ++++++++++++++++++++++
subroutine Force_Harmonic
  Use Parameters, & 
    only: r, fr, Natom, Nbead, Eenergy, potential, &
          alabel, dp_inv, address, istepsv, MyRank, &
          Lsave_force, beta, physmass, dipoler, &
          AUtoAng, AngtoAU, &
          ista, iend
  use utility, only: program_abort
  implicit none
  integer :: i, j, imode, xyz
  integer :: Udis, Ucent
  real(8), parameter :: qh = 1.0d0, qo = -1.0d0
  real(8) :: f_two(3), power, dis
  real(8) :: rij(3), temp!, rij2
  real(8) :: dis_beads(Nbead)
  real(8) :: rcent(3,Natom), dipo_cent(3)

  ! Parameters in Atomic unit
  ! Hydrogen molecule with U(x) = 1/2 * mass * omega^2 * x^2
  !real(8), parameter :: r0    = 1.41014d0
  real(8), parameter :: r0    = 1.0d0 * AngtoAU
  real(8), parameter :: omega  = 3000 * 0.0000045564 ! cm-1 to Hartree
  real(8), parameter :: kb  = 0.49536 ! Force constant of OH
  real(8) :: cons, reduce_mass

  if ( Natom /= 2 ) then
    call program_abort('Natom should be 2')
  end if

!  reduce_mass = (physmass(1)*physmass(2)) / (physmass(1)+physmass(2))
!  cons = reduce_mass * omega * omega
  cons = kb

  ! +++ Calculating Force in which atom (i) feels from atom (j) +++
  fr(:,:,:) = 0.0d0
  Eenergy(:) = 0.0d0
  do imode = Ista, Iend
    do i = 1, Natom
      do j = i+1, Natom
        rij(:) = r(:,i,imode)-r(:,j,imode)
        dis = norm2(rij(:))
        f_two(:) = (-1) * cons * (dis - r0) * rij(:) / dis
        fr(:,i,imode) = fr(:,i,imode) + f_two(:)
        fr(:,j,imode) = fr(:,j,imode) - f_two(:)
      end do
    end do
    dipoler(:,imode) = rij(:)
    Eenergy(imode) = Eenergy(imode) + 0.5d0 * cons * (dis - r0)**2
    dis_beads(imode) = dis
  end do
  fr(:,:,:) = fr(:,:,:) * dp_inv
  dis_beads(:) = dis_beads(:) * AUtoAng
  ! +++ End Calculating Force in which atom (i) feels from atom (j) +++

  !! +++ Centroid +++
  !do i = 1, Natom
  !  do xyz = 1, 3
  !    rcent(xyz,i) = sum(r(xyz,i,:))/dble(Nbead)
  !  end do
  !end do
  !rij(:) = rcent(:,1) - rcent(:,2)
  !dipo_cent(:) = rij(:)
  !dis = norm2(rij(:)*AUtoAng)
  !! +++ Centroid +++
  !! +++ Print distance +++
  !if ( MyRank == 0 ) then
  !  open(newunit=Udis,file=trim(address)//'/distance.out',status='unknown',position='append')
  !    write(Udis,*) "# ", istepsv
  !    do imode = 1, Nbead
  !      write(Udis,*) dis_beads(imode)
  !    end do
  !  close(Udis)

  !  open(newunit=Ucent,file=trim(address)//'/cent_dipo.out',position='append')
  !    write(Ucent,*) istepsv, dis, dipo_cent(:)
  !  close(Ucent)
  !end if
  ! +++ End Print distance +++

  9998 format(3E23.15)
  9999 format(a2,1x,E15.9,1x,E15.9,1x,E15.9)
  return
end subroutine Force_Harmonic
! ++++++++++++++++++++++++++
! +++ End Force_Harmonic +++
! ++++++++++++++++++++++++++


