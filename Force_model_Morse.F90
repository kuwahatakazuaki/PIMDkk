! Obtaining force(fx) from position(x)

! Model potential of double morse
! Atoms must be order of "O", "O", "H"
subroutine Force_Double_Morse
  Use Parameters, & 
    only: r, fr, Natom, Nbead, Eenergy, potential, &
          alabel, dp_inv, address, istepsv, MyRank, &
          Lsave_force, physmass, dipoler, &
          AUtoAng, AngtoAU, Ista, Iend
  use utility, only: program_abort
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

  integer :: i, j, imode, xyz
  real(8) :: f31(3), f32(3)

  fr(:,:,:) = 0.0d0
  Eenergy(:) = 0.0d0
  do imode = Ista, Iend
  end do
  !f31(:) = Fmorse(r(:,3),r(:,1))
  !f32(:) = Fmorse(r(:,3),r(:,2))
  !f(:,1) = Fharmo(r(:,1),rO1(:)) - f31(:)
  !f(:,2) = Fharmo(r(:,2),rO2(:)) - f32(:)
  !f(:,3) = f31(:) + f32(:)

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
end subroutine Force_Double_Morse


! ++++++++++++++++++++++
! +++ Force_Harmonic +++
! ++++++++++++++++++++++
subroutine Force_Harmonic
  Use Parameters, & 
    only: r, fr, Natom, Nbead, Eenergy, potential, &
          alabel, dp_inv, address, istepsv, MyRank, &
          Lsave_force, physmass, dipoler, &
          AUtoAng, AngtoAU, &
          Ista, Iend
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
  fr(:,:,Ista:Iend) = fr(:,:,Ista:Iend) * dp_inv
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
  if ( MyRank == 0 ) then
    open(newunit=Udis,file=trim(address)//'/distance.out',status='unknown',position='append')
      write(Udis,*) "# ", istepsv
      do imode = 1, Nbead
        write(Udis,*) dis_beads(imode)
      end do
    close(Udis)

  !  open(newunit=Ucent,file=trim(address)//'/cent_dipo.out',position='append')
  !    write(Ucent,*) istepsv, dis, dipo_cent(:)
  !  close(Ucent)
  end if
  ! +++ End Print distance +++

  9998 format(3E23.15)
  9999 format(a2,1x,E15.9,1x,E15.9,1x,E15.9)
  return
end subroutine Force_Harmonic
! ++++++++++++++++++++++++++
! +++ End Force_Harmonic +++
! ++++++++++++++++++++++++++


! +++++++++++++++++++++++++
! +++ Force_model_Morse +++
! +++++++++++++++++++++++++
subroutine Force_model_Morse
  Use Parameters, & 
    only: r, fr, Natom, Nbead, Eenergy, potential, &
          alabel, dp_inv, address, istepsv, &
          Lsave_force, AUtoAng, AngtoAU
  implicit none
  integer :: i, j
  integer :: Udis, Ucoor, Ufor, Uene, imode
  real(8) :: f_two(3), power, dis !, Epoten
  real(8) :: rij(3), temp

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
          !dis = dsqrt( dot_product( rij(:),rij(:) ) )
          dis = norm2(rij(:))
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

