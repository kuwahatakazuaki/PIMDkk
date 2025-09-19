! Obtaining force(fx) from position(x)
! Model potential of double morse
! Atoms must be order of "O", "O", "H"
subroutine Force_Double_Morse
  use Parameters, &
    only: r, fr, Natom, Nbead, pot_bead, potential, &
          alabel, dp_inv, dir_result, istepsv, MyRank, &
          Lsave_force, physmass, dipoler, &
          AU2Ang, Ang2AU, Ista, Iend
  use utility, only: program_abort
  use mod_model_force, &
    only : harmonic, Fharmo, morse, Fmorse

  implicit none
  real(8), parameter :: lx = 0.5d0                     ! Angstrom
  real(8), parameter :: rO1(3) = [ lx, 0.0d0, 0.0d0]  * Ang2AU ! AU
  real(8), parameter :: rO2(3) = [-lx, 0.0d0, 0.0d0]  * Ang2AU ! AU

  integer :: i, j, imode, xyz
  real(8) :: f31(3), f32(3)

  fr(:,:,:) = 0.0d0
  pot_bead(:) = 0.0d0
  do imode = Ista, Iend
    f31(:)  = Fmorse(r(:,3,imode),r(:,1,imode))
    f32(:)  = Fmorse(r(:,3,imode),r(:,2,imode))
    fr(:,1,imode) = Fharmo(r(:,1,imode),rO1(:)) - f31(:)
    fr(:,2,imode) = Fharmo(r(:,2,imode),rO2(:)) - f32(:)
    fr(:,3,imode) = f31(:) + f32(:)

    pot_bead(imode) &
          = harmonic(r(:,1,imode),rO1(:)) + &
            harmonic(r(:,2,imode),rO2(:)) + &
            morse(r(:,3,imode),r(:,1,imode)) + morse(r(:,3,imode),r(:,2,imode))
  end do
  fr(:,:,Ista:Iend) = fr(:,:,Ista:Iend) * dp_inv

end subroutine Force_Double_Morse


! ++++++++++++++++++++++
! +++ Force_Harmonic +++
! ++++++++++++++++++++++
subroutine Force_Harmonic
  Use Parameters, &
    only: r, fr, Natom, Nbead, pot_bead, potential, &
          alabel, dp_inv, dir_result, istepsv, MyRank, &
          Lsave_force, physmass, dipoler, &
          AU2Ang, Ang2AU, Ista, Iend
  use utility, only: program_abort
  implicit none
  integer :: i, j, imode, xyz
  integer :: Udis, Ucent
  !real(8), parameter :: qh = 1.0d0, qo = -1.0d0
  real(8) :: f_two(3), power, dis
  real(8) :: rij(3), temp!, rij2
  real(8) :: dis_beads(Nbead)
  real(8) :: rcent(3,Natom), dipo_cent(3)

  ! Parameters in Atomic unit
  ! Hydrogen molecule with U(x) = 1/2 * mass * omega^2 * x^2
  !real(8), parameter :: r0    = 1.41014d0
  real(8), parameter :: r0    = 1.0d0 * Ang2AU
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
  pot_bead(:) = 0.0d0
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
    pot_bead(imode) = pot_bead(imode) + 0.5d0 * cons * (dis - r0)**2
    dis_beads(imode) = dis
  end do
  fr(:,:,Ista:Iend) = fr(:,:,Ista:Iend) * dp_inv
  dis_beads(:) = dis_beads(:) * AU2Ang
  ! +++ End Calculating Force in which atom (i) feels from atom (j) +++

  !!! +++ Print distance +++
  !if ( MyRank == 0 ) then
  !  open(newunit=Udis,file=trim(dir_result)//'/distance.out',status='unknown',position='append')
  !    write(Udis,*) "# ", istepsv
  !    do imode = 1, Nbead
  !      write(Udis,*) dis_beads(imode)
  !    end do
  !  close(Udis)
  !end if
  !! +++ End Print distance +++

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
    only: r, fr, Natom, Nbead, pot_bead, potential, &
          alabel, dp_inv, dir_result, istepsv, &
          Lsave_force, AU2Ang, Ang2AU
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
  open(newunit=Udis,file=trim(dir_result)//'/distance.dat',status='unknown',position='append')
    pot_bead(:) = 0.0d0
    write(Udis,*) "# ", istepsv
    do imode = 1, Nbead
      do i = 1, Natom
        do j = i+1, Natom
          rij(:) = r(:,i,imode)-r(:,j,imode)
          !dis = dsqrt( dot_product( rij(:),rij(:) ) )
          dis = norm2(rij(:))
          temp = 1 - exp(-curvat * (dis-r0))
          pot_bead(imode) = pot_bead(imode) + De * temp * temp
          write(Udis,*) dis * AU2Ang
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

