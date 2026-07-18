subroutine PIHMC_normal
  use Parameters
  use utility, only: program_abort, ranf1
  implicit none
  integer :: iref, Uout
  integer :: Idyn

  call Setup_time_mass
  call set_pallarel
  call Normal_Mode
  call Init_Mass

  call set_Iforce(Iforce)
  if ( Ldual ) call set_Iforce(dual_Iforce)

  if ( Lrestart .eqv. .True. ) then
    call restart_read
    call Broad3
    call nmtrans_ur2r
    call Force_New_MPI(Iforce)
    call nmtrans_fr2fur
    if ( Ldual ) call eval_potential_hmc
    call getenergy_hmc
  else
    call print_ini
    if ( MyRank == 0 ) call init_pihmc_output
    call NM_Position
    if (MyRank == 0) then
      call Init_Velocity
    end if
    call Broad3
    call nmtrans_ur2r ! x(i) = x(i) + sum_j tnm(i,j)*u(j)
    call Force_New_MPI(Iforce)
    if ( mod(istepsv,out_step) == 0 ) call print_result
    call nmtrans_fr2fur     !call Getfnm  ! fu(i) = fu(i) + sum_j fx(j)*tnm(j,i)

    if ( Ldual ) call eval_potential_hmc
    call Ham_Temp
    call print_ham(Irestep)
    call getenergy_hmc
  end if

  call getforce_ref_nor

  call save_hmc

  if ( MyRank == 0 ) then

    loop_main: &
    do istepsv = Irestep+1, Nstep

      loop_dyn: &
      do Idyn = 1, Ndyn

        call update_vel_nor
        !do iref=1, Nref   ! Nref = 5
        !  call update_vel_ref_nor
        !  call update_pos_nor
        !  call getforce_ref_nor
        !  call update_vel_ref_nor
        !end do
        call update_pos_vel_analy

        call nmtrans_ur2r       ! x(i) = x(i) + sum_j tnm(i,j)*u(j)
        call Force_New_MPI(Iforce)      ! Obtaining fx
        call nmtrans_fr2fur     ! fu(i) = fu(i) + sum_j fx(j)*tnm(j,i) !call Getfnm
        call update_vel_nor

      end do loop_dyn

      if ( Ldual ) call eval_potential_hmc
      call getenergy_hmc
      call judge_hmc

      if ( mod(istepsv,out_step) == 0 ) call print_result
      call Ham_Temp
      call print_ham(istepsv)
      call print_pihmc(istepsv)
      if ( mod(istepsv,out_step)==0 ) call restart_write(istepsv)

      if (mod(istepsv,10) == 0) then
        call exit_program
      end if

    end do loop_main

  else
    do istepsv = Irestep+1, Nstep
      do Idyn = 1, Ndyn
        call Force_New_MPI(Iforce)
      end do
      if ( Ldual ) call eval_potential_hmc
    end do
  end if

  if (MyRank == 0) then
    open(newunit=Uout,file=Fout,status='old',position='append')
      write(Uout,'(" ",a)')  repeat('*',121)
    close(Uout)
  end if

  return

contains

  subroutine init_pihmc_output
    implicit none
    integer :: Uout

    open(newunit=Uout,file=trim(dir_result)//'/pihmc.out',status='replace',form='formatted')
      write(Uout,'(a)') '# 1 Step  2 AcceptanceRatio  3 ProposalPotential  4 TargetPotential'
    close(Uout)
  end subroutine init_pihmc_output

  subroutine eval_potential_hmc
    implicit none
    real(8) :: fr_sv(Ndim,Natom,Nbead), pot_bead_sv(Nbead), potential_sv

    fr_sv = fr
    pot_bead_sv = pot_bead
    potential_sv = potential

    call Force_New_MPI(dual_Iforce)
    potential_hmc = potential
    call Virial_Estimator

    fr = fr_sv
    pot_bead = pot_bead_sv
    potential = potential_sv
  end subroutine eval_potential_hmc

  subroutine print_pihmc(Istep)
    implicit none
    integer, intent(in) :: Istep
    integer :: Uout
  if ( MyRank == 0 ) then
  if ( mod(Istep,out_step) == 0 ) then
    open(newunit=Uout,file=trim(dir_result)//'/pihmc.out',position='append')
      write(Uout,'(i10,3es24.16)') Istep, ratio, potential, potential_hmc
    close(Uout)
  end if
  end if
  end subroutine print_pihmc

  subroutine judge_hmc
    implicit none
    real(8) :: bfactor
    real(8) :: temp_rand

    bfactor = beta * (hamiltonian - hamiltonian_old)

    if ( bfactor < 75.d0 ) then
      if ( bfactor <= 0.0d0 ) then
        Naccept = Naccept + 1
      else
        temp_rand = ranf1()
        if ( exp(-bfactor) >= temp_rand ) then
          Naccept = Naccept + 1
        else
          Nreject = Nreject + 1
          call recover_hmc
        end if
      end if
    else
      Nreject = Nreject + 1
      call recover_hmc
    end if
    ratio = dble(Naccept) / dble(Naccept+Nreject)
    if (MyRank ==0 ) call Init_Velocity
    call getenergy_hmc
    call save_hmc
  end subroutine judge_hmc

  subroutine recover_hmc
    r(:,:,:)         = r_old(:,:,:)
    ur(:,:,:)        = ur_old(:,:,:)
    vur(:,:,:)       = vur_old(:,:,:)
    fr(:,:,:)        = fr_old(:,:,:)
    fur(:,:,:)       = fur_old(:,:,:)
    fur_ref(:,:,:)   = fur_ref_old(:,:,:)
    pot_bead(:)      = pot_old(:)
    potential        = potential_old
    potential_hmc    = potential_hmc_old
    E_Virial         = E_Virial_old
    dkinetic         = dkinetic_old
    qkinetic         = qkinetic_old
    hamiltonian      = hamiltonian_old
  end subroutine recover_hmc

  subroutine save_hmc
    r_old(:,:,:)       = r(:,:,:)
    ur_old(:,:,:)      = ur(:,:,:)
    vur_old(:,:,:)     = vur(:,:,:)
    fr_old(:,:,:)      = fr(:,:,:)
    fur_old(:,:,:)     = fur(:,:,:)
    fur_ref_old(:,:,:) = fur_ref(:,:,:)
    pot_old(:)         = pot_bead(:)
    potential_old      = potential
    potential_hmc_old  = potential_hmc
    E_Virial_old       = E_Virial
    dkinetic_old       = dkinetic
    qkinetic_old       = qkinetic
    hamiltonian_old    = hamiltonian
  end subroutine save_hmc

  subroutine getenergy_hmc
    use utility, only: norm_seq
    integer :: Imode, Iatom
    real(8) :: factqk
    real(8) :: get_kinetic_ene

    dkinetic = get_kinetic_ene()

    qkinetic = 0.0d0
    do Imode = 2, nbead
    do Iatom = 1, natom
      factqk = 0.5d0*dnmmass(Iatom,Imode)*omega_p2
      qkinetic = qkinetic + factqk * norm_seq( ur(:,Iatom,Imode) )
    end do
    end do

    if ( .not. Ldual ) potential_hmc = potential
    hamiltonian = dkinetic + qkinetic + potential_hmc
  end subroutine getenergy_hmc


end subroutine PIHMC_normal
