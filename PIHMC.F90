subroutine PIHMC_normal
  use Parameters
  use utility, only: program_abort, ranf1
  implicit none
  integer :: iref, Uout
  integer :: Idyn, Ndyn = 20

  call Setup_time_mass
  call set_pallarel
  call Normal_Mode
  call Init_Mass

  call set_Iforce

  if ( Lrestart .eqv. .True. ) then
    !if (MyRank == 0) Then
    call restart_read
    !end if
    call Broad3
  else
    call print_ini
    call NM_Position
    if (MyRank == 0) then
      call Init_Velocity
    end if
    call Broad3
    call Temp_ctr
    call nmtrans_ur2r ! x(i) = x(i) + sum_j tnm(i,j)*u(j)
    call Force_New_MPI
    if ( mod(istepsv,out_step) == 0 ) call print_result_qm
    call nmtrans_fr2fur     !call Getfnm  ! fu(i) = fu(i) + sum_j fx(j)*tnm(j,i)

    !if ( MyRank == 0 ) then
    call Ham_Temp
    call print_ham(Irestep)
    !end if
  end if

  call getforce_ref_nor

  call save_hmc

  if ( MyRank == 0 ) then

    loop_main: &
    do istepsv = Irestep+1, nstep

      loop_dyn: &
      do Idyn = 1, Ndyn

        call update_vel_nor
        do iref=1, Nref   ! Nref = 5
          call update_vel_ref_nor
          call update_pos_nor
          call getforce_ref_nor
          call update_vel_ref_nor
        end do
        call nmtrans_ur2r       ! x(i) = x(i) + sum_j tnm(i,j)*u(j)
        call Force_New_MPI      ! Obtaining fx
        call nmtrans_fr2fur     ! fu(i) = fu(i) + sum_j fx(j)*tnm(j,i) !call Getfnm
        call update_vel_nor

      end do loop_dyn

      call getenergy_hmc
      call judge_hmc

      if ( mod(istepsv,out_step) == 0 ) call print_result_qm
      call Ham_Temp
      call print_ham(istepsv)
      call print_pihmc(istepsv)
      if ( mod(istepsv,out_step)==0 ) call Restart_Write(istepsv)

      if (mod(istepsv,10) == 0) then
        call exit_program
      end if

    end do loop_main

  else
    do istepsv = Irestep+1, nstep
      call Force_New_MPI
    end do
  end if

  if (MyRank == 0) then
    open(newunit=Uout,file=Fout,status='old',position='append')
      write(Uout,'(" ",a)')  repeat('*',121)
    close(Uout)
  end if

  return

contains

  subroutine print_pihmc(Istep)
    implicit none
    integer, intent(in) :: Istep
    integer :: Uout
  if ( MyRank == 0 ) then
    open(newunit=Uout,file=trim(dir_result)//'/pihmc.out',position='append')
      write(Uout,*) Istep, ratio
    close(Uout)
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
!print *, istepsv, bfactor, exp(-bfactor), temp_rand
        !if ( exp(-bfactor) >= ranf1() ) then
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
!print *, istepsv, bfactor, Naccept, Nreject, ratio
    if (MyRank ==0 ) call Init_Velocity
    !call Broad_velocity
    call Temp_ctr
    call getenergy_hmc
    call save_hmc
  end subroutine judge_hmc

  subroutine recover_hmc
    ur(:,:,:)        = ur_old(:,:,:)
    vur(:,:,:)       = vur_old(:,:,:)
    fur(:,:,:)       = fur_old(:,:,:)
    fur_ref(:,:,:)   = fur_ref_old(:,:,:)
    pot_bead(:)      = pot_old(:)
    hamiltonian      = hamiltonian_old
  end subroutine recover_hmc

  subroutine save_hmc
    ur_old(:,:,:)      = ur(:,:,:)
    vur_old(:,:,:)     = vur(:,:,:)
    fur_old(:,:,:)     = fur(:,:,:)
    fur_ref_old(:,:,:) = fur_ref(:,:,:)
    pot_old(:)         = pot_bead(:)
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

    hamiltonian = dkinetic + qkinetic + potential
  end subroutine getenergy_hmc


end subroutine PIHMC_normal
