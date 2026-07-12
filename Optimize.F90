subroutine Optimize
  use Parameters
  use mod_lbfgs, only: lbfgs_init, lbfgs_step, lbfgs_free
  implicit none

  real(8) :: dx(Ndim,Natom)
  real(8) :: fmax, dmax
  integer :: iatom, niter, Uopt, Uout
  logical :: converged

  call set_pallarel
  call set_Iforce(Iforce)
  call lbfgs_init(Ndim*Natom, opt_history)

  call Force_Classical
  istepsv = 0
  niter = 0
  dmax = 0.0d0
  fmax = maximum_force()
  converged = fmax < opt_fmax

  if (MyRank == 0) then
    open(newunit=Uopt, file=trim(dir_result)//'/opt.dat', status='replace', &
         action='write', form='formatted')
    write(Uopt,'(a)') '# istep potential[Ha] fmax[Ha/bohr] dmax[Ang]'
    write(Uopt,'(i10,3es24.15)') istepsv, potential, fmax, dmax
  end if
  call print_result

  do while (niter < Nstep .and. .not. converged)
    niter = niter + 1
    istepsv = niter

    call lbfgs_step(Ndim, Natom, ur(:,:,1), -fr(:,:,1), opt_maxstep, dx)
    ur(:,:,1) = ur(:,:,1) + dx(:,:)
    dmax = 0.0d0
    do iatom = 1, Natom
      dmax = max(dmax, norm2(dx(:,iatom)))
    end do
    dmax = dmax*AU2Ang

    call Force_Classical
    fmax = maximum_force()
    converged = fmax < opt_fmax

    if (MyRank == 0) then
      write(Uopt,'(i10,3es24.15)') istepsv, potential, fmax, dmax
    end if
    if (mod(istepsv,out_step) == 0) call print_result
    if (mod(istepsv,10) == 0) call exit_program
  end do

  if (MyRank == 0) then
    close(Uopt)
    call write_final_xyz(converged)

    open(newunit=Uout, file=Fout, status='old', position='append')
    write(Uout,'(" ",a)') repeat('*',95)
    if (converged) then
      write(Uout,'(a)') ' Geometry optimization converged'
    else
      write(Uout,'(a)') ' Geometry optimization NOT converged'
    end if
    write(Uout,'(a,i0)') ' Optimization iterations: ', niter
    write(Uout,'(a,es23.15)') ' Final potential [Ha]: ', potential
    write(Uout,'(a,es23.15)') ' Final fmax [Ha/bohr]: ', fmax
    write(Uout,'(" ",a)') repeat('*',95)
    close(Uout)
  end if

  call lbfgs_free()

contains

  real(8) function maximum_force()
    integer :: i

    maximum_force = 0.0d0
    do i = 1, Natom
      maximum_force = max(maximum_force, norm2(fr(:,i,1)))
    end do
  end function maximum_force


  subroutine write_final_xyz(is_converged)
    logical, intent(in) :: is_converged
    real(8) :: xyz(3)
    integer :: i, Uxyz
    character(len=5) :: pbc_xyz

    if (Lperiodic) then
      pbc_xyz = 'T T T'
    else
      pbc_xyz = 'F F F'
    end if

    open(newunit=Uxyz, file=trim(dir_result)//'/opt_final.xyz', &
         status='replace', action='write', form='formatted')
    write(Uxyz,'(i0)') Natom
    write(Uxyz,'(a,a,a,es23.15,a,a)') &
      'Properties=species:S:1:pos:R:3 pbc="', pbc_xyz, '" energy=', &
      potential, ' converged=', merge('T','F',is_converged)
    do i = 1, Natom
      xyz(:) = 0.0d0
      xyz(1:Ndim) = ur(:,i,1)*AU2Ang
      write(Uxyz,'(a,3(1x,es23.15))') trim(alabel(i)), xyz(:)
    end do
    close(Uxyz)
  end subroutine write_final_xyz

end subroutine Optimize
