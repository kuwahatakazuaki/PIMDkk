module mod_lbfgs
  implicit none
  private

  public :: lbfgs_init, lbfgs_free, lbfgs_step

  real(8), parameter :: h0_diag = 1.0d0 / 0.720d0

  real(8), allocatable :: s_hist(:,:), y_hist(:,:), rho_hist(:)
  real(8), allocatable :: x_prev(:), g_prev(:)
  integer :: n_alloc = 0
  integer :: nhist = 0
  integer :: m_max = 0
  integer :: hist_start = 1
  integer :: ncall = 0

contains

  subroutine lbfgs_init(n, m)
    integer, intent(in) :: n, m

    if (n <= 0) error stop 'lbfgs_init: n must be positive'
    if (m <= 0) error stop 'lbfgs_init: m must be positive'

    call lbfgs_free()

    allocate(s_hist(n,m), y_hist(n,m), rho_hist(m))
    allocate(x_prev(n), g_prev(n))

    s_hist(:,:) = 0.0d0
    y_hist(:,:) = 0.0d0
    rho_hist(:) = 0.0d0
    x_prev(:) = 0.0d0
    g_prev(:) = 0.0d0

    n_alloc = n
    m_max = m
    nhist = 0
    hist_start = 1
    ncall = 0
  end subroutine lbfgs_init


  subroutine lbfgs_free()
    if (allocated(s_hist)) deallocate(s_hist)
    if (allocated(y_hist)) deallocate(y_hist)
    if (allocated(rho_hist)) deallocate(rho_hist)
    if (allocated(x_prev)) deallocate(x_prev)
    if (allocated(g_prev)) deallocate(g_prev)

    n_alloc = 0
    m_max = 0
    nhist = 0
    hist_start = 1
    ncall = 0
  end subroutine lbfgs_free


  subroutine lbfgs_step(Ndim, Natom, x, g, maxstep, dx)
    integer, intent(in) :: Ndim, Natom
    real(8), intent(in) :: x(Ndim,Natom), g(Ndim,Natom)
    real(8), intent(in) :: maxstep
    real(8), intent(out) :: dx(Ndim,Natom)

    real(8), allocatable :: x_flat(:), g_flat(:), q(:), r(:)
    real(8), allocatable :: s_new(:), y_new(:), alpha(:)
    integer, allocatable :: slots(:)
    real(8) :: beta, gamma, sy, yy, max_atom_step, scale
    integer :: iatom, ihist, idx, n

    if (.not. allocated(s_hist)) then
      error stop 'lbfgs_step: lbfgs_init has not been called'
    end if
    if (Ndim <= 0 .or. Natom <= 0) then
      error stop 'lbfgs_step: dimensions must be positive'
    end if
    if (maxstep <= 0.0d0) then
      error stop 'lbfgs_step: maxstep must be positive'
    end if

    n = Ndim*Natom
    if (n /= n_alloc) then
      error stop 'lbfgs_step: dimensions differ from lbfgs_init'
    end if

    allocate(x_flat(n), g_flat(n), q(n), r(n))
    x_flat(:) = reshape(x, [n])
    g_flat(:) = reshape(g, [n])

    if (ncall >= 1) then
      allocate(s_new(n), y_new(n))
      s_new(:) = x_flat(:) - x_prev(:)
      y_new(:) = g_flat(:) - g_prev(:)
      sy = dot_product(s_new, y_new)

      if (sy > 0.0d0) then
        if (nhist < m_max) then
          idx = modulo(hist_start + nhist - 1, m_max) + 1
          nhist = nhist + 1
        else
          idx = hist_start
          hist_start = modulo(hist_start, m_max) + 1
        end if
        s_hist(:,idx) = s_new(:)
        y_hist(:,idx) = y_new(:)
        rho_hist(idx) = 1.0d0/sy
      end if
    end if

    q(:) = g_flat(:)
    if (nhist > 0) then
      allocate(alpha(nhist), slots(nhist))
      do ihist = 1, nhist
        slots(ihist) = modulo(hist_start + ihist - 2, m_max) + 1
      end do

      do ihist = nhist, 1, -1
        idx = slots(ihist)
        alpha(ihist) = rho_hist(idx)*dot_product(s_hist(:,idx), q)
        q(:) = q(:) - alpha(ihist)*y_hist(:,idx)
      end do

      idx = slots(nhist)
      sy = dot_product(s_hist(:,idx), y_hist(:,idx))
      yy = dot_product(y_hist(:,idx), y_hist(:,idx))
      if (yy > 0.0d0) then
        gamma = sy/yy
      else
        gamma = h0_diag
      end if
      r(:) = gamma*q(:)

      do ihist = 1, nhist
        idx = slots(ihist)
        beta = rho_hist(idx)*dot_product(y_hist(:,idx), r)
        r(:) = r(:) + s_hist(:,idx)*(alpha(ihist) - beta)
      end do
      r(:) = -r(:)
    else
      r(:) = -h0_diag*g_flat(:)
    end if

    if (dot_product(r, g_flat) >= 0.0d0) then
      nhist = 0
      hist_start = 1
      r(:) = -h0_diag*g_flat(:)
    end if

    dx(:,:) = reshape(r, [Ndim,Natom])
    max_atom_step = 0.0d0
    do iatom = 1, Natom
      max_atom_step = max(max_atom_step, norm2(dx(:,iatom)))
    end do
    if (max_atom_step > maxstep) then
      scale = maxstep/max_atom_step
      dx(:,:) = scale*dx(:,:)
    end if

    x_prev(:) = x_flat(:)
    g_prev(:) = g_flat(:)
    ncall = ncall + 1
  end subroutine lbfgs_step

end module mod_lbfgs
