module utility
  implicit none

contains

  function calc_inv_mat33(mat) result(inv)
    real(8), intent(in) :: mat(3,3)
    real(8) :: inv(3,3), det
    real(8) :: a11, a12, a13
    real(8) :: a21, a22, a23
    real(8) :: a31, a32, a33
    a11 = mat(1,1); a12 = mat(1,2); a13 = mat(1,3)
    a21 = mat(2,1); a22 = mat(2,2); a23 = mat(2,3)
    a31 = mat(3,1); a32 = mat(3,2); a33 = mat(3,3)
    inv(1,1) =  a22*a33 - a23*a32
    inv(1,2) = -a12*a33 + a13*a32
    inv(1,3) =  a12*a23 - a13*a22
    inv(2,1) = -a21*a33 + a23*a31
    inv(2,2) =  a11*a33 - a13*a31
    inv(2,3) = -a11*a23 + a13*a21
    inv(3,1) =  a21*a32 - a22*a31
    inv(3,2) = -a11*a32 + a12*a31
    inv(3,3) =  a11*a22 - a12*a21
    det = calc_determinant33(mat)
    inv(:,:) = inv(:,:) / det
  end function calc_inv_mat33

  function calc_determinant33(mat) result(det)
    real(8), intent(in) :: mat(3,3)
    real(8) :: det
    det = mat(1,1)*mat(2,2)*mat(3,3) - mat(1,1)*mat(2,3)*mat(3,2) &
        + mat(1,2)*mat(2,3)*mat(3,1) - mat(1,2)*mat(2,1)*mat(3,3) &
        + mat(1,3)*mat(2,1)*mat(3,2) - mat(1,3)*mat(2,2)*mat(3,1)
  end function calc_determinant33

  function outer_product(a,b) result(vec)
    real(8), intent(in) :: a(3), b(3)
    real(8) :: vec(3)

    vec(1) = a(2) * b(3) - a(3) * b(2)
    vec(2) = a(3) * b(1) - a(1) * b(3)
    vec(3) = a(1) * b(2) - a(2) * b(1)
  end function outer_product

  real(8) function norm_seq(x)
    implicit none
    real(8), intent(in) :: x(:)
    norm_seq = dot_product(x(:),x(:))
  end function

  function atom2mass(cha_in) result(mass)
    implicit none
    character(*), intent(in) :: cha_in
    character(len(cha_in)) :: cha
    real(8) :: mass
    cha = lowerchr(trim(cha_in))
    mass = 0
    if     ( trim(cha) == 'h' ) then
      mass = 1.00782503223d0
    elseif ( trim(cha) == 'd' ) then
      mass = 2.01410177812d0
    elseif ( trim(cha) == 'mu' ) then
      mass = 0.1134289257d0
    elseif ( trim(cha) == 'he' ) then
      mass = 4.00260325413d0
    elseif ( trim(cha) == 'li' ) then
      mass = 7.0160034366d0
    elseif ( trim(cha) == 'be' ) then
      mass = 9.012183065d0
    elseif ( trim(cha) == 'b' ) then
      mass = 6.9675d0
    elseif ( trim(cha) == 'c' ) then
      mass = 12.0d0
    elseif ( trim(cha) == 'n' ) then
      mass = 14.00307400443d0
    elseif ( trim(cha) == 'o' ) then
      mass = 15.99491461957d0
    elseif ( trim(cha) == 'f' ) then
      mass = 18.99840316273d0
    elseif ( trim(cha) == 'ne' ) then
      mass = 19.9924401762d0
    elseif ( trim(cha) == 'na' ) then
      mass = 22.9897692820d0
    elseif ( trim(cha) == 'mg' ) then
      mass = 24.3055d0
    elseif ( trim(cha) == 'si' ) then
      mass = 27.97692653465d0
    elseif ( trim(cha) == 'p' ) then
      mass = 30.97376199842d0
    elseif ( trim(cha) == 's' ) then
      mass = 31.9720711744d0
    elseif ( trim(cha) == 'cl' ) then
      mass = 35.4515d0
    else
      call program_abort('ERROR!! "'//cha_in//'" is not exist in atom2mass')
    end if
  end function atom2mass

  !function lowerchr(str)
  !  character(*), intent(in) :: str
  !  character(len(str)) :: lowerchr
  !  integer :: i
  !  do i = 1, len_trim(str)
  !    if ( str(i:i) >= 'A' .and. str(i:i) <= 'Z' ) then
  !      lowerchr(i:i) = char(ichar(str(i:i))+32)
  !    else
  !      lowerchr(i:i) = str(i:i)
  !    end if
  !  end do
  !end function lowerchr

  function lowerchr(str_in) result(str_out)
    character(len=*), intent(in) :: str_in
    character(len=len(str_in))   :: str_out
    integer, parameter :: ilowerA = ichar('a')
    integer, parameter :: iupperA = ichar('A')
    integer, parameter :: iupperZ = ichar('Z')

    integer :: i, ichr, nchr, iconv

    iconv = ilowerA - iupperA

    nchr = len(str_in)
    do i = 1, nchr
       ichr = ichar(str_in(i:i))
       if ((ichr >= iupperA) .and. (ichr <= iupperZ)) then
          str_out(i:i) = char(ichr + iconv)
       else
          str_out(i:i) = str_in(i:i)
       end if
    end do
  end function lowerchr

  !--------------------------------------------------------------------!

  function upperchr(str_in) result(str_out)
    character(len=*), intent(in) :: str_in
    character(len=len(str_in))   :: str_out
    integer, parameter :: ilowerA = ichar('a')
    integer, parameter :: ilowerZ = ichar('z')
    integer, parameter :: iupperA = ichar('A')

    integer :: i, ichr, nchr, iconv

    iconv = iupperA - ilowerA

    nchr = len(str_in)
    do i = 1, nchr
       ichr = ichar(str_in(i:i))
       if ((ichr >= ilowerA) .and. (ichr <= ilowerZ)) then
          str_out(i:i) = char(ichr + iconv)
       else
          str_out(i:i) = str_in(i:i)
       end if
    end do
  end function upperchr

  function get_time() result(date_time)
    character(len=19) :: date_time
    integer :: newtime(8)
    character :: b(3)
    call date_and_time(b(1),b(2),b(3),newtime)
    write(date_time,'(i4,a1,i2.2,a1,i2.2,a1,i2.2,a1,i2.2,a1,i2.2)') &
      newtime(1),'/',newtime(2),'/',newtime(3),' ',newtime(5),':',newtime(6),':',newtime(7)
  end function get_time

  subroutine program_abort(message)
#ifdef _mpi_
    use mpi
    character(*) :: message
    integer :: ierr
    print *, message
    call mpi_abort(MPI_COMM_WORLD, -1, ierr)
    stop
#else
    character(*) :: message
    print *, message
    stop
#endif
  end subroutine program_abort

  !subroutine get_inv_mat(mat,inv,n)
  !  integer :: n
  !  real(8), intent(in)  :: mat(n,n)
  !  real(8), intent(out) :: inv(n,n)
  !  integer :: lwork, lda, info
  !  real(8), allocatable :: work(:)
  !  integer, allocatable :: ipiv(:)
  !  inv(:,:) = mat(:,:)
  !  lda = n
  !  lwork = 64*n
  !  allocate(work(lwork),ipiv(n))
  !  call dgetrf(N, N, inv, lda, ipiv, info)
  !  call dgetri(N, inv, lda, ipiv, work, lwork, info)
  !end subroutine get_inv_mat

end module utility

!  subroutine program_abort(message)
!#ifdef _mpi_
!    use mpi
!#endif
!    character(*) :: message
!    integer :: ierr
!    print *, message
!#ifdef _mpi_
!    call mpi_abort(MPI_COMM_WORLD, -1, ierr)
!#else
!    stop
!#endif
!  end subroutine program_abort

  !character(len=2) function itoc(i)
  !  integer :: i
  !  select case(i)
  !    case(1)
  !      itoc = 'H'
  !    case(3)
  !      itoc = 'Li'
  !    case(5)
  !      itoc = 'B'
  !    case(6)
  !      itoc = 'C'
  !    case(7)
  !      itoc = 'N'
  !    case(8)
  !      itoc = 'O'
  !    case(9)
  !      itoc = 'F'
  !  end select
  !end function itoc

!  subroutine gasdev(gasd)
!    real(8), intent(inout) :: gasd
!    integer, save :: iset = 0
!    real(8) :: v1, v2, rsq, fac
!    real(8), save :: gset
!    real(8) :: ran
!    integer :: i
!
!    if ( iset == 0 ) then
!      do
!        call random_generator(1,ran)
!        v1 = 2.0 * ran - 1.0
!        call random_generator(1,ran)
!        v2 = 2.0 * ran - 1.0
!        rsq = v1 * v1 + v2 * v2
!!        if ( rsq > 0.0d0 .and. rsq < 1.0d0 ) then
!        if ( rsq < 1.0d0 ) then
!          fac = dsqrt( -2.0 * dlog(rsq) / rsq )
!          gset = v1 * fac
!          gasd = v2 * fac
!          exit
!        end if
!      end do
!      iset = 1
!    else
!      gasd = gset
!      iset = 0
!    end if
!  end subroutine gasdev

!  subroutine set_input_file(Ifile,Def_file)
!    integer :: leng
!    character(:), allocatable, intent(out):: Ifile
!    character(*), optional :: Def_file
!    if ( command_argument_count() == 0) then
!      Ifile = Def_file
!!      print *, "Reading from default file as ", Ifile
!!      print *, 'ERROR!! There is no input file'
!    else
!      call get_command_argument(1, length=leng)
!        allocate(character(leng) :: Ifile)
!        call get_command_argument(1, Ifile)
!!      print *, "Reading from ", Ifile
!    end if
!    return
!  end subroutine set_input_file

!  subroutine random_seed_clock()
!    integer :: nseed, clock
!    integer, allocatable :: seed(:)
!
!    call system_clock(clock)
!
!    call random_seed(size=nseed)
!    allocate(seed(nseed))
!    seed = clock
!    call random_seed(put=seed)
!  end subroutine random_seed_clock



