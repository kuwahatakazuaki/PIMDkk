module utility
  use Parameters, only : Ferr
  implicit none

contains

real(8) function gasdev()
  !use utility, only : ranf1
  implicit none
  real(8) :: R_Number,v1,v2,rsq,gset,fac!,gasd
  Integer          :: iset
  save iset, gset
  data iset /0/

  if (iset==0) Then
     do
       call RandomG (1,R_Number)
       v1  = 2.d0 * R_Number - 1.d0
       call RandomG (1,R_Number)
       v2  = 2.d0 * R_Number - 1.d0

       !v1  = 2.d0 * ranf1() - 1.d0
       !v2  = 2.d0 * ranf1() - 1.d0

       rsq = v1*v1 + v2*v2
       if (rsq < 1.d0 .and. rsq  > 0.d0) Then
          fac  = dsqrt(-2.d0*dlog(rsq)/rsq)
          gset   = v1 * fac
          gasdev = v2 * fac
          exit
       end if
     end do
     iset   = 1
  else
     gasdev = gset
     iset   = 0
  end if
  return
!end subroutine
end function gasdev


subroutine set_random_seed(Irand)
  integer :: Irand, i
  integer :: Nseeds
  integer, allocatable :: seeds(:)

  call random_seed(size=Nseeds)
  if ( .not. allocated(seeds) ) allocate(seeds(Nseeds))

  if ( Irand == 0 ) then
  do i = 1, Nseeds
    call system_clock(count=seeds(i))
  end do
  end if
  call random_seed(put=seeds(:))
end subroutine set_random_seed

! +++ Module for IO +++
  !function exist_file(name_file) result(bool)
  !  character(len=*), intent(in) :: name_file
  !  logical :: bool
  !end function exist_file

  subroutine read_val_next(line,key,val)
    character(len=*), intent(in) :: key, line
    real(8), intent(out) :: val
    integer :: pos
    pos = index(line,key)
    read( line(pos+1:),* ) val
  end subroutine read_val_next

  subroutine search_line(Iunit,key,line)
    integer, intent(in) :: Iunit
    character(len=*), intent(in) :: key
    character(len=*), intent(out) :: line

    do
      read(Iunit,'(a)',end=401) line
      if (index(line,key) > 0 ) exit
    end do
    return
    401 print *, 'ERROR!!: There is no key: "'//key//'"' ; stop
  end subroutine search_line
! +++ Module for IO +++

! +++ Calculation of Matrix +++
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
! +++ Calculation of Matrix +++

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
    elseif ( trim(cha) == 'al' ) then
      mass = 26.9815386d0
    elseif ( trim(cha) == 'si' ) then
      mass = 28.0855d0
    elseif ( trim(cha) == 'p' ) then
      mass = 30.97376199842d0
    elseif ( trim(cha) == 's' ) then
      mass = 32.065d0
    elseif ( trim(cha) == 'cl' ) then
      mass = 35.453d0
    elseif ( trim(cha) == 'ar' ) then
      mass = 39.948d0

    elseif ( trim(cha) == 'k' ) then
      mass = 39.0983d0
    elseif ( trim(cha) == 'ca' ) then
      mass = 40.078d0
    elseif ( trim(cha) == 'sc' ) then
      mass = 44.955912d0
    elseif ( trim(cha) == 'ti' ) then
      mass = 47.867d0
    elseif ( trim(cha) == 'v' ) then
      mass = 50.9415d0
    elseif ( trim(cha) == 'cr' ) then
      mass = 51.9961d0
    elseif ( trim(cha) == 'mn' ) then
      mass = 54.938045d0
    elseif ( trim(cha) == 'fe' ) then
      mass = 55.845d0
    elseif ( trim(cha) == 'co' ) then
      mass = 58.933195d0
    elseif ( trim(cha) == 'ni' ) then
      mass = 58.6934d0
    elseif ( trim(cha) == 'cu' ) then
      mass = 63.546d0
    elseif ( trim(cha) == 'zn' ) then
      mass = 65.38d0
    elseif ( trim(cha) == 'ga' ) then
      mass = 69.723d0
    elseif ( trim(cha) == 'ge' ) then
      mass = 72.64d0
    elseif ( trim(cha) == 'as' ) then
      mass = 74.92160d0
    elseif ( trim(cha) == 'se' ) then
      mass = 78.96d0
    elseif ( trim(cha) == 'br' ) then
      mass = 79.904d0
    elseif ( trim(cha) == 'kr' ) then
      mass = 83.798d0

    elseif ( trim(cha) == 'rb' ) then
      mass = 85.4678d0
    elseif ( trim(cha) == 'sr' ) then
      mass = 87.62d0
    elseif ( trim(cha) == 'y'  ) then
      mass = 88.90585d0
    elseif ( trim(cha) == 'zr' ) then
      mass = 91.224d0
    elseif ( trim(cha) == 'nb' ) then
      mass = 92.90638d0
    elseif ( trim(cha) == 'mo' ) then
      mass = 95.96d0
    elseif ( trim(cha) == 'tc' ) then
      mass = 98d0
    elseif ( trim(cha) == 'ru' ) then
      mass = 101.07d0
    elseif ( trim(cha) == 'rh' ) then
      mass = 102.90550d0
    elseif ( trim(cha) == 'pd' ) then
      mass = 106.42d0
    elseif ( trim(cha) == 'ag' ) then
      mass = 107.8682d0
    elseif ( trim(cha) == 'cd' ) then
      mass = 112.411d0
    elseif ( trim(cha) == 'in' ) then
      mass = 114.818d0
    elseif ( trim(cha) == 'sn' ) then
      mass = 118.710d0
    elseif ( trim(cha) == 'sb' ) then
      mass = 121.760d0
    elseif ( trim(cha) == 'te' ) then
      mass = 127.60d0
    elseif ( trim(cha) == 'i' ) then
      mass = 126.90447d0
    elseif ( trim(cha) == 'xe' ) then
      mass = 131.294d0

    elseif ( trim(cha) == 'cs' ) then
      mass = 132.9054519d0
    elseif ( trim(cha) == 'ba' ) then
      mass = 137.327d0
    elseif ( trim(cha) == 'la' ) then
      mass = 138.90547d0
    elseif ( trim(cha) == 'ce' ) then
      mass = 140.116d0

    elseif ( trim(cha) == 'pt' ) then
      mass = 195.084d0
    elseif ( trim(cha) == 'au' ) then
      mass = 196.966569d0
    elseif ( trim(cha) == 'hg' ) then
      mass = 200.59d0

    else
      call program_abort('ERROR!! "'//cha_in//'" is not exist in atom2mass')
    end if
  end function atom2mass

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
    character(len=10) :: b(3)
    call date_and_time(b(1),b(2),b(3),newtime)
    write(date_time,'(i4,a1,i2.2,a1,i2.2,a1,i2.2,a1,i2.2,a1,i2.2)') &
      newtime(1),'/',newtime(2),'/',newtime(3),' ',newtime(5),':',newtime(6),':',newtime(7)
  end function get_time

  subroutine write_err(err_messa)
    character(*), intent(in) :: err_messa
    integer :: Uerr
    open(newunit=Uerr,file=Ferr,position='append')
      write(Uerr,*) err_messa
    close(Uerr)
  end subroutine write_err

  subroutine program_abort(message)
#ifdef _mpi_
    use mpi
    character(*), intent(in) :: message
    integer :: ierr
    print *, message
    call write_err(message)
    call mpi_abort(MPI_COMM_WORLD, -1, ierr)
    stop
#else
    character(*) :: message
    print *, message
    !write_err(message)
    stop
#endif
  end subroutine program_abort

  subroutine makedir(outdir)
    character(len=*), intent(in) :: outdir
    character(len=256) command
    write(command,*) 'mkdir -p ', trim(outdir)
    write(*,*) trim(command)
    call system(command)
  end subroutine makedir

  real(8) function ranf1()
    call random_number(ranf1)
  end function ranf1


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


