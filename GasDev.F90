!subroutine gasdev(gasd)
real(8) function gasdev()
  use utility, only : ranf1
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

