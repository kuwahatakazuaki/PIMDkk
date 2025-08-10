subroutine gasdev(gasd)
  use utility, only : ranf1
  implicit none
  Double Precision :: R_Number,v1,v2,rsq,gset,fac,gasd
  Integer          :: iset
  save iset, gset
  data iset /0/

  if (iset==0) Then
     do
        Call RandomG (1,R_Number)
        v1  = 2.d0 * R_Number - 1.d0
        Call RandomG (1,R_Number)
        v2  = 2.d0 * R_Number - 1.d0

        !v1  = 2.d0 * ranf1() - 1.d0
        !v2  = 2.d0 * ranf1() - 1.d0

        rsq = v1*v1 + v2*v2
        If (rsq < 1.d0 .and. rsq  > 0.d0) Then
           fac  = dsqrt(-2.d0*dlog(rsq)/rsq)
           gset   = v1 * fac
           gasd   = v2 * fac
           Exit
        Endif
     Enddo
     iset   = 1
  Else
     gasd   = gset
     iset   = 0
  End if
  Return
End Subroutine


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

