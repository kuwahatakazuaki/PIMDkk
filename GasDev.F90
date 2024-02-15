Subroutine gasdev(gasd)
  implicit none
  Double Precision :: R_Number,v1,v2,rsq,gset,fac,gasd
  Integer          :: iset
  save iset, gset
  data iset /0/

  If (iset==0)Then
     Do
        Call RandomG (1,R_Number)
        v1  = 2.d0 * R_Number - 1.d0
        Call RandomG (1,R_Number)
        v2  = 2.d0 * R_Number - 1.d0
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
