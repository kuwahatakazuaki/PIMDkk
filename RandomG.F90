SUBROUTINE RANDOMG(Number,RNDM1)

! This function calls the random number generator
! routine "DLARNV" in LAPACK!
!
! ISEED must be odd and in between 0 and 4095 !!!
! ISEED must have 4 dimensions!!!
!
! When N=1, generates random uniform number between (0,1)
! When N=2, generates random uniform number between (-1,1)
!
! IR is the number of the random numbers to be generated

   use Parameters
   implicit none

   integer:: Number, Number1, IRall
   double precision :: RNDM1

   IRall=1

   If(Number==0) then
      Number1=1
   Else
      Number1=Number
   End If

   Call dlarnv(Number1,Iseeds,IRall,RNDM1)

  return
end subroutine RANDOMG
