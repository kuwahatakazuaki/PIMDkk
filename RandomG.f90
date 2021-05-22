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
   integer, dimension(4) :: ISEED 
!   save ISEED
     
   IRall=1

!  WRITE(*,*) 'HERE?'

   ISEED(1)=ISEED1
   ISEED(2)=ISEED2
   ISEED(3)=ISEED3
   ISEED(4)=ISEED4

   If(Number==0) then
      Number1=1
   Else
      Number1=Number
   End If

   Call dlarnv(Number1,ISEED,IRall,RNDM1)

   ISEED1 = ISEED(1)
   ISEED2 = ISEED(2)
   ISEED3 = ISEED(3)
   ISEED4 = ISEED(4)

  return
end subroutine RANDOMG
