subroutine init_bath
  use global_variable
  use utility
  use mpi
  implicit none
  integer :: i, j, inhc, icolor, xyz
  real(8) :: vsigma, gasd

!     /*  massive Nose-Hoover chain  */
!YK Remove the initial velocity for the centroid which is not used
  vbath(:,1:Natom,1:Nnhc,1      ) = 0.0d0
   bath(:,1:Natom,1:Nnhc,1:Nbead) = 0.0d0

  do j = 2, Nbead
    vsigma = dsqrt(1.0 / beta / qmass(j))
    do inhc = 1, Nnhc
      do i = 1, Natom
        do xyz = 1, 3
          call gasdev(gasd)
          vbath(xyz,i,inhc,j) = vsigma * gasd
        end do
      end do
    end do
  end do

  if ( Ncent == 1 ) then

    if ( Ncolor == 1) then
      do inhc = 1, Nnhc
        vsigma = dsqrt(1.0/beta/qmcent11(inhc))
        call gasdev(gasd)
        vbc11(inhc) = vsigma * gasd
        rbc11(inhc) = 0.0d0
      end do
    else if ( Ncolor > 1 ) then
      do icolor = 1, Ncolor
        do inhc = 1, Nnhc
          vsigma = dsqrt(1.0d0/beta/qmcent1(inhc,icolor))
          call gasdev(gasd)
          vbc1(inhc,icolor) = vsigma * gasd
          rbc1(inhc,icolor) = 0.0d0
        end do
      end do
    end if

  else if ( Ncent == 3 ) then

    if ( Ncolor == 1 ) then
      do inhc = 1, Nnhc
        do i = 1, Natom
          vsigma = dsqrt(1.0d0/beta/qmcent31(inhc))
          do xyz = 1, 3
            call gasdev(gasd)
            vbc31(xyz,i,inhc) = vsigma * gasd
          end do
          rbc31(:,i,inhc) = 0.0d0
        end do
      end do
    else if (Ncolor > 1) then
      do icolor = 1, Ncolor
        do inhc = 1, Nnhc
          do i = 1, Natom
            vsigma = dsqrt(1.0d0/beta/qmcent3(inhc,icolor))
            do xyz = 1, 3
              call gasdev(gasd)
              vbc3(xyz,i,inhc,icolor) = vsigma * gasd
            end do
            rbc3(:,i,inhc,icolor) = 0.0d0
          end do
        end do
      end do
    end if

  else
    call program_abort("Ncent should be 1 or 3")
  end if

!  print *, "Ncent", Ncent
!  print *, vbc11(:)
!  do inhc = 1, Nnhc
!    do i = 1, Natom
!      print *,i, vbc31(:,i,inhc)
!    end do
!    print *, ""
!  end do


end subroutine init_bath

