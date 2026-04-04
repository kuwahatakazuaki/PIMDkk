subroutine restart_read
  use Parameters
  implicit none
  integer :: Uinp
  integer :: i, j, Inhc
  real(8) :: pur(3)

  if ( MyRank == 0 ) then
    open(newunit=Uinp, file=trim(dir_result)//'/restart.dat', status = 'unknown')
      read(Uinp,*) Irestep
      do j = 1, Nbead
        do i = 1, Natom
          read(Uinp,*) ur(:,i,j)
        end do
      end do

      do j = 1, Nbead
        do i = 1, Natom
          read(Uinp,*) pur(:)
          vur(:,i,j) = pur(:) / sqrt(fictmass(i,j))
        end do
      end do

      do j = 1, Nbead
        do i = 1, Natom
          read(Uinp,*) fur(:,i,j)
        end do
      end do

      if ( allocated(rbath) ) then
        do j = 1, Nbead
          do Inhc = 1, Nnhc
            do i = 1, Natom
              read(Uinp,*) rbath(:,i,Inhc,j)
            enddo
          enddo
        enddo

        do j = 1, Nbead
          do Inhc  = 1, Nnhc
            do i = 1, Natom
              read(Uinp,*) vrbath(:,i,Inhc,j)
            enddo
          enddo
        enddo

        do j = 1, Nbead
          do Inhc  = 1, Nnhc
            do i = 1, Natom
              read(Uinp,*) frbath(:,i,Inhc,j)
            enddo
          enddo
        enddo
      end if

      select case(Ncent)
        case(0)
          continue
        case(1)
          do Inhc=1,Nnhc
            read(Uinp,*) rbc11(Inhc),vbc11(Inhc),fbc11(Inhc)
          enddo
        case(3)
          do Inhc=1,Nnhc
            do i=1,Natom
              read(Uinp,*) rbc31(:,i,Inhc)
            enddo
          enddo
          do Inhc=1,Nnhc
            do i=1,Natom
              read(Uinp,*) vrbc31(:,i,Inhc)
            enddo
          enddo
          do Inhc=1,Nnhc
            do i=1,Natom
              read(Uinp,*) frbc31(:,i,Inhc)
            enddo
          enddo
      end select
    close(Uinp)
  end if
  return
end subroutine restart_read

subroutine restart_write(istep)
  use Parameters
  implicit none
  integer :: Istep, Uout
  integer :: i, j, Inhc
  logical :: ex_restart
  character(len=256) :: frestart, frestart1, frestart_tmp

  if ( MyRank == 0 ) then
    frestart     = trim(dir_result)//'/restart.dat'
    frestart1    = trim(dir_result)//'/restart1.dat'
    frestart_tmp = trim(dir_result)//'/restart.tmp'

    open(newunit=Uout, file=trim(frestart_tmp), status='replace')
      write(Uout,'(i10)') Istep
      do j = 1, Nbead
        do i = 1, Natom
          write(Uout,*) ur(:,i,j)
        end do
      end do

      do j = 1, Nbead
        do i = 1, Natom
          write(Uout,*) vur(:,i,j) * sqrt(fictmass(i,j))
        end do
      end do

      do j = 1, Nbead
        do i = 1, Natom
          write(Uout,*) fur(:,i,j)
        end do
      end do

      if ( allocated(rbath) ) then
        do j = 1, Nbead
          do Inhc = 1, Nnhc
            do i = 1, Natom
              write(Uout,*) rbath(:,i,Inhc,j)
            enddo
          enddo
        enddo
        do j = 1, Nbead
          do Inhc  = 1, Nnhc
            do i = 1, Natom
              write(Uout,*) vrbath(:,i,Inhc,j)
            enddo
          enddo
        enddo
        do j = 1, Nbead
          do Inhc  = 1, Nnhc
            do i = 1, Natom
              write(Uout,*) frbath(:,i,Inhc,j)
            enddo
          enddo
        enddo
      end if

      select case(Ncent)
        case(1)
          do Inhc=1,Nnhc
            write(Uout,*) rbc11(Inhc),vbc11(Inhc),fbc11(Inhc)
          enddo
        case(3)
          do Inhc=1,Nnhc
            do i=1,Natom
              write(Uout,*) rbc31(:,i,Inhc)
            enddo
          enddo
          do Inhc=1,Nnhc
            do i=1,Natom
              write(Uout,*) vrbc31(:,i,Inhc)
            enddo
          enddo
          do Inhc=1,Nnhc
            do i=1,Natom
              write(Uout,*) frbc31(:,i,Inhc)
            enddo
          enddo
      end select
    close(Uout)

    inquire(file=trim(frestart), exist=ex_restart)
    if (ex_restart .and. istep > out_step) then
      call system('mv -f "'//trim(frestart)//'" "'//trim(frestart1)//'"')
    end if
    call system('mv -f "'//trim(frestart_tmp)//'" "'//trim(frestart)//'"')
  end if

return
end subroutine restart_write


subroutine restart_write_bin(istep)
  use Parameters
  implicit none
  integer :: Istep, Uout
  integer :: i, j, Inhc
  logical :: ex_restart
  character(len=256) :: frestart, frestart1, frestart_tmp

  if ( MyRank == 0 ) then
    frestart     = trim(dir_result)//'/restart.dat'
    frestart1    = trim(dir_result)//'/restart1.dat'
    frestart_tmp = trim(dir_result)//'/restart.tmp'

    !open(newunit=Uout, file=trim(frestart_tmp), status='replace')
    open(newunit=Uout, file=trim(frestart_tmp), &
       form='unformatted', access='stream', status='replace')
      write(Uout) Istep
      write(Uout) ur
      write(Uout) vur * spread(sqrt(fictmass), dim=1, ncopies=3)
      write(Uout) fur

      if (allocated(rbath)) then
        write(Uout) rbath
        write(Uout) vrbath
        write(Uout) frbath
      end if

      select case(Ncent)
        case(1)
          write(Uout) rbc11
          write(Uout) vbc11
          write(Uout) fbc11
        case(3)
          write(Uout) rbc31
          write(Uout) vrbc31
          write(Uout) frbc31
      end select
    close(Uout)

    inquire(file=trim(frestart), exist=ex_restart)
    if (ex_restart .and. istep > out_step) then
      call system('mv -f "'//trim(frestart)//'" "'//trim(frestart1)//'"')
    end if
    call system('mv -f "'//trim(frestart_tmp)//'" "'//trim(frestart)//'"')
  end if

return
end subroutine restart_write_bin


!subroutine restart_read_Classical
!  use Parameters
!  implicit none
!  integer :: Uinp
!  integer :: i, j, Inhc
!  real(8) :: pur(3)
!
!  open(newunit=Uinp, file=trim(dir_result)//'/restart.dat', status = 'unknown')
!    read(Uinp,*) Irestep
!    do i = 1, Natom
!      read(Uinp,*) ur(:,i,1)
!    end do
!    do i = 1, Natom
!      read(Uinp,*) pur(:)
!      vur(:,i,1) = pur(:) / sqrt(fictmass(i,1))
!    end do
!    do i = 1, Natom
!      read(Uinp,*) fur(:,i,1)
!    end do
!
!    select case(Ncent)
!      case(1)
!        do Inhc=1,Nnhc
!          read(Uinp,*) rbc11(Inhc),vbc11(Inhc),fbc11(Inhc)
!        enddo
!      case(3)
!        do Inhc=1,Nnhc
!          do i=1,Natom
!            read(Uinp,*) rbc31(:,i,Inhc)
!          enddo
!        enddo
!        do Inhc=1,Nnhc
!          do i=1,Natom
!            read(Uinp,*) vrbc31(:,i,Inhc)
!          enddo
!        enddo
!        do Inhc=1,Nnhc
!          do i=1,Natom
!            read(Uinp,*) frbc31(:,i,Inhc)
!          enddo
!        enddo
!    end select
!  close(Uinp)
!return
!end subroutine restart_read_Classical

!subroutine restart_write_Classical(istep)
!  use Parameters
!  implicit none
!  integer :: Istep, Uout
!  integer :: i, j, Inhc
!
!  if (istep > out_step) call system('cp '//trim(dir_result)//'/restart.dat '//trim(dir_result)//'/restart1.dat')
!
!  open(newunit=Uout, file=trim(dir_result)//'/restart.dat', status = 'unknown')
!    write(Uout,'(i10)') Istep
!    do i = 1, Natom
!      write(Uout,*) ur(:,i,1)
!    end do
!    do i = 1, Natom
!      write(Uout,*) vur(:,i,1) * sqrt(fictmass(i,1))
!    end do
!    do i = 1, Natom
!      write(Uout,*) fur(:,i,1)
!    end do
!
!    select case(Ncent)
!      case(1)
!        do Inhc=1,Nnhc
!          write(Uout,*) rbc11(Inhc),vbc11(Inhc),fbc11(Inhc)
!        enddo
!      case(3)
!        do Inhc=1,Nnhc
!          do i=1,Natom
!            write(Uout,*) rbc31(:,i,Inhc)
!          enddo
!        enddo
!
!        do Inhc=1,Nnhc
!          do i=1,Natom
!            write(Uout,*) vrbc31(:,i,Inhc)
!          enddo
!        enddo
!
!        do Inhc=1,Nnhc
!          do i=1,Natom
!            write(Uout,*) frbc31(:,i,Inhc)
!          enddo
!        enddo
!    end select
!  close(Uout)
!
!return
!end subroutine restart_write_Classical

