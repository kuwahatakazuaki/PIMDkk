subroutine restart_read
  use Parameters
  implicit none
  integer :: Uinp
  integer :: i, j, inhc, icolor
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
        !read(Uinp,*) vur(:,i,j)
      end do
    end do

    do j = 1, Nbead
      do i = 1, Natom
        read(Uinp,*) fur(:,i,j)
      end do
    end do

    do j = 1, Nbead
      do inhc = 1, Nnhc
        do i = 1, Natom
          read(Uinp,*) rbath(:,i,inhc,j)
        enddo
      enddo
    enddo

    do j = 1, Nbead
      do inhc  = 1, Nnhc
        do i = 1, Natom
          read(Uinp,*) vrbath(:,i,inhc,j)
        enddo
      enddo
    enddo

    do j = 1, Nbead
      do inhc  = 1, Nnhc
        do i = 1, Natom
          read(Uinp,*) frbath(:,i,inhc,j)
        enddo
      enddo
    enddo

    select case(Ncent)
      case(0)
        continue
      case(1)
        do inhc=1,nnhc
          read(Uinp,*) rbc11(inhc),vbc11(inhc),fbc11(inhc)
        enddo
      case(3)
        do inhc=1,nnhc
          do i=1,natom
            read(Uinp,*) rbc31(:,i,inhc)
          enddo
        enddo
        do inhc=1,nnhc
          do i=1,natom
            read(Uinp,*) vrbc31(:,i,inhc)
          enddo
        enddo
        do inhc=1,nnhc
          do i=1,natom
            read(Uinp,*) frbc31(:,i,inhc)
          enddo
        enddo
    end select
  close(Uinp)
end if
  return
end subroutine restart_read

subroutine Restart_Write(istep)
  use Parameters
  implicit none
  integer :: Istep, Uout
  integer :: i, j, inhc, icolor

if ( MyRank == 0 ) then
  if (istep > out_step) then
    call system('cat '//trim(dir_result)//'/restart.dat >'//trim(dir_result)//'/restart1.dat')
  end if

  open(newunit=Uout, file=trim(dir_result)//'/restart.dat', status = 'unknown')
    write(Uout,'(i10)') Istep
    do j = 1, Nbead
      do i = 1, Natom
        write(Uout,*) ur(:,i,j)
      end do
    end do

    do j = 1, Nbead
      do i = 1, Natom
        write(Uout,*) vur(:,i,j) * sqrt(fictmass(i,j))
        !write(Uout,*) vur(:,i,j)
      end do
    end do

    do j = 1, Nbead
      do i = 1, Natom
        write(Uout,*) fur(:,i,j)
      end do
    end do

    do j = 1, Nbead
      do inhc = 1, Nnhc
        do i = 1, Natom
          write(Uout,*) rbath(:,i,inhc,j)
        enddo
      enddo
    enddo
    do j = 1, Nbead
      do inhc  = 1, Nnhc
        do i = 1, Natom
          write(Uout,*) vrbath(:,i,inhc,j)
        enddo
      enddo
    enddo
    do j = 1, Nbead
      do inhc  = 1, Nnhc
        do i = 1, Natom
          write(Uout,*) frbath(:,i,inhc,j)
        enddo
      enddo
    enddo

    select case(Ncent)
      case(1)
        do inhc=1,nnhc
          write(Uout,*) rbc11(inhc),vbc11(inhc),fbc11(inhc)
        enddo
      case(3)
        do inhc=1,nnhc
          do i=1,natom
            write(Uout,*) rbc31(:,i,inhc)
          enddo
        enddo
        do inhc=1,nnhc
          do i=1,natom
            write(Uout,*) vrbc31(:,i,inhc)
          enddo
        enddo
        do inhc=1,nnhc
          do i=1,natom
            write(Uout,*) frbc31(:,i,inhc)
          enddo
        enddo
    end select
  close(Uout)
end if

return
end subroutine Restart_Write


subroutine restart_read_Classical
  use Parameters
  implicit none
  integer :: Uinp
  integer :: i, j, inhc, icolor

  open(newunit=Uinp, file=trim(dir_result)//'/restart.dat', status = 'unknown')
    read(Uinp,*) Irestep
    do i = 1, Natom
      read(Uinp,*) ur(:,i,1)
    end do
    do i = 1, Natom
      read(Uinp,*) vur(:,i,1)
    end do
    do i = 1, Natom
      read(Uinp,*) fur(:,i,1)
    end do

    select case(Ncent)
      case(1)
        do inhc=1,nnhc
          read(Uinp,*) rbc11(inhc),vbc11(inhc),fbc11(inhc)
        enddo
      case(3)
        do inhc=1,nnhc
          do i=1,natom
            read(Uinp,*) rbc31(:,i,inhc)
          enddo
        enddo
        do inhc=1,nnhc
          do i=1,natom
            read(Uinp,*) vrbc31(:,i,inhc)
          enddo
        enddo
        do inhc=1,nnhc
          do i=1,natom
            read(Uinp,*) frbc31(:,i,inhc)
          enddo
        enddo
    end select
  close(Uinp)
return
end subroutine restart_read_Classical

subroutine Restart_Write_Classical(istep)
  use Parameters
  implicit none
  integer :: Istep, Uout
  integer :: i, j, inhc, icolor

  if (istep > 1) call system('cp '//trim(dir_result)//'/restart.dat '//trim(dir_result)//'/restart1.dat')

  open(Uout, file=trim(dir_result)//'/restart.dat', status = 'unknown')
    write(Uout,'(i10)') IStep
    do i = 1, Natom
      write(Uout,*) ur(:,i,1)
    end do
    do i = 1, Natom
      write(Uout,*) vur(:,i,1)
    end do
    do i = 1, Natom
      write(Uout,*) fur(:,i,1)
    end do

    select case(Ncent)
      case(1)
        do inhc=1,nnhc
          write(Uout,*) rbc11(inhc),vbc11(inhc),fbc11(inhc)
        enddo
      case(3)
        do inhc=1,nnhc
          do i=1,natom
            write(Uout,*) rbc31(:,i,inhc)
          enddo
        enddo

        do inhc=1,nnhc
          do i=1,natom
            write(Uout,*) vrbc31(:,i,inhc)
          enddo
        enddo

        do inhc=1,nnhc
          do i=1,natom
            write(Uout,*) frbc31(:,i,inhc)
          enddo
        enddo
    end select
  close(Uout)

return
end subroutine Restart_Write_Classical

