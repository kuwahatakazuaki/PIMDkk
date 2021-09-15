Subroutine restart_read_old
  use Parameters
  Implicit none
!
!  /*   read in restart file   */
!

  open(irst, file=trim(address)//'/restart.dat', status = 'unknown')
  read(irst,*) nrstep
  do imode = 1, NBEAD
    do iatom = 1, NATOM
      read(irst,*) ux(iatom,imode), uy(iatom,imode), uz(iatom,imode)
    end do
  end do

  do imode = 1, NBEAD
    do iatom = 1, NATOM
      read(irst,*) vux(iatom,imode), vuy(iatom,imode), vuz(iatom,imode)
    end do
  end do

  do imode = 1, NBEAD
    do iatom = 1, NATOM
      read(irst,*) fux(iatom,imode), fuy(iatom,imode), fuz(iatom,imode)
    end do
  end do

  do imode = 1, NBEAD
    do inhc = 1, NNHC
      do iatom = 1, NATOM
        read(irst,*) xbath(iatom,inhc,imode),ybath(iatom,inhc,imode),zbath(iatom,inhc,imode)
      enddo
    enddo
  enddo

  do imode = 1, NBEAD
    do inhc  = 1, NNHC
      do iatom = 1, NATOM
        read(irst,*) vxbath(iatom,inhc,imode),vybath(iatom,inhc,imode),vzbath(iatom,inhc,imode)
      enddo
    enddo
  enddo

  do imode = 1, NBEAD
    do inhc  = 1, NNHC
      do iatom = 1, NATOM
        read(irst,*) fxbath(iatom,inhc,imode),fybath(iatom,inhc,imode),fzbath(iatom,inhc,imode)
      enddo
    enddo
  enddo

  If(NCent==1) Then

    If(NColor==1) Then

      do inhc=1,nnhc
        read(irst,*) rbc11(inhc),vbc11(inhc),fbc11(inhc)
      enddo

    Else

      do icolor=1,ncolor
        do inhc=1,nnhc
          read(irst,*) rbc1(inhc,icolor),vbc1(inhc,icolor),fbc1(inhc,icolor)
        enddo
      enddo

    EndIf

  EndIf

  If(NCent==3) Then

    If(NColor==1) Then

      do inhc=1,nnhc
        do iatom=1,natom
          read(irst,*) xbc31(iatom,inhc),ybc31(iatom,inhc),zbc31(iatom,inhc)
        enddo
      enddo

      do inhc=1,nnhc
        do iatom=1,natom
          read(irst,*) vxbc31(iatom,inhc),vybc31(iatom,inhc),vzbc31(iatom,inhc)
        enddo
      enddo

      do inhc=1,nnhc
        do iatom=1,natom
          read(irst,*) fxbc31(iatom,inhc),fybc31(iatom,inhc),fzbc31(iatom,inhc)
        enddo
      enddo

    Else

      do icolor=1,ncolor
        do inhc=1,nnhc
          do iatom=1,natom
            read(irst,*) xbc3(iatom,inhc,icolor),ybc3(iatom,inhc,icolor),zbc3(iatom,inhc,icolor)
          enddo
        enddo
      enddo

      do icolor=1,ncolor
        do inhc=1,nnhc
          do iatom=1,natom
            read(irst,*) vxbc3(iatom,inhc,icolor),vybc3(iatom,inhc,icolor),vzbc3(iatom,inhc,icolor)
          enddo
        enddo
      enddo

      do icolor=1,ncolor
        do inhc=1,nnhc
          do iatom=1,natom
            read(irst,*) fxbc3(iatom,inhc,icolor),fybc3(iatom,inhc,icolor),fzbc3(iatom,inhc,icolor)
          enddo
        enddo
      enddo

    end if
  end if

  close(irst)
  return
End Subroutine

Subroutine Restart_Write_old(istep)
  use Parameters
  Implicit none
  Integer :: IStep
!
!        /*   write out restart file   */
!

! Kuwahata 2020/01/23
if (istep > 1) call system('cat '//trim(address)//'/restart.dat >'//trim(address)//'/restart1.dat')
  open(irst, file=trim(address)//'/restart.dat', status = 'unknown')
!  open(irst, file=trim(address)//'/restart1.dat', status = 'unknown')
! End Kuwahata 2020/01/23

  write(irst,'(i10)') IStep
  do imode = 1, NBEAD
    do iatom = 1, NATOM
      write(irst,'(d23.15,1x,d23.15,1x,d23.15)') &
        ux(iatom,imode), uy(iatom,imode), uz(iatom,imode)
    end do
  end do

  do imode = 1, NBEAD
    do iatom = 1, NATOM
      write(irst,'(d23.15,1x,d23.15,1x,d23.15)') &
        vux(iatom,imode), vuy(iatom,imode), vuz(iatom,imode)
    end do
  end do

  do imode = 1, NBEAD
    do iatom = 1, NATOM
      write(irst,'(d23.15,1x,d23.15,1x,d23.15)') &
        fux(iatom,imode), fuy(iatom,imode), fuz(iatom,imode)
    end do
  end do

  do imode = 1, NBEAD
    do inhc = 1, NNHC
      do iatom = 1, NATOM
        write(irst,'(d23.15,1x,d23.15,1x,d23.15)') &
          xbath(iatom,inhc,imode),ybath(iatom,inhc,imode),zbath(iatom,inhc,imode)
      enddo
    enddo
  enddo

  do imode = 1, NBEAD
    do inhc  = 1, NNHC
      do iatom = 1, NATOM
        write(irst,'(d23.15,1x,d23.15,1x,d23.15)') &
          vxbath(iatom,inhc,imode),vybath(iatom,inhc,imode),vzbath(iatom,inhc,imode)
      enddo
    enddo
  enddo

  do imode = 1, NBEAD
    do inhc  = 1, NNHC
      do iatom = 1, NATOM
        write(irst,'(d23.15,1x,d23.15,1x,d23.15)') &
          fxbath(iatom,inhc,imode),fybath(iatom,inhc,imode),fzbath(iatom,inhc,imode)
      enddo
    enddo
  enddo

  If(NCent==1) Then
    If(NColor==1) Then
      do inhc=1,nnhc
        write(irst,'(d23.15,1x,d23.15,1x,d23.15)') rbc11(inhc),vbc11(inhc),fbc11(inhc)
      enddo
    Else
      do icolor=1,ncolor
        do inhc=1,nnhc
          write(irst,'(d23.15,1x,d23.15,1x,d23.15)') rbc1(inhc,icolor),vbc1(inhc,icolor),fbc1(inhc,icolor)
        enddo
      enddo
    EndIf
  EndIf

  If(NCent==3) Then
    If(NColor==1) Then

      do inhc=1,nnhc
        do iatom=1,natom
          write(irst,'(d23.15,1x,d23.15,1x,d23.15)') &
            xbc31(iatom,inhc),ybc31(iatom,inhc),zbc31(iatom,inhc)
        enddo
      enddo

      do inhc=1,nnhc
        do iatom=1,natom
          write(irst,'(d23.15,1x,d23.15,1x,d23.15)') &
            vxbc31(iatom,inhc),vybc31(iatom,inhc),vzbc31(iatom,inhc)
        enddo
      enddo

      do inhc=1,nnhc
        do iatom=1,natom
          write(irst,'(d23.15,1x,d23.15,1x,d23.15)') &
            fxbc31(iatom,inhc),fybc31(iatom,inhc),fzbc31(iatom,inhc)
        enddo
      enddo

    Else

      do icolor=1,ncolor
        do inhc=1,nnhc
          do iatom=1,natom
            write(irst,'(d23.15,1x,d23.15,1x,d23.15)') &
              xbc3(iatom,inhc,icolor),ybc3(iatom,inhc,icolor),zbc3(iatom,inhc,icolor)
          enddo
        enddo
      enddo

      do icolor=1,ncolor
        do inhc=1,nnhc
          do iatom=1,natom
            write(irst,'(d23.15,1x,d23.15,1x,d23.15)') &
              vxbc3(iatom,inhc,icolor),vybc3(iatom,inhc,icolor),vzbc3(iatom,inhc,icolor)
          enddo
        enddo
      enddo

      do icolor=1,ncolor
        do inhc=1,nnhc
          do iatom=1,natom
            write(irst,'(d23.15,1x,d23.15,1x,d23.15)') &
              fxbc3(iatom,inhc,icolor),fybc3(iatom,inhc,icolor),fzbc3(iatom,inhc,icolor)
          enddo
        enddo
      enddo
    end if
  end if
  close(irst)

return
End Subroutine


Subroutine restart_read_Classical
  use Parameters 
  Implicit none
!
!     /*   read in restart file   */
!
  open(irst, file=trim(address)//'/restart.dat', status = 'unknown')
  read(irst,*) nrstep
  do iatom = 1, NATOM
    read(irst,*) ux(iatom,1), uy(iatom,1), uz(iatom,1)
  end do

  do iatom = 1, NATOM
    read(irst,*) vux(iatom,1), vuy(iatom,1), vuz(iatom,1)
  end do

  do iatom = 1, NATOM
    read(irst,*) fux(iatom,1), fuy(iatom,1), fuz(iatom,1)
  end do

  If(NCent==1) Then

    If(NColor==1) Then

      do inhc=1,nnhc
        read(irst,*) rbc11(inhc),vbc11(inhc),fbc11(inhc)
      enddo

    Else

      do icolor=1,ncolor
        do inhc=1,nnhc
          read(irst,*) rbc1(inhc,icolor),vbc1(inhc,icolor),fbc1(inhc,icolor)
        enddo
      enddo

    EndIf

  EndIf

  If(NCent==3) Then

    If(NColor==1) Then

      do inhc=1,nnhc
        do iatom=1,natom
          read(irst,*) xbc31(iatom,inhc),ybc31(iatom,inhc),zbc31(iatom,inhc)
        enddo
      enddo

      do inhc=1,nnhc
        do iatom=1,natom
          read(irst,*) vxbc31(iatom,inhc),vybc31(iatom,inhc),vzbc31(iatom,inhc)
        enddo
      enddo

      do inhc=1,nnhc
        do iatom=1,natom
          read(irst,*) fxbc31(iatom,inhc),fybc31(iatom,inhc),fzbc31(iatom,inhc)
        enddo
      enddo

    Else

      do icolor=1,ncolor
        do inhc=1,nnhc
          do iatom=1,natom
            read(irst,*) xbc3(iatom,inhc,icolor),ybc3(iatom,inhc,icolor),zbc3(iatom,inhc,icolor)
          enddo
        enddo
      enddo

      do icolor=1,ncolor
        do inhc=1,nnhc
          do iatom=1,natom
            read(irst,*) vxbc3(iatom,inhc,icolor),vybc3(iatom,inhc,icolor),vzbc3(iatom,inhc,icolor)
          enddo
        enddo
      enddo

      do icolor=1,ncolor
        do inhc=1,nnhc
          do iatom=1,natom
            read(irst,*) fxbc3(iatom,inhc,icolor),fybc3(iatom,inhc,icolor),fzbc3(iatom,inhc,icolor)
          enddo
        enddo
      enddo

    EndIf
  EndIf

close(irst)
return
End Subroutine

Subroutine Restart_Write_Classical(istep)
  use Parameters
  Implicit none
  Integer :: IStep
!
!        /*   write out restart file   */
!

! Kuwahata 2020/01/23
if (istep > 1) call system('cp '//trim(address)//'/restart.dat '//trim(address)//'/restart1.dat')
  open(irst, file=trim(address)//'/restart.dat', status = 'unknown')
!  open(irst, file=trim(address)//'/restart1.dat', status = 'unknown')
! End Kuwahata 2020/01/23
  write(irst,'(i10)') IStep
  do iatom = 1, NATOM
    write(irst,'(d23.15,1x,d23.15,1x,d23.15)') &
    ux(iatom,1), uy(iatom,1), uz(iatom,1)
  end do

  do iatom = 1, NATOM
    write(irst,'(d23.15,1x,d23.15,1x,d23.15)') &
    vux(iatom,1), vuy(iatom,1), vuz(iatom,1)
  end do

  do iatom = 1, NATOM
    write(irst,'(d23.15,1x,d23.15,1x,d23.15)') &
    fux(iatom,1), fuy(iatom,1), fuz(iatom,1)
  end do

  If(NCent==1) Then
    If(NColor==1) Then
      do inhc=1,nnhc
        write(irst,'(d23.15,1x,d23.15,1x,d23.15)') rbc11(inhc),vbc11(inhc),fbc11(inhc)
      enddo
    Else
      do icolor=1,ncolor
        do inhc=1,nnhc
          write(irst,'(d23.15,1x,d23.15,1x,d23.15)') rbc1(inhc,icolor),vbc1(inhc,icolor),fbc1(inhc,icolor)
        enddo
      enddo
    EndIf
  EndIf

  If(NCent==3) Then
    If(NColor==1) Then

      do inhc=1,nnhc
        do iatom=1,natom
          write(irst,'(d23.15,1x,d23.15,1x,d23.15)') &
          xbc31(iatom,inhc),ybc31(iatom,inhc),zbc31(iatom,inhc)
        enddo
      enddo

      do inhc=1,nnhc
        do iatom=1,natom
          write(irst,'(d23.15,1x,d23.15,1x,d23.15)') &
          vxbc31(iatom,inhc),vybc31(iatom,inhc),vzbc31(iatom,inhc)
        enddo
      enddo

      do inhc=1,nnhc
        do iatom=1,natom
          write(irst,'(d23.15,1x,d23.15,1x,d23.15)') &
          fxbc31(iatom,inhc),fybc31(iatom,inhc),fzbc31(iatom,inhc)
        enddo
      enddo

    Else

      do icolor=1,ncolor
        do inhc=1,nnhc
          do iatom=1,natom
            write(irst,'(d23.15,1x,d23.15,1x,d23.15)') &
            xbc3(iatom,inhc,icolor),ybc3(iatom,inhc,icolor),zbc3(iatom,inhc,icolor)
          enddo
        enddo
      enddo

      do icolor=1,ncolor
        do inhc=1,nnhc
          do iatom=1,natom
            write(irst,'(d23.15,1x,d23.15,1x,d23.15)') &
            vxbc3(iatom,inhc,icolor),vybc3(iatom,inhc,icolor),vzbc3(iatom,inhc,icolor)
          enddo
        enddo
      enddo

      do icolor=1,ncolor
        do inhc=1,nnhc
          do iatom=1,natom
            write(irst,'(d23.15,1x,d23.15,1x,d23.15)') &
            fxbc3(iatom,inhc,icolor),fybc3(iatom,inhc,icolor),fzbc3(iatom,inhc,icolor)
          enddo
        enddo
      enddo

    EndIf
  EndIf

  close(irst)

! Kuwahata 2020/01/23
!  call system('cp '//trim(address)//'/restart1.dat '//trim(address)//'/restart.dat')
! End Kuwahata 2020/01/23
return
End Subroutine

    Subroutine restart_read_Shoot
      use Parameters 
      Implicit none
!
!     /*   read in restart file   */
!
      open(irst, file=trim(address)//'/restart.dat', status = 'unknown')
      read(irst,*) nrstep
      read(irst,*) iseed1
      read(irst,*) iseed2
      read(irst,*) iseed3
      read(irst,*) iseed4
      read(irst,*) natom
      read(irst,*) nhydtot
      read(irst,*) nmshot
      read(irst,*) nhydrfl
      read(irst,*) nhyddel
      read(irst,*) nhydthr
      read(irst,*) gnkt
      do iatom=1,nshoot
         read(irst,*) index_hyd(iatom),ncontact(iatom),nhspin(iatom)
      enddo
      do iatom=1,nshoot
         read(irst,*) &
         hvelin(iatom),hvelout(iatom),hinitx(iatom),hinity(iatom)
      enddo
      do imode = 1, NBEAD
         do iatom = 1, NATOM0+NSHOOT
            read(irst,*) ux(iatom,imode), uy(iatom,imode), uz(iatom,imode)
         end do
      end do
!
      do imode = 1, NBEAD
         do iatom = 1, NATOM0+NSHOOT
            read(irst,*) vux(iatom,imode), vuy(iatom,imode), vuz(iatom,imode)
         end do
      end do
!
      do imode = 1, NBEAD
         do iatom = 1, NATOM0+NSHOOT
            read(irst,*) fux(iatom,imode), fuy(iatom,imode), fuz(iatom,imode)
         end do
      end do
!
      do imode = 1, NBEAD
         do inhc = 1, NNHC
            do iatom = 1, NATOM0+NSHOOT
               read(irst,*) xbath(iatom,inhc,imode),ybath(iatom,inhc,imode),zbath(iatom,inhc,imode)
            enddo
         enddo
      enddo
!
      do imode = 1, NBEAD
         do inhc  = 1, NNHC
            do iatom = 1, NATOM0+NSHOOT
               read(irst,*) vxbath(iatom,inhc,imode),vybath(iatom,inhc,imode),vzbath(iatom,inhc,imode)
            enddo
         enddo
      enddo
!
      do imode = 1, NBEAD
         do inhc  = 1, NNHC
            do iatom = 1, NATOM0+NSHOOT
               read(irst,*) fxbath(iatom,inhc,imode),fybath(iatom,inhc,imode),fzbath(iatom,inhc,imode)
            enddo
         enddo
      enddo
!


      close(irst)
      return
    End Subroutine

    Subroutine Restart_Write_Shoot(istep)
      use Parameters 
      Implicit none
      Integer :: IStep
!
!        /*   write out restart file   */
!
      open(irst, file=trim(address)//'/restart1.dat', status = 'unknown')
      write(irst,'(i10)') IStep
      write(irst,'(I10)') iseed1
      write(irst,'(I10)') iseed2
      write(irst,'(I10)') iseed3
      write(irst,'(I10)') iseed4
      write(irst,'(i10)') natom
      write(irst,'(i10)') nhydtot
      write(irst,'(i10)') nmshot
      write(irst,'(i10)') nhydrfl
      write(irst,'(i10)') nhyddel
      write(irst,'(i10)') nhydthr
      write(irst,'(d23.15)') gnkt
      do iatom=1,nshoot
         write(irst,'(i10,1x,i10,1x,i10)') &
         index_hyd(iatom),ncontact(iatom),nhspin(iatom)
      enddo
      do iatom=1,nshoot
         write(irst,'(d23.15,1x,d23.15,1x,d23.15,1x,d23.15)') &
         hvelin(iatom),hvelout(iatom),hinitx(iatom),hinity(iatom)
      enddo
      do imode = 1, NBEAD
         do iatom = 1, NATOM0+NSHOOT
            write(irst,'(d23.15,1x,d23.15,1x,d23.15)') &
            ux(iatom,imode), uy(iatom,imode), uz(iatom,imode)
         end do
      end do
!
      do imode = 1, NBEAD
         do iatom = 1, NATOM0+NSHOOT
            write(irst,'(d23.15,1x,d23.15,1x,d23.15)') &
            vux(iatom,imode), vuy(iatom,imode), vuz(iatom,imode)
         end do
      end do
!
      do imode = 1, NBEAD
         do iatom = 1, NATOM0+NSHOOT
            write(irst,'(d23.15,1x,d23.15,1x,d23.15)') &
            fux(iatom,imode), fuy(iatom,imode), fuz(iatom,imode)
         end do
      end do
!
      do imode = 1, NBEAD
         do inhc = 1, NNHC
            do iatom = 1, NATOM0+NSHOOT
               write(irst,'(d23.15,1x,d23.15,1x,d23.15)') &
               xbath(iatom,inhc,imode),ybath(iatom,inhc,imode),zbath(iatom,inhc,imode)
            enddo
         enddo
      enddo
!
      do imode = 1, NBEAD
         do inhc  = 1, NNHC
            do iatom = 1, NATOM0+NSHOOT
               write(irst,'(d23.15,1x,d23.15,1x,d23.15)') &
               vxbath(iatom,inhc,imode),vybath(iatom,inhc,imode),vzbath(iatom,inhc,imode)
            enddo
         enddo
      enddo
!
      do imode = 1, NBEAD
         do inhc  = 1, NNHC
            do iatom = 1, NATOM0+NSHOOT
               write(irst,'(d23.15,1x,d23.15,1x,d23.15)') &
               fxbath(iatom,inhc,imode),fybath(iatom,inhc,imode),fzbath(iatom,inhc,imode)
            enddo
         enddo
      enddo
!

!
      close(irst)

      call system('cp '//trim(address)//'/restart1.dat '//trim(address)//'/restart.dat')

      return
    End Subroutine

