  Subroutine Force_Gaussian_MPI

!$  use omp_lib
    Use Parameters
    Use MPI
    Implicit None

    Character    (Len=90)              :: key1, key2, key3, key4, key5, line
    Character    (Len=5)               :: nobeadtmp
    Integer                            :: iline, numeach, ista, iend, itemp,id,imode2,iatom2,nhmod,icount
    Logical                            :: lopen
    Double Precision                   :: enetemp
    Double Precision,dimension(natom)  :: x1,y1,z1
    Character(len=256) :: g0xtmp
    ! koba add 2018/6/4
    integer :: chl, chr
    key1  = ('SCF Done')
    key2  = ('Number     Number              X              Y              Z')
    key3  = ('Dipole        =')
    key4  = ('Mulliken atomic charges')

!    Call system('echo $OMP_NUM_THREADS')
!    Call system('echo $GAUSS_SCRDIR')

    id=0
    nhmod=mod(nbead,nprocs)
    numeach=(nbead-nhmod)/nprocs
    if(nhmod/=0) then
       if(myrank<nhmod) then
          ista=myrank*numeach+1+myrank
          iend=ista+numeach
       else
          ista=myrank*numeach+nhmod+1
          iend=ista+numeach-1
       endif
    else
       ista=myrank*numeach+1
       iend=ista+numeach-1
    endif

    Call MPI_BARRIER(MPI_COMM_WORLD,IERR)

!$omp parallel private(nobeadtmp,line,iatom2,id,iline,enetemp)
!$omp do
    do imode2=ista,iend
!$     id=omp_get_thread_num()
 200   continue
       write(nobeadtmp,'(i5.5)') imode2
       call system('mkdir -p '//trim(address)//'/'//nobeadtmp)
       call system('cp gauss.tmp '//trim(address)//'/'//nobeadtmp)
       If(NGenGau==1) Then
          call system('cp gauss.bss '//trim(address)//'/'//nobeadtmp)
       EndIf
       open(igauss+id,file=trim(address)//'/'//nobeadtmp//'/gauss.tmp1',status='unknown')
       write(igauss+id,'(a)') '%Chk='//trim(address)//'/'//nobeadtmp//'/gauss.chk'
!tkawatsu
       write(igauss+id,'(a)') '%RWF='//trim(address)//'/'//nobeadtmp//'/'
       write(igauss+id,'(a)') '%Int='//trim(address)//'/'//nobeadtmp//'/'
       write(igauss+id,'(a)') '%D2E='//trim(address)//'/'//nobeadtmp//'/'
!tkawatsu end
       close(igauss+id)
       open(igauss+id,file=trim(address)//'/'//nobeadtmp//'/gauss.xyz',status='unknown')
       do iatom2=1,natom
          write(igauss+id,9999) alabel(iatom2),x(iatom2,imode2)/bohr, &
                             y(iatom2,imode2)/bohr,z(iatom2,imode2)/bohr
       enddo
       write(igauss+id,*) 
       close(igauss+id)
       call system('cat '//trim(address)//'/'//nobeadtmp//'/gauss.tmp1 > '//trim(address)//'/'//nobeadtmp//'/gauss.com')
       call system('cat '//trim(address)//'/'//nobeadtmp//'/gauss.tmp >> '//trim(address)//'/'//nobeadtmp//'/gauss.com')
       call system('cat '//trim(address)//'/'//nobeadtmp//'/gauss.xyz >> '//trim(address)//'/'//nobeadtmp//'/gauss.com')

       If(NGenGau==1) Then
          call system('cat '//trim(address)//'/'//nobeadtmp//'/gauss.bss >> '//trim(address)//'/'//nobeadtmp//'/gauss.com')
       EndIf

!tkawatsu
!       call system('g09 < '//trim(address)//'/'//nobeadtmp//'/gauss.com > '//trim(address)//'/'//nobeadtmp//'/gauss.log')
!       g0xtmp='g03run_p '//trim(address)//'/'//nobeadtmp//'/gauss.com '//trim(address)//'/'//nobeadtmp//'/gauss.log'
       g0xtmp=trim(address)//'/g0xrun_p '//trim(address)//'/'//nobeadtmp//'/gauss.com'
!       print*,trim(g0xtmp)//' '//trim(address)//'/'//nobeadtmp//'/gauss.log '//trim(address)//'/'//nobeadtmp 
       call system(trim(g0xtmp)//' '//trim(address)//'/'//nobeadtmp//'/gauss.log '//trim(address)//'/'//nobeadtmp)
       open(igauss+id,file=trim(address)//'/'//nobeadtmp//'/gauss.log')

 10    continue
       Read(igauss+id,'(a)') line

!       iline=index(trim(line),trim(key1))
       iline=index(line(2:9),trim(key1))
       if(iline==0) goto 10
!       read(line(25:40),*) enetemp
       !read(line(21:40),*) enetemp
       chl = index(line, '=')
       chr = index(line(chl+1:len(line)-chl), '.') + chl
       chr = index(line(chr+1:len(line)-chr), ' ') + chr
       read(line(chl+1:chr-1),*) enetemp
       Eenergy(imode2)=enetemp
       Rewind(igauss+id)

 20    continue
       Read(igauss+id,'(a)') line
       iline=index(trim(line),trim(key2))
       if(iline==0) goto 20
       Read(igauss+id,*) 
       do iatom2=1,natom
          read(igauss+id,'(a)') line
          read(line(24:38),*) fx(iatom2,imode2)
          read(line(39:53),*) fy(iatom2,imode2)
          read(line(54:68),*) fz(iatom2,imode2)
       enddo
       Rewind(igauss+id)
 30    continue
       if(nodipole==0) then
          Read(igauss+id,'(a)') line
          iline=index(trim(line),trim(key3))
          if(iline==0) goto 30
          Read(line(17:31),*) dipolex(imode2)
          Read(line(32:46),*) dipoley(imode2)
          Read(line(47:61),*) dipolez(imode2)
          Rewind(igauss+id)
       endif
 40    continue
       if(nocharge==0) then
          Read(igauss+id,'(a)') line
          iline=index(trim(line),trim(key4))
          if(iline==0) goto 40
          Read(igauss+id,*) 
          do iatom2=1,natom
             Read(igauss+id,'(a)') line
             Read(line(11:21),*) charge(iatom2,imode2)
          enddo
       endif
       close(igauss+id)
       call system('rm -rf '//trim(address)//'/'//nobeadtmp)
    enddo
!$omp end do
!$omp end parallel
  call system('rm -rf '//trim(address)//'/'//nobeadtmp//' /Gau*')
  nhmod=mod(nbead,nprocs)
  numeach=(nbead-nhmod)/nprocs
  do imode=1,nbead
     if(nhmod/=0) then
        if(numeach*nhmod+nhmod>=imode) then
           itemp=int(dble(imode-1)/dble(numeach+1))
        else
           itemp=int(dble(imode-1-numeach*nhmod-nhmod)/dble(numeach))
           itemp=itemp+nhmod
        endif
     else
        itemp=int(dble(imode-1)/dble(numeach))
     endif
     call MPI_BCAST(Eenergy(imode),1,MPI_DOUBLE_PRECISION,itemp,MPI_COMM_WORLD,IERR)
     do iatom=1,natom
        x1(iatom)=fx(iatom,imode)
        y1(iatom)=fy(iatom,imode)
        z1(iatom)=fz(iatom,imode)
     enddo
     call MPI_BCAST(x1,natom,MPI_DOUBLE_PRECISION,itemp,MPI_COMM_WORLD,IERR)
     call MPI_BCAST(y1,natom,MPI_DOUBLE_PRECISION,itemp,MPI_COMM_WORLD,IERR)
     call MPI_BCAST(z1,natom,MPI_DOUBLE_PRECISION,itemp,MPI_COMM_WORLD,IERR)
     do iatom=1,natom
        fx(iatom,imode)=x1(iatom)
        fy(iatom,imode)=y1(iatom)
        fz(iatom,imode)=z1(iatom)
     enddo
     if(nodipole==0) then
        call MPI_BCAST(dipolex(imode),1,MPI_DOUBLE_PRECISION,itemp,MPI_COMM_WORLD,IERR)
        call MPI_BCAST(dipoley(imode),1,MPI_DOUBLE_PRECISION,itemp,MPI_COMM_WORLD,IERR)
        call MPI_BCAST(dipolez(imode),1,MPI_DOUBLE_PRECISION,itemp,MPI_COMM_WORLD,IERR)
     endif
     if(nocharge==0) then
        do iatom=1,natom
           x1(iatom)=charge(iatom,imode)
        enddo
        call MPI_BCAST(x1,natom,MPI_DOUBLE_PRECISION,itemp,MPI_COMM_WORLD,IERR)
        do iatom=1,natom
           charge(iatom,imode)=x1(iatom)
        enddo
     endif
  enddo
  Call MPI_BARRIER(MPI_COMM_WORLD,IERR)

  If(MyRank==0) Then
     Open(igetx,file=trim(address)//'/coor.dat',status='unknown',form='formatted',position='append')
     Open(igetf,file=trim(address)//'/force.dat',status='unknown',form='formatted',position='append')
     if(nodipole==0) then
        Open(igetd,file=trim(address)//'/dipole.dat',status='unknown',form='formatted',position='append')
        DO imode=1,nbead
           write(igetd,9998) dipolex(imode),dipoley(imode),dipolez(imode)
        ENDDO
        Close(igetd)
     endif
     Open(igete,file=trim(address)//'/ene.dat',status='unknown',form='formatted',position='append')
     if(nocharge==0) then
        Open(igetc,file=trim(address)//'/charge.dat',status='unknown',form='formatted',position='append')
        DO imode=1,nbead
           DO iatom=1,natom
              write(igetc,9996) charge(iatom,imode)
           ENDDO
        ENDDO
        Close(igetc)
     endif
     Open(igetxyz,file=trim(address)//'/cent.xyz',status='unknown',form='formatted',position='append')
!    Open(igethl,file=trim(address)//'/homolumo.dat',status='unknown',form='formatted',position='append')


     write(igetxyz,'(I5)') natom
     write(igetxyz,'(I10)') istepsv 
     DO iatom=1,natom
        write(igetxyz,9999) alabel(iatom),ux(iatom,1)/bohr,uy(iatom,1)/bohr,uz(iatom,1)/bohr
     ENDDO
     write(igetx,'(I5)') natom*nbead
     write(igetx,'(I10)') istepsv 
     DO imode=1,nbead
        DO iatom=1,natom
           write(igetx,9999) alabel(iatom),x(iatom,imode)/bohr,y(iatom,imode)/bohr,z(iatom,imode)/bohr
        ENDDO
     ENDDO
     DO imode=1,nbead
        DO iatom=1,natom
           write(igetf,9998) fx(iatom,imode),fy(iatom,imode),fz(iatom,imode)
        ENDDO
     ENDDO
     DO imode=1,nbead
        write(igete,9996) Eenergy(imode)
     ENDDO
!      DO imode=1,nbead
!         write(igethl,9997) homo(imode),lumo(imode)
!      ENDDO
     Close(igetx)
     Close(igetf)
     Close(igete)
     Close(igetxyz)
!    Close(igethl)
  EndIf

  do imode=1,nbead
     do iatom=1,natom
        fx(iatom,imode)=fx(iatom,imode)/dp
        fy(iatom,imode)=fy(iatom,imode)/dp
        fz(iatom,imode)=fz(iatom,imode)/dp
     enddo
  enddo
  potential=0.0D+00
  do imode=1,nbead
     potential=potential+Eenergy(imode)
  enddo
  potential=potential/dp

9999 format(a2,1x,E15.9,1x,E15.9,1x,E15.9) 
9998 format(3E23.15) 
9997 format(2E23.15) 
9996 format(E23.15) 
9995 format(4E23.15) 

Return
End Subroutine
