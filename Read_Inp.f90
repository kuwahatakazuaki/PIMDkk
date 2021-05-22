Subroutine Read_Inp1

  Use Parameters
  Use Parameter_tk
! Kuwahata 2019/11/20 temp
!    Use MPI
! End Kuwahata 2019/11/20

  Implicit None

  Integer                       :: th
  Character    (Len=90)         :: line

  Character    (Len=90)         :: key1, key2, key3, key4, key5
  Character    (Len=90)         :: key6, key7, key8, key9, key10
!YK Added for increasing input parameters
  Character    (Len=90)         :: key11, key12, key13, key14
  Character    (Len=90)         :: key15, key16, key17, key18
  Character    (Len=90)         :: key19, key20, key21, key22
  Character    (Len=90)         :: key23, key24, key25, key26
  Character    (Len=90)         :: key27, key28, key29, key30
  Character    (Len=90)         :: key31, key32, key33, key34
  Character    (Len=90)         :: key35, key36, key37, key38
  Character    (Len=90)         :: key39, key40, key41, key42
  Character    (Len=90)         :: key43, key44, key45, key46
  Character    (Len=90)         :: key47, key48, key49, key50

!!YK Initial Setting for Some Variables
    NColor=0
    NRestart=0
    Nrstep=0
    gamma=1.0
    angstrom=0
    NSavevel=999999
!    NSavevel=100
    nrandomc=0
    nodipole=1
    nohomo=1
    nocharge=1
    nohfcc=1
    freq1=10.0d0
    npop=0

    !! Setting For Keywords
    key1  = ('$natom')
    key2  = ('$Temperature')
    key3  = ('$nbead')
    key4  = ('$nstep')
    key5  = ('$dt')
!YK added
! Simulation=3, CMD test. Simulation=4, CMD ab initio 
!YK
    key6  = ('$Simulation')
    key7  = ('$nref')
    key8  = ('$End')
    key9  = ('$iys')
    key10 = ('$nnhc')
    key11 = ('$Order')
!YK added adiabaticity parameters to read
    key12 = ('$gamma')
!YK
!YK Type of Ensemble NVE=0, NVT with Nose-Hoover Chain=1
    key13 = ('$nensemble')
!YK 
!YK Type of Electronic Structure Theory Used For Force Calculation
    key14 = ('$theory')
!YK Switch for QM/MM Calculation
    key15 = ('$nqmmm')
!YK
!YK Switch for Additional Basis Sets in Gaussian
    key16 = ('$ngengau')
!YK
!YK Frequency of saving gaussian com files
    key17 = ('$nsavevel')
!YK
!YK Frequency of saving gaussian log files
    key18 = ('$nsavelog')
!YK
!YK Frequency of saving gaussian chk files
    key19 = ('$nsavechk')
!YK
!YK Color of Nose-Hoover Chain on the centroid 
    key20 = ('$ncolor')
!YK
!YK Switch for Restarting from old MD
    key21 = ('$nrestart')
!YK
!YK Specify how to calculate force
    key22 = ('$nforce')
!YK
!YK Specify how to calculate force
! Kuwahata 2019/11/16 Specify the result directory
    key23 = ('$address')
! End Kuwahata 2019/11/16
!YK
!YK Specify random number seed
!YK four odd numbers between 0 and 4095
    key24 = ('$seed')
!YK
!YK Specify Method of Centroid NHC
    key25 = ('$ncent')
!YK Input Coordinates in Angstrom?
    key26 = ('$ang')
!YK Number of Hydrogen to Shoot
    key27 = ('$nshoot')
!YK Shoot Hydrogen at every Ntshoot steps
    key28 = ('$tshoot')
!YK Criterion to delete hydrogen (z-axis)
    key29 = ('$hdel')
!YK Z Coordinates to shoot hydrogen (z-axis)
    key30 = ('$zshoot')
!YK Minimum and Maximum X Coordinates to shoot hydrogen (x-axis)
    key31 = ('$xshoot')
!YK Minimum and Maximum Y Coordinates to shoot hydrogen (y-axis)
    key32 = ('$yshoot')
!YK Energy for shooting hydrogen (eV) 
    key33 = ('$eshoot')
!YK Random bead coordinate?
    key34 = ('$randomc')
!YK Specify the address to copy files to home
! Kuwahata 2019/11/16 Specify the scratch directory
    key35 = ('$address2')
! End Kuwahata 2019/11/16 Specify the scratch directory
!YK Specify if you need Charge data from QC calculation
    key36 = ('$charge')
!YK Specify if you need HOMO-LUMO data from QC calculation
    key37 = ('$homolumo')
!YK Specify if you need dipole data from QC calculatin
    key38 = ('$dipole')
!YK Specify the frequency of the system (used for omega_system)
    key39 = ('$freq')
    key40 = ('$pop')
!tkawatsu MPI process for extra output files
    kproc=0
    key47 = ('$kproc')
    key48 = ('$Coords')
    key49 = ('$hfcc')

    Do
       Read (5,*, Iostat=th) line
       If(th <0) Exit

       If     (Index(Trim(Adjustl(line)), Trim(Adjustl(key1))) /=0)Then
          Read (5, *) natom
       ElseIf (Index(Trim(Adjustl(line)), Trim(Adjustl(key2))) /=0)Then
          Read (5, *) temperature 
       ElseIf (Index(Trim(Adjustl(line)), Trim(Adjustl(key3))) /=0)Then
          Read (5, *) nbead 
       ElseIf (Index(Trim(Adjustl(line)), Trim(Adjustl(key4))) /=0)Then
          Read (5, *) nstep 
       ElseIf (Index(Trim(Adjustl(line)), Trim(Adjustl(key5))) /=0)Then
          Read (5, *) dt
       ElseIf (Index(Trim(Adjustl(line)), Trim(Adjustl(key6))) /=0)Then
          Read (5, *) Simulation 
       ElseIf (Index(Trim(Adjustl(line)), Trim(Adjustl(key7))) /=0)Then
          Read (5, *) nref 
       ElseIf (Index(Trim(Adjustl(line)), Trim(Adjustl(key9))) /=0)Then
          Read (5, *) nys 
       ElseIf (Index(Trim(Adjustl(line)), Trim(Adjustl(key10))) /=0)Then
          Read (5, *) nnhc 
       ElseIf (Index(Trim(Adjustl(line)), Trim(Adjustl(key11))) /=0)Then
          Read (5, *) Order 
!YK added adiabaticity parameters to read
       ElseIf (Index(Trim(Adjustl(line)), Trim(Adjustl(key12))) /=0)Then
          Read (5, *) gamma
!
!YK added to switch the type of ensemble
       ElseIf (Index(Trim(Adjustl(line)), Trim(Adjustl(key13))) /=0)Then
          Read (5, *) Nensemble
!
!YK added to select the electronic structure calculation used for ab initio MD
       ElseIf (Index(Trim(Adjustl(line)), Trim(Adjustl(key14))) /=0)Then
          Read (5, *) theory
!YK added to switch QM/MM calculation
       ElseIf (Index(Trim(Adjustl(line)), Trim(Adjustl(key15))) /=0)Then
          Read (5, *) NQMMM
!
!YK added for additional basis sets (Gen) for gaussian
       ElseIf (Index(Trim(Adjustl(line)), Trim(Adjustl(key16))) /=0)Then
          Read (5, *) NGenGau
!
!YK added for determining the frequency to save gaussian input files
       ElseIf (Index(Trim(Adjustl(line)), Trim(Adjustl(key17))) /=0)Then
          Read (5, *) NSavevel
!
!YK added for determining the frequency to save gaussian output files
       ElseIf (Index(Trim(Adjustl(line)), Trim(Adjustl(key18))) /=0)Then
          Read (5, *) NSavelog
!
!YK added for determining the frequency to save gaussian chk files
       ElseIf (Index(Trim(Adjustl(line)), Trim(Adjustl(key19))) /=0)Then
          Read (5, *) NSavechk
!
!YK added for determining the color of centroid NHC
       ElseIf (Index(Trim(Adjustl(line)), Trim(Adjustl(key20))) /=0)Then
          Read (5, *) NColor
!YK added to switch restart flag from previous MD
       ElseIf (Index(Trim(Adjustl(line)), Trim(Adjustl(key21))) /=0)Then
          Read (5, *) NRestart
!YK added to specify how to calculate force
       ElseIf (Index(Trim(Adjustl(line)), Trim(Adjustl(key22))) /=0)Then
          Read (5, *) NForce
       ElseIf (Index(Trim(Adjustl(line)), Trim(Adjustl(key35))) /=0)Then
          Read (5,'(a)') address2
          address2=trim(address2)
!YK added to specify scratch directory
       ElseIf (Index(Trim(Adjustl(line)), Trim(Adjustl(key23))) /=0)Then
          read(5,'(a)') address
          address=trim(address)
!YK added for specifying random number generator seed
       ElseIf (Index(Trim(Adjustl(line)), Trim(Adjustl(key24))) /=0)Then
          Read (5, *) ISEED1
          Read (5, *) ISEED2
          Read (5, *) ISEED3
          Read (5, *) ISEED4
!
!YK added for determining the color of centroid NHC
       ElseIf (Index(Trim(Adjustl(line)), Trim(Adjustl(key25))) /=0)Then
          Read (5, *) NCent
!
!YK added for determining the unit of the coordinates
       ElseIf (Index(Trim(Adjustl(line)), Trim(Adjustl(key26))) /=0)Then
          Read (5, *) angstrom
!YK added for determining the number of hydrogen to shoot
       ElseIf (Index(Trim(Adjustl(line)), Trim(Adjustl(key27))) /=0)Then
          Read (5, *) nshoot
!YK added for determining the time lag (number of steps) of hydrogen shot
       ElseIf (Index(Trim(Adjustl(line)), Trim(Adjustl(key28))) /=0)Then
          Read (5, *) ntshoot
!YK added for determining the Z coordinate criterion of deleting hydrogen
       ElseIf (Index(Trim(Adjustl(line)), Trim(Adjustl(key29))) /=0)Then
          Read (5, *) hdel
!YK added for determining the Z coordinate to shoot hydrogen
       ElseIf (Index(Trim(Adjustl(line)), Trim(Adjustl(key30))) /=0)Then
          Read (5, *) zshoot
!YK added for determining the X coordinate to shoot hydrogen
       ElseIf (Index(Trim(Adjustl(line)), Trim(Adjustl(key31))) /=0)Then
          Read (5, *) xshootmin,xshootmax
!YK added for determining the Y coordinate to shoot hydrogen
       ElseIf (Index(Trim(Adjustl(line)), Trim(Adjustl(key32))) /=0)Then
          Read (5, *) yshootmin,yshootmax
!YK added for determining the energy of hydrogen shot
       ElseIf (Index(Trim(Adjustl(line)), Trim(Adjustl(key33))) /=0)Then
          Read (5, *) eshoot
!YK added for determining to set random bead coordinate or not
     ElseIf (Index(Trim(Adjustl(line)), Trim(Adjustl(key34))) /=0)Then
        Read (5, *) nrandomc
     ElseIf (Index(Trim(Adjustl(line)), Trim(Adjustl(key36))) /=0)Then
        Read (5, *) nocharge
     ElseIf (Index(Trim(Adjustl(line)), Trim(Adjustl(key37))) /=0)Then
        Read (5, *) nohomo
     ElseIf (Index(Trim(Adjustl(line)), Trim(Adjustl(key38))) /=0)Then
        Read (5, *) nodipole
     ElseIf (Index(Trim(Adjustl(line)), Trim(Adjustl(key39))) /=0)Then
        Read (5, *) freq1
     ElseIf (Index(Trim(Adjustl(line)), Trim(Adjustl(key40))) /=0)Then
        Read (5, *) npop
!     ElseIf (Index(Trim(Adjustl(line)), Trim(Adjustl(key47))) /=0)Then
!        Read (5, *) kproc
!        If(kproc >= nprocs) then
!          print*,"Error: kproc is larger than nprocs. ",kproc," >= ", nprocs
!          Exit
!        Endif
     ElseIf (Index(Trim(Adjustl(line)), Trim(Adjustl(key48))) /=0)Then
        print*,"Error_tk: tk version requires $End before $Coords."
        Stop
     ElseIf (Index(Trim(Adjustl(line)), Trim(Adjustl(key49))) /=0)Then
        Read (5, *) nohfcc
     ElseIf (Index(Trim(Adjustl(line)), Trim(Adjustl('$version_gaussian'))) /=0)Then
        Read (5, *) version_gaussian
     ElseIf (Index(Trim(Adjustl(line)), Trim(Adjustl('$Save_force'))) /=0)Then
        Read (5, *, iostat=ierr) Save_force
          if (ierr .ne. 0) then
            print *, '"Save_force" should be logical'
            print *,  'Type True(T) or False(F)'
            call mpi_abort(mpi_comm_world,-1,ierr)
          end if
     ElseIf (Index(Trim(Adjustl(line)), Trim(Adjustl('$umbrella_sampling'))) /=0)Then
        Read (5, *, iostat=ierr) umbrella_sampling
          if (ierr .ne. 0) then
            print *, '"umbrella_sampling" should be logical'
            print *,  'Type True(T) or False(F)'
            call mpi_abort(mpi_comm_world,-1,ierr)
          end if
     ElseIf (Index(Trim(Adjustl(line)), Trim(Adjustl('$umbrella_width'))) /=0)Then
        Read (5, *) umbrella_width
     ElseIf (Index(Trim(Adjustl(line)), Trim(Adjustl('$umbrella_height'))) /=0)Then
        Read (5, *) umbrella_height
     ElseIf (Index(Trim(Adjustl(line)), Trim(Adjustl('$umbrella_atom1'))) /=0)Then
        Read (5, *) umbrella_atom1
     ElseIf (Index(Trim(Adjustl(line)), Trim(Adjustl('$umbrella_atom2'))) /=0)Then
        Read (5, *) umbrella_atom2
     ElseIf (Index(Trim(Adjustl(line)), Trim(Adjustl(key8))) /=0)Then ! Stop after'$End' keyword
        Exit 
     EndIf
  Enddo

print *, "version_gaussian is ", version_gaussian

  Return
901 print *, 'Type True(T) or False(F)'
End Subroutine

Subroutine Read_Inp2

  Use Parameters
  Implicit None
!  Character    (Len=90)    :: key1, key2
  Character    (Len=90)    :: line
  Integer                  :: th,i

!  key1  = ('$Nuclear_Mass')
!  key2  = ('$Coords')

  Do
    Read (5,*, Iostat=th) line
    If(th <0) Exit

    If     (Index(Trim(Adjustl(line)), Trim(Adjustl('$Coords'))) /=0)Then
      If(NForce==3) Then ! For DFTB
        Do i=1, natom
          Read (5, *) alabel(i), no_atom(i),Physmass(i), ux(i,1), uy(i,1), uz(i,1)
          If(angstrom==1) Then
             ux(i,1)=ux(i,1)*bohr
             uy(i,1)=uy(i,1)*bohr
             uz(i,1)=uz(i,1)*bohr
          EndIf
        Enddo 
      Else
        Do i=1, natom
          Read (5, *) alabel(i), Physmass(i), ux(i,1), uy(i,1), uz(i,1)
        Enddo

        If(angstrom==1) Then
           ux(:,1)=ux(:,1)*bohr
           uy(:,1)=uy(:,1)*bohr
           uz(:,1)=uz(:,1)*bohr
        EndIf

      EndIf
      Exit
    EndIf
!    If     (Index(Trim(Adjustl(line)), Trim(Adjustl(key1))) /=0)Then
!      Do i=1, natom
!        Read (5, *) Physmass(i) 
!      Enddo 
!    ElseIf (Index(Trim(Adjustl(line)), Trim(Adjustl(key2))) /=0)Then
!      Do i=1, natom
!        Read (5, *) ux(i,1), uy(i,1), uz(i,1) 
!      Enddo 
!      Exit 
!    EndIf

  Enddo
  Return
End Subroutine


