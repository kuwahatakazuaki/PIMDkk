subroutine read_parameter
  use,intrinsic :: iso_fortran_env
  use Parameters
  use utility, only: program_abort
  implicit none
  integer :: i, j, Uin, ios
  character(len=100) :: line

  open(newunit=Uin,file=Finp,status='old',iostat=ios)
    if ( ios /= 0) then
      call program_abort('ERROR!!: There is no input file of "input.inp"')
    end if

    InputFile:do
      read(Uin,'(a)',iostat=ios) line
        if ( ios == IOSTAT_END ) then
          print *, "Reaching the bottom of the input file"
          call program_abort('ERROR!!: There is no "$end parameter"')
        end if

      if     (index(line,"$Natom" )         == 1) then; read(Uin,*) Natom
      elseif (index(line,"$Nbead")          == 1) then; read(Uin,*) Nbead
      elseif (index(line,"$Nstep")          == 1) then; read(Uin,*) Nstep
      elseif (index(line,"$temperature")    == 1) then; read(Uin,*) temperature
      elseif (index(line,"$dt")             == 1) then; read(Uin,*) dt
      elseif (index(line,"$Isimulation")    == 1) then; read(Uin,*) Isimulation
      elseif (index(line,"$Nref")           == 1) then; read(Uin,*) Nref
      elseif (index(line,"$Nys")            == 1) then; read(Uin,*) Nys
      elseif (index(line,"$Nnhc")           == 1) then; read(Uin,*) Nnhc
      elseif (index(line,"$Ncent")          == 1) then; read(Uin,*) Ncent
      elseif (index(line,"$gamma")          == 1) then; read(Uin,*) gamma
      elseif (index(line,"$Iforce")         == 1) then; read(Uin,*) Iforce
      elseif (index(line,"$freq1")          == 1) then; read(Uin,*) freq1
      elseif (index(line,"$version")        == 1) then; read(Uin,*) version
      elseif (index(line,"$Lsave_force")    == 1) then; read(Uin,*) Lsave_force
      elseif (index(line,"$Lsave_npa")      == 1) then; read(Uin,*) Lsave_npa
      elseif (index(line,"$Lsave_charge")   == 1) then; read(Uin,*) Lsave_charge
      elseif (index(line,"$Lsave_dipole")   == 1) then; read(Uin,*) Lsave_dipole
      elseif (index(line,"$Lsave_hfcc")     == 1) then; read(Uin,*) Lsave_hfcc
      elseif (index(line,"$Lsave_energy")   == 1) then; read(Uin,*) Lsave_energy
      elseif (index(line,"$Langstrom")      == 1) then; read(Uin,*) Langstrom
      elseif (index(line,"$Lrestart")       == 1) then; read(Uin,*) Lrestart
      elseif (index(line,"$Lrandom_coor")   == 1) then; read(Uin,*) Lrandom_coor
      elseif (index(line,"$seed")           == 1) then
        do i = 1, 4
          read(Uin,*) Iseeds(i)
        end do
      elseif (index(line,"$address_result") == 1) then; read(Uin,'(a)') address
      elseif (index(line,"$address_scr")    == 1) then; read(Uin,'(a)') address2

      elseif (index(line,"$Iumb")           == 1) then; read(Uin,*) Iumb
      elseif (index(line,"$umb_cons")       == 1) then; read(Uin,*) umb_cons
      elseif (index(line,"$umb_atom1")      == 1) then; read(Uin,*) umb_atom1
      elseif (index(line,"$umb_atom2")      == 1) then; read(Uin,*) umb_atom2
      elseif (index(line,"$umb_atom3")      == 1) then; read(Uin,*) umb_atom3

      elseif (index(line,"$end parameter")  == 1) then; exit
      end if
    end do InputFile

    Search_Coords:do
      read(Uin,'(a)',iostat=ios) line
        if ( ios == IOSTAT_END ) exit
        if ( index(line,"$Coords") == 1 ) then
          Natom = 0
          CountNatom:do
            read(Uin,'(a)',iostat=ios) line
            if ( index(line,"$end Coords") == 1 ) then
              exit
            else
              Natom = Natom + 1
            end if
          end do CountNatom
        end if
    end do Search_Coords

    if ( Nbead == 1 ) then
      Isimulation = 10
    end if

    if ( Isimulation == 1 .and. Ncent > 0 ) then
      call program_abort("You should be NVE for RPMD")
    end if

  close(Uin)
contains
  !subroutine sub_read_umbrella
  !  do
  !  end do
  !end subroutine sub_read_umbrella

end subroutine read_parameter

subroutine read_structure
  use,intrinsic :: iso_fortran_env
  use Parameters
  use utility, only: program_abort, atom2mass, lowerchr
  implicit none
  integer :: i, Uin, ios
  character(len=100) :: line

  open(newunit=Uin,file=Finp,status='old',iostat=ios)
    if ( ios /= 0) then
      call program_abort('ERROR!!: There is no input file of "input.inp"')
    end if

    do
      read(Uin,'(a)',iostat=ios) line
        if ( ios == IOSTAT_END ) then
          print *, "Reaching the bottom of the input file"
          call program_abort('ERROR!!: There is no "$Coords"')
        end if

        if ( index(line,"$Coords") == 1 ) then
          do i = 1, Natom
            read(Uin, *) alabel(i), ur(:,i,1)
          end do
          exit
        end if
    end do
  close(Uin)

  do i = 1, Natom
    physmass(i) = atom2mass( alabel(i) )
  end do

  do i = 1, Natom
    if      ( lowerchr(trim(alabel(i))) == "mu" ) then
      alabel(i) = "H "
    else if ( lowerchr(trim(alabel(i))) == "d" ) then
      alabel(i) = "H "
    end if
  end do

  if ( Langstrom .eqv. .True. ) then
    ur(:,:,1) = ur(:,:,1) * AngtoAU
  end if

end subroutine read_structure

