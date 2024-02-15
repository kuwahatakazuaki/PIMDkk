subroutine exit_program
  use,intrinsic :: iso_fortran_env
  use Parameters
  use utility, only: program_abort
  implicit none
  integer :: i, j, Uin, ios
  character(len=16) :: line
  logical :: Lexit

  if ( MyRank == 0 ) then
    open(newunit=Uin,file=Finp,status='old',iostat=ios)
      if ( ios /= 0) then
        call program_abort('ERROR!!: There is no input file of "input.inp"')
      end if

      InputFile:do
        read(Uin,'(a)',iostat=ios,end=100) line
        if (index(line,"$Lexit" ) == 1) then
          read(Uin,*) Lexit
          if ( Lexit .eqv. .True. ) then
            call program_abort('Exit!!: "exit_program"')
          end if
        end if
      end do InputFile
      100 continue
    close(Uin)
  end if

return
end subroutine exit_program

