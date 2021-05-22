Subroutine Close_files_MPI_tk

  Use MPI
  Use Parameters
  Use Parameter_tk
  Implicit None

  if(myrank==0) then
     Close(igetx)
     Close(igetf)
     Close(igete)
     Close(igetxyz)
!     Close(igethl)

  Write(*,'(a)') &
    &'***********************************************************************&
    &************************************************************************'

!  endif
!  if(myrank==kproc) then

    if(nodipole==0) then
      Close(igetd)
    endif
    if(nocharge==0) then
      Close(igetc)
    endif
  endif
Return
End Subroutine

