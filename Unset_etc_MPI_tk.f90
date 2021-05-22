Subroutine Unset_etc_MPI_tk

  Use MPI
  Use Parameters
  Use Parameter_tk
  Implicit None

!#if 0
!  Deallocate(ireqe)
!  Deallocate(ireqfx)
!  Deallocate(ireqfy)
!  Deallocate(ireqfz)
!  Deallocate(ireqc)
!  Deallocate(ireqdx)
!  Deallocate(ireqdy)
!  Deallocate(ireqdz)
!#endif

  Deallocate (listeach)
  Deallocate (listeachtmp)
  Deallocate (work)

  Deallocate (recvlist)
  Deallocate (recvlist1)
  Deallocate (recvlist2)
!  Deallocate (recvlistx)
!  Deallocate (recvlistx1)
!  Deallocate (recvlistx2)
  Deallocate (ireqa)
  Deallocate (ireqb)


Return
End Subroutine

