! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
! This subroutine is Takahashi-Imada 4th-order correction. 
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 

  Subroutine Takahashi_Imada

    Use Parameters
    Implicit None

    Double Precision                   :: fact_ti
    Integer                            :: kx, ky, kz
    Integer                            :: jx, jy, jz
    Integer                            :: i, j, k 


        pot_ti (1:nbead)   = 0.d0
        fact_ti      = beta*beta/24.d0/dble(nbead)/dble(nbead)
    Do i=1, nbead
      Do j=1, natom
        pot_ti(i) = pot_ti(i)                                           &
                  + (-fx(j,i))*(-fx(j,i))/physmass(j)*fact_ti           &
                  + (-fy(j,i))*(-fy(j,i))/physmass(j)*fact_ti           &
                  + (-fz(j,i))*(-fz(j,i))/physmass(j)*fact_ti   
      Enddo
!       Write(*,'(f21.14)') pot_ti(i)
    Enddo

        fx_ti (1:natom,1:nbead)  = 0.d0
        fy_ti (1:natom,1:nbead)  = 0.d0
        fz_ti (1:natom,1:nbead)  = 0.d0
        fact_ti      = beta*beta/12.d0/dble(nbead)/dble(nbead)


    Do i=1, nbead
      Do k=1, natom
        kx = 3*(k-1) + 1
        ky = 3*(k-1) + 2
        kz = 3*(k-1) + 3
        Do j=1, natom
          jx = 3*(j-1) + 1
          jy = 3*(j-1) + 2
          jz = 3*(j-1) + 3
          fx_ti(k,i)  = fx_ti(k,i)                                                 & 
                      + hess(kx, jx, i) * (-fx(j,i)) / physmass(j) * fact_ti       &
                      + hess(kx, jy, i) * (-fy(j,i)) / physmass(j) * fact_ti       &
                      + hess(kx, jz, i) * (-fz(j,i)) / physmass(j) * fact_ti        
          fy_ti(k,i)  = fy_ti(k,i)                                                 & 
                      + hess(ky, jx, i) * (-fx(j,i)) / physmass(j) * fact_ti       &
                      + hess(ky, jy, i) * (-fy(j,i)) / physmass(j) * fact_ti       &
                      + hess(ky, jz, i) * (-fz(j,i)) / physmass(j) * fact_ti        
          fz_ti(k,i)  = fz_ti(k,i)                                                 & 
                      + hess(kz, jx, i) * (-fx(j,i)) / physmass(j) * fact_ti       &
                      + hess(kz, jy, i) * (-fy(j,i)) / physmass(j) * fact_ti       &
                      + hess(kz, jz, i) * (-fz(j,i)) / physmass(j) * fact_ti        
        Enddo
!       Write(*,"(3(f14.7,2x))") fx_ti(k,i), fy_ti(k,i), fz_ti(k,i)
      Enddo
    Enddo

  End Subroutine
