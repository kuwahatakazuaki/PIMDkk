  Subroutine Force_New

    Use Parameters
  
    Implicit None

    If(NForce==1) Call Force_MOPAC
    If(NForce==2) Call Force_Gaussian
    If(NForce==3) Call Force_DFTB
    If(NForce==4) Call Force_Turbomole

    Return
  End Subroutine

