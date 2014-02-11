subroutine CalcstuffBFGS(X,ndvart,fobj,dfdD,fct)
  use dimpce,only:fctindx
  implicit none

  integer  :: ndvart,fct
  double precision :: X(ndvart),fobj,dfdD(ndvart),x3
  double precision ::  rho, L, sigmay, pi, p, E, Fs  

  fctindx=fct

  call get_f(ndvart,12,x,fobj)

  call get_df(ndvart,12,x,dfDD)


  return
end subroutine CalcstuffBFGS

!++++++++++++++++++++++++++++++++++
subroutine CalcExact(X,ndvart,fobj,dfdD,fct,DATA)
  use dimpce,only:fctindx,DAT,mainprog
  implicit none

  real*8::DATA(20)
  integer  :: ndvart,fct
  double precision :: X(ndvart),fobj,dfdD(ndvart),x3
  double precision ::  rho, L, sigmay, pi, p, E, Fs  

  fctindx=fct
  mainprog=.false.
  DAT(1:20)=DATA(1:20)
  
  call get_f(ndvart,12,x,fobj)

  call get_df(ndvart,12,x,dfDD)


  return
end subroutine CalcExact
