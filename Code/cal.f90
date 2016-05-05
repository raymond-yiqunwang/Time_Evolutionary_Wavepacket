program cal
  implicit none
  integer     :: i
  complex*16  :: temp(2,2)
  complex*16  :: result(2)
  complex*16  :: H(2,2),Ope(2,2),psi_0(2)
  complex*16  :: ci,phi_0(2,2),phi_1(2,2),phi_2(2,2)

  H = reshape((/1,0,0,4/),(/2,2/))
  ci = (0d0,1d0)
  psi_0 = (/1,2/)
  !phi_0 = 0.5*reshape((/1,0,0,1/),(/2,2/))
  !phi_1 = reshape((/-1,0,0,1/),(/2,2/))
  !phi_2 = reshape((/1,0,0,1/),(/2,2/))
  !result = 2d0*BESSEL_JN(0,1.5*0.01)*phi_0 + 2d0*BESSEL_JN(1,1.5*0.01)*phi_1*(-ci) + 2d0*BESSEL_JN(2,1.5*0.01)*(-1)*phi_2 
  !result = exp(-ci*2.5*0.01)*result
  temp = expo(-ci*H*0.01)
  result =  MATMUL(temp,psi_0)
  print *, 'result:',result

contains
function expo(x)
  implicit none
  integer       :: i,k
  complex*16    :: expo(2,2),e(2,2)
  complex*16    :: H_temp(2,2)
  complex*16    :: x(2,2)
  
  e=0d0
  do i=1,2
    e(i,i)=1
  enddo

  expo = e
  H_temp = x
  k=1
  do i=1,30
    if(i>1)then
      H_temp = MATMUL(x,H_temp)
      k = k*i
    endif
    expo = expo + H_temp/k
  enddo
  return
end function expo
  
end program cal
