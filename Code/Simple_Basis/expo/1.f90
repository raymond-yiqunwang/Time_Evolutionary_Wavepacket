program main
  implicit none
  integer             :: i,k,n
  integer*8           :: j
  integer             :: contin !should be modified
  complex*8,parameter    :: torr=1d-4
  complex*8,allocatable  :: E(:,:),X(:,:),Xk(:,:),expo1(:,:),expo2(:,:)

  n = 2; contin = 0
  allocate(X(n,n),Xk(n,n),expo1(n,n),expo2(n,n),E(n,n))
  
  E = 0d0
  do i=1,n
    E(i,i) = 1d0
  enddo
  
  X = reshape((/(0d0,1d0),(0d0,2d0),(0d0,3d0),(0d0,4d0)/),(/n,n/))
  expo1 = E
  Xk = E

  i=1d0;j=1d0
  do while( contin == 0 )
    if(i > 15)  exit
    Xk = matmul(Xk,X)
    expo2 = expo1 + 1d0 / j * Xk
    i = i + 1
    j = j * i
    do k=1,n
      write(*,*) expo2(k,:)
    enddo
    write(*,*)
    expo1 = expo2
  enddo
  
end program main
