program main
  implicit none
  integer       :: i,k
  integer       :: A(2,2),B(1,2),C(2,1),D(2,1)

  A = reshape((/1,3,2,4/),(/2,2/))
  B = reshape((/1,2/),(/1,2/))
  do i=1,2
    write(*,*) A(i,:)
  enddo
  C = transpose(B)
  D = MATMUL(A,C)
write(*,*) 'matmul:',D
end program main
