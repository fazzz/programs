character(80):: input,output
integer(4)::n,i
real(8)::pi
real(8),allocatable::w(:)

pi=dacos(-1.d0)
call getarg(1,input)
call getarg(2,output)

open(10,file=input)
open(20,file=output)

read(10,*)n

allocate(w(n))
read(10,*)(w(i),i=1,n)

do i=1,n
  w(i) = dsqrt(4.1844*w(i))*1.d3/(2.d0*pi*2.9979d0)
  write(20,'(i5,f10.3)')i,w(i)
  write(*,'(i5,f10.3)')i,w(i)
enddo
end
