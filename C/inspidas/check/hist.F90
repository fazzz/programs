integer(4)::n
integer(4)::i,j
real(8),allocatable::hist(:)
real(8)::a
real(8)::bin

n=50
bin=1.d1
allocate(hist(n))
hist=0

10 read(5,*,end=99)i,a
j=a/bin+1
if(j.gt.0 .and. j.le.n)hist(j) = hist(j)+1.d0
goto 10

99 continue
do i=1,n
a=(i-1)*bin+bin*0.5
write(6,*)a,hist(i)
enddo

end
