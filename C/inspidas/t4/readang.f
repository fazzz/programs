	real*8 a(10)
	integer n(10)
	character*3 r

 	open(10,file='angreg0')
c	open(10,file='angfile')

	read(10,*)
	read(10,*)

	do 10 i=1,164

		read(10,11)j,r,(n(k),a(k),k=1,10)
		write(6,21)(a(k),k=1,10)

10	continue
11	format(i5,1x,a3,1x,10(i3,f10.4))
21	format(10f8.3)

	end
