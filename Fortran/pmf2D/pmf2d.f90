program pmftwod

  implicit none

  character(100) :: filename

  integer i, j, k, l, N
  integer num_bin_x, num_bin_y
  double precision bin_sizex, bin_sizey, maxx, minx, maxy, miny, maxpm
  double precision,allocatable,dimension(:,:) :: data, hist

  character(10) :: argv, programname
  integer :: argc, iargc

  NAMELIST/htinfo/bin_sizex, bin_sizey, N
  open(16,file='parameters_hist')
  read(16,htinfo)
  close(16)

  argc=iargc()
  if (argc < 1) then
     call getarg(0,argv)
     call usage(argv)
  else
     call getarg(1,filename)
  end if

  allocate(data(2, N))

  open(17,file=filename,status='old')
  do i=1,N,1
     read(17,'(F8.4,1x,F8.4)'),data(1,i), data(2,i)
  end do
  close(17)

  maxx=data(1,1); minx=data(1,1)
  do i=2,N,1
     if ( data(1,i) > maxx) then
        maxx=data(1,i)
     end if
     if ( data(1,i) < minx) then
        minx=data(1,i)
     end if
  end do

  maxy=data(2,1); miny=data(2,1)
  do i=2,N,1
     if ( data(2,i) > maxy) then
        maxy=data(2,i)
     end if
     if ( data(2,i) < miny) then
        miny=data(2,i)
     end if
  end do

  num_bin_x=int((maxx-minx)/bin_sizex) + 1
  num_bin_y=int((maxy-miny)/bin_sizey) + 1
  
  allocate(hist(num_bin_x, num_bin_y))

  hist=0.0d0

  do i=1,N,1
     k=int((data(1, i)-minx)/bin_sizex)+1
     l=int((data(2, i)-miny)/bin_sizey)+1
     hist(k, l) = hist(k, l) + 1
  end do

  maxpm = -1000.0d0
  do i=1,num_bin_x,1
     do j=1,num_bin_y,1
        if (hist(i,j) .ne. 0.0d0) then
           hist(i,j) = -1.0*log(hist(i,j)/bin_sizex/bin_sizey/N)

           if (maxpm < hist(i,j)) then
              maxpm = hist(i,j)
           end if
        end if
     end do
  end do

  do i=1,num_bin_x,1
     do j=1,num_bin_y,1
        if (hist(i,j) .ne. 0.0d0) then
           hist(i,j) = hist(i,j) - maxpm
        end if
     end do
  end do
  
  do i=1,num_bin_x,1
     do j=1,num_bin_y,1
        print *,minx+bin_sizex*i-bin_sizex*5.0e-1, miny+bin_sizey*j-bin_sizey*5.0e-1 ,hist(i,j)
     end do
     print *,' '
  end do
  
  deallocate(data)
  deallocate(hist)

  contains

  subroutine usage(argv)
    implicit none
    character(10) :: argv

    
    print "('Usage: ',A10,' filename')", argv
    stop
  end subroutine usage

end program pmftwod
