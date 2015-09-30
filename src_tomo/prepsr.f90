! This program generates random points with giving a 
! retangle area.
! Input: pnum: the number of points
!        minx,maxx,miny,maxy: boundaries of retangle area
! Output: x,y: the coordinates of random points
program prepsr

implicit none
real(kind=8) :: minx,maxx
real(kind=8) :: miny,maxy
real(kind=8) :: xscale,yscale
real(kind=8) :: xcenter,ycenter
real(kind=8) :: xrandom,yrandom
real(kind=8) :: x,y
real(kind=8),allocatable :: r(:,:)
integer,allocatable :: seed(:)
integer :: pnum
integer :: i,n,clock
character(len=70) :: outfile



open(21,file='prepsr.inp')
read(21,*)outfile
read(21,*)pnum
read(21,*)minx,maxx
read(21,*)miny,maxy
close(21)

allocate(r(pnum,2))

xscale=maxx-minx
yscale=maxy-miny


call random_seed(size=n)
allocate(seed(n))
!write(*,*)n
call system_clock(count=clock)
!write(*,*)clock
seed=clock+(/(i-1,i=1,n)/)
call random_seed(put=seed)
deallocate(seed)
call random_number(r)
!write(*,*)r


open(23,file=outfile)
write(23,*)pnum
do i=1,pnum
    x=r(i,1)*xscale+minx
    y=r(i,2)*yscale+miny
    write(23,*)x,y
end do

stop
end

