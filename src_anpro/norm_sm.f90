program norm_sm

implicit none

integer,parameter :: nd=30000000
real*4 :: Fhead
real*4 :: sdata(nd)
real*4 :: wei(nd),weitemp,smwei(nd)
integer :: Nhead,npts
integer :: num,p1,p2,smnum,smph,smp
character*4 Khead
character*100 :: arg(10)
character*70 inpfile,outpfile
integer :: i,j,k,m,n,status1,cnvrtNum


do i=1,5
    call getarg(i,arg(i))
end do
inpfile=arg(1)
outpfile=arg(2)
read(arg(3),*)n
read(arg(4),*)smp
read(arg(5),*)cnvrtNum

write(*,*)"input file: ",inpfile
write(*,*)"output file: ",outpfile 

if ( cnvrtNum < 1) then
open(unit=22,file=inpfile,form="unformatted",&
    access="direct",recl=4,convert="swap")
else
open(unit=22,file=inpfile,form="unformatted",&
    access="direct",recl=4)
end if

open(unit=23,file=outpfile,form="unformatted",&
    access="direct",recl=4)

open(unit=24,file="log.txt")

do i=1,70
    read(22,rec=i)Fhead
    write(23,rec=i)Fhead
!    write(*,*)i,Fhead
end do

do i=71,110
    read(22,rec=i)Nhead
    if (i==80) then
        npts=Nhead
    end if
    write(23,rec=i)Nhead
!    write(*,*)i,Nhead
end do

do i=111,158,2
    read(22,rec=i)Khead
    write(23,rec=i)Khead
!    write(*,*)i,Khead
end do

i=159
write(*,*)inpfile,npts
! Read data
do while(.true.)
    read(22,rec=i,iostat=status1)sdata(i-158)
    if(status1/=0)exit
    i=i+1
end do

write(*,*)i,npts,"number of data points"
! Calculate weight
num=npts/n + 1
write(*,*)"line65",num
m=0
do j=1,num-1
    weitemp=0.0
    do k=1,n
        m=m+1
        weitemp=weitemp+abs(sdata(m))
    end do
    weitemp=weitemp/n

    p1=(j-1)*n+1
    p2=j*n
    wei(p1:p2)=weitemp
end do
m=(num-1)*n+1
weitemp=0.0
do j=m,npts
    weitemp=weitemp+abs(sdata(m))
end do
weitemp=weitemp/(npts-m+1)
wei(m:npts)=weitemp

! Smooth weight
smph=smp/2
do j=1,npts
    p1=j-smph
    p2=j+smph
    smnum=smp
    if(p1<1)then
        p1=1
        p2=j-p1+j
        smnum=p2-p1+1
    else if(p2 > npts)then
        p2=npts
        p1=j-(npts-j)
        smnum=p2-p1+1
    end if
    smwei(j)=sum(wei(p1:p2))/real(smnum)
end do


write(*,*)"line 103"
! Normalize data and write file
i=159
do j=1,npts
    !write(*,*)"line108",j,sdata(j),smwei(j),wei(j)
    sdata(j)=sdata(j)/smwei(j)
    !if ( abs(sdata(j)) > 10 ) then
    !    write(*,*)"i datus: ",i,sdata(j),smwei(j),wei(j)
   !     sdata3=sdata3/abs(sdata3)
    !end if

    write(23,rec=i)sdata(j)
    i=i+1
end do

!write(*,*)i

stop
end
