program cacuns

! This program calculates the signal noise ratio of ambient
! noise signal.

implicit none

real*4 :: sacdata(5000)
real*4 :: pi
real*4 :: Fhead
real*4 :: datasum,dataave,datarms
real*4 :: nlv,nlvzero
integer :: ivt,ipf,idata,j,status1,cnvrtNum
integer :: Nhead
character*4 Khead
character*70 inpfile,outpfile
parameter(pi=3.1415926535)



! Read input parameters. Open sac file and output file
!------------------------------------------------------
open(unit=21,file='cacuns.inp')
read(21,*)cnvrtNum
read(21,*)inpfile

if ( cnvrtNum < 1) then
    open(unit=22,file=inpfile,form="unformatted",&
        &   access="direct",recl=4,convert="swap")
else
    open(unit=22,file=inpfile,form="unformatted",&
        &   access="direct",recl=4)
end if

!---------------------------------------------------------

! Read headers of sacfile
!----------------------------------------------
do idata=1,70
    read(22,rec=idata)Fhead
end do

do idata=71,110
    read(22,rec=idata)Nhead
end do

do idata=111,158,2
    read(22,rec=idata)Khead
end do

!----------------------------------------------


j=1    
datasum=0
do idata=6000+148,6000+168
    read(22,rec=idata)sacdata(j)
    datasum=datasum+sacdata(j)
    j=j+1
end do

datarms=0
dataave=datasum/21.0
do idata=1,21
    datarms=datarms+abs(sacdata(idata)-dataave)
end do
nlvzero=datarms/21.0

j=1
datasum=0
do idata=159,958
    read(22,rec=idata)sacdata(j)
    datasum=datasum+sacdata(j)
    j=j+1
end do

datarms=0
dataave=datasum/800.0
do idata=1,800
    datarms=datarms+abs(sacdata(idata)-dataave)
end do
nlv=datarms/800.0

open(23,file="noise.out")
write(23,*)inpfile,nlv,nlvzero

close(23)

close(22)

stop
end
