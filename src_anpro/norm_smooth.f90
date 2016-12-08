program norm_smooth

implicit none

real*4 :: Fhead
real*4 :: sacdata(5000)
real*4 :: abssume
integer :: Nhead,npts
character*4 Khead
character*8 KEVNM
character*70 inpfile,outpfile
integer :: i,j,k,m,n,status1,cnvrtNum

read(*,*)inpfile,outpfile,n,cnvrtNum
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
    if ( i .eq. 80) then
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
do while(.true.)
  abssume=0
  m=i
  do j=1,n
    read(22,rec=m,iostat=status1)sacdata(j)
    abssume=abssume+abs(sacdata(j))/n
    m=m+1
    if(m>npts)exit
  end do

  do k=1,2
    sacdata(k)=sacdata(k)/abssume
    !if ( abs(sacdata(k)) > 10 ) then
    !    write(*,*)"i datus: ",i,abs(sacdata(k))
    !    sacdata(k)=sacdata(k)/abs(sacdata(k))
    !end if
    write(23,rec=i)sacdata(k)
    i=i+1
  end do
  if(i>npts)exit
end do

!write(*,*)i

stop
end
