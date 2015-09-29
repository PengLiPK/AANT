program norm

implicit none

real*4 :: Fhead
real*4 :: sacdata(5000)
real*4 :: abssume
integer :: Nhead
character*4 Khead
character*8 KEVNM
character*70 inpfile,outpfile
integer :: i,j,k,n,status1,cnvrtNum

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
    write(23,rec=i)Nhead
!    write(*,*)i,Nhead
end do

do i=111,158,2
    read(22,rec=i)Khead
    write(23,rec=i)Khead
!    write(*,*)i,Khead
end do

i=159
!write(*,*)inpfile
do while(.true.)
  j=1
  abssume=0
  do while(.true.)
    read(22,rec=i,iostat=status1)sacdata(j)
    abssume=abssume+abs(sacdata(j))/n
!    write(*,"(i4,1x,G16.7)")i,sacdata
    j=j+1
    i=i+1
    if(j>n)exit
    if(status1/=0)exit
  end do
!  write(24,*)"j",j,"i",i,"weight",abssume
  i=i-j+1
!  write(24,*)"i",i
  do k=1,j-1
!    write(24,*)"i",sacdata(k)
    sacdata(k)=sacdata(k)/abssume
!One-bit normalization
!    sacdata(k)=sacdata(k)/abs(sacdata(k))
    if ( abs(sacdata(k)) > 10 ) then
        write(*,*)"i datus: ",i,abs(sacdata(k))
        sacdata(k)=sacdata(k)/abs(sacdata(k))
    end if
    write(23,rec=i)sacdata(k)
    i=i+1
  end do
!  write(24,*)"i",i
  if(status1/=0)exit
end do

!write(*,*)i

stop
end
