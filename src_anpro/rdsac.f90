program rdsac

implicit none

real*4 :: Fhead
real*4 :: sacdata
integer :: Nhead
character*4 Khead
character*8 KEVNM
character*70 sacfile
integer :: i,j,k,n,status1

open(21,file='rdsac.inp',status='old')
read(*,*)sacfile
read(*,*)n
read(*,*)k


if ( k < 1 ) then
!If data block in the binary file end in low bit,use convert
!model to read
open(unit=22,file=sacfile,form='unformatted',&
    access='direct',recl=4,convert='swap')
else
!If data block in the binary file end in high bit,use non-convert
!model to read
open(unit=22,file=sacfile,form='unformatted',&
    access='direct',recl=4)
end if

do i=1,70
    read(22,rec=i)Fhead
    write(*,*)i,Fhead
end do

do i=71,110
    read(22,rec=i)Nhead
    write(*,*)i,Nhead
end do

do i=111,158,2
    read(22,rec=i)Khead
    write(*,*)i,Khead
end do

do i=159,n+158
    read(22,rec=i)sacdata
    write(*,*)i,sacdata
end do

stop
end
