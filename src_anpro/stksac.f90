program stksac

!This Prog is used for stacking two station's waveform in time domain
!The number of two data must be same

implicit none

integer, parameter :: maxfile=1000
real*4 :: Fhead
real*4 :: sacdt(maxfile)
real*4 :: sacdtsum
integer :: Nhead
character*4 Khead
character*8 KEVNM
character*70 sacfile(maxfile)
character*70 ipf
integer :: i,j,cvtpara,filepara,n,status1,filenum

write(*,*)"inputfile1,inputfile2,integer k1,k2(k<1,convert):"
read(*,*)ipf,cvtpara

!#########################################################
!Input sacfiles' names which need to be stacked
open(unit=21,file=ipf)

filenum=0
do while(.true.)
    filenum=filenum+1
    read(21,*,iostat=status1)sacfile(filenum)
    if(status1/=0)exit
end do
filenum=filenum-1
write(*,*)filenum," sac files are going to be stacked!"
!#############################################################


!##############################################################
!Open sacfiles, depend on the cvtpara
if ( cvtpara < 1 ) then
!If data block in the binary file end in low bit,use convert
!model to read
    do i=22,21+filenum 
        open(unit=i,file=sacfile(i-21),form='unformatted',&
            access='direct',recl=4,convert='swap')
    end do
else
!If data block in the binary file end in high bit,use non-convert
!model to read
    do i=22,21+filenum
        open(unit=i,file=sacfile(i-21),form='unformatted',&
            access='direct',recl=4)
    end do
end if
!##############################################################



open(unit=22+filenum,file='sum.sac',form='unformatted',&
    access='direct',recl=4)

!#############################################################
!Read head variables of sacfiles, use the first sacfile's head
!var as the sumfile's.
do i=1,70
    read(22,rec=i)Fhead
    write(22+filenum,rec=i)Fhead
end do

do i=71,110
    read(22,rec=i)Nhead
    write(22+filenum,rec=i)Nhead
end do

do i=111,157,2
    read(22,rec=i)Khead
    write(22+filenum,rec=i)Khead
end do
!###############################################################

i=159
do while(.true.)
    do j=1,filenum
         read(21+j,rec=i,iostat=status1)sacdt(j)
    end do
    if(status1/=0)exit
    sacdtsum=0
    do j=1,filenum
        sacdtsum=sacdtsum+sacdt(j)
    end do
    write(22+filenum,rec=i)sacdtsum
    i=i+1
end do
    
stop
end
