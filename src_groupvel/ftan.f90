program ftan

! This program does time-frequency analysis with Surface waves.
! The basic method is Gabor transform, but using group velocity
! c and period tp instead of frequency and time.

implicit none

real*4 :: Vel(1000),t(1000)
real*4 :: Period(1000),w(1000)
real*4 :: sacdata(15000)
real*4 :: TfanReal(1000,1000),TfanImg(1000,1000)
real*4 :: Tfanamp(1000,1000)
real*4 :: TfanampN(1000,1000)
real*4 :: NormP
real*4 :: Realpart,Imgpart,Gpart
real*4 :: Dist,DeltaT,Npts
real*4 :: MinPeriod,MaxPeriod,DPeriod
real*4 :: DVel
real*4 :: pi
real*4 :: a,BeginT
real*4 :: Fhead
real*4 :: vgroup
real*4 :: tempAmp
integer :: ivtr
integer :: Maxipf
integer :: Maxivt
integer :: ivt,ipf,idata,j,status1,cnvrtNum
integer :: Nhead
character*4 Khead
character*70 inpfile,outwhole,outgroup
parameter(pi=3.1415926535)



! Read input parameters. Open sac file and output file
!---------------------------------------------
open(unit=21,file='InputPrm.txt')
read(21,*)Dist
read(21,*)BeginT
read(21,*)cnvrtNum
read(21,*)MinPeriod,MaxPeriod,DPeriod
read(21,*)DVel
read(21,*)inpfile
read(21,*)outwhole
read(21,*)outgroup

if ( cnvrtNum < 1) then
open(unit=22,file=inpfile,form="unformatted",&
    access="direct",recl=4,convert="swap")
else
open(unit=22,file=inpfile,form="unformatted",&
    access="direct",recl=4)
end if

open(unit=23,file=outwhole)
open(unit=25,file=outgroup)
!open(unit=24,file='log.txt')
!----------------------------------------------


! Parameters: group velocity, time, periods, frequency
!---------------------------------------------

Maxivt=nint(((6.5-2.0)/DVel)+1)
do ivt=1,Maxivt
    Vel(ivt)=2.0+(ivt-1)*DVel
    t(ivt)=Dist/Vel(ivt)-BeginT
!    write(*,*)Vel(ivt)
end do

Maxipf=nint(((MaxPeriod-MinPeriod)/DPeriod)+1)
do ipf=1,Maxipf
    Period(ipf)=(ipf-1)*DPeriod+MinPeriod
    w(ipf)=2*pi/Period(ipf)
end do

!---------------------------------------------


! Read headers of sacfile
!----------------------------------------------
do idata=1,70
    read(22,rec=idata)Fhead
    if ( idata == 1 ) then
        DeltaT=Fhead
    end if
end do

do idata=71,110
    read(22,rec=idata)Nhead
    if ( idata == 80 ) then
        Npts=Nhead
    end if
end do

do idata=111,158,2
    read(22,rec=idata)Khead
end do

!----------------------------------------------



! Time-Frenquency analysis
!-----------------------------------------------------------------

    
do ipf=1,Maxipf
    write(*,*)"Period: ",ipf
    a=sqrt(log(2.0))/Period(ipf)
    do ivt=1,Maxivt   
        idata=159
        j=1
        TfanReal(ivt,ipf)=0
        TfanImg(ivt,ipf)=0
!        write(*,*)"ivt",ivt
        do while(.true.)
            read(22,rec=idata,iostat=status1)sacdata(j)
            Realpart=cos(w(ipf)*(j-1)*DeltaT)
            Imgpart=-sin(w(ipf)*(j-1)*DeltaT)
            Gpart=exp(-(a*(t(ivt)-(j-1)*DeltaT))**2)
            TfanReal(ivt,ipf)=TfanReal(ivt,ipf)+sacdata(j)*RealPart*Gpart
            TfanImg(ivt,ipf)=TfanImg(ivt,ipf)+sacdata(j)*ImgPart*Gpart
!            write(24,*)j,ivt,ipf,Realpart,Imgpart,Gpart,TfanReal(ivt,ipf),&
!                       TfanImg(ivt,ipf),sacdata(j)
            j=j+1
            idata=idata+1
            if(status1/=0)exit
        end do
        TfanAmp(ivt,ipf)=sqrt(TfanReal(ivt,ipf)**2+TfanImg(ivt,ipf)**2)
    end do
end do

!------------------------------------------------------------------

do  ipf=Maxipf,1,-1
    vgroup=2.0
    tempAmp=TfanAmp(1,ipf)
    do ivt=1,Maxivt-1
        if(TfanAmp(ivt+1,ipf) .gt. tempAmp)then
            tempAmp=TfanAmp(ivt+1,ipf)
            vgroup=ivt*DVel+2.0
        end if
    end do
    write(25,*)MinPeriod+(ipf-1)*DPeriod,vgroup,tempAmp
end do

! Write results.
!------------------------------------------------------------------
do ipf=1,Maxipf
    NormP=maxval(TfanAmp(:,ipf))
    write(*,*)NormP
    do ivt=1,Maxivt
        TfanAmpN(ivt,ipf)=TfanAmp(ivt,ipf)/NormP
        write(23,*)Vel(ivt),Period(ipf),TfanAmpN(ivt,ipf),&
            TfanAmp(ivt,ipf),TfanReal(ivt,ipf),TfanImg(ivt,ipf)
    end do
end do

!------------------------------------------------------------------
close(23)
close(25)

stop
end
