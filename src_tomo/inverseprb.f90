! This module defines some struct type variables which will be used in 
! main program and subroutine.
!----------------------------------------------------------------------

module typedef
    implicit none
    real(kind=8),parameter :: pi=3.14159265
    real(kind=8),parameter :: radii=6371.0
    real(kind=8),parameter :: e=1D-12
    integer,parameter :: imodel=500
    integer,parameter :: idata=500
end module typedef
!-----------------------------------------------------------------------



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

program inverseprb


! This program uses conjugate gradient method solving
! least square problem.
! The equation is [G'Cd^(-1)G+lambda*Cm^(-1)]m=G'Cd^(-1)d+lambda*Cm^(-1)m0.
! Inputs are travel time, grids and frech derivative.
! Outputs are velocity models.




use typedef
implicit none
real(kind=8) :: G(idata,imodel)
real(kind=8) :: Gedge(idata,imodel)
real(kind=8) :: MinLon, MaxLon
real(kind=8) :: MinLan, MaxLan
real(kind=8) :: InterLon
real(kind=8) :: InterLan
real(kind=8) :: m(imodel)
real(kind=8) :: d(idata)
real(kind=8) :: s(idata)
real(kind=8) :: r(imodel)
real(kind=8) :: rabs(imodel)
real(kind=8) :: p(imodel)
real(kind=8) :: Gp(idata)
real(kind=8) :: Alpha, AlphaUp, AlphaDown1,AlphaDown2
real(kind=8) :: Beta, BetaUp, BetaDown
real(kind=8) :: rmax
real(kind=8) :: m0(imodel),nodex(imodel),nodey(imodel)
real(kind=8) :: edgev(imodel),edgex(imodel),edgey(imodel)
real(kind=8) :: lambda
real(kind=8) :: cdval,cmval
real(kind=8) :: cdinv(idata),cminv(imodel)
real(kind=8) :: minm,maxm
real(kind=8) :: vconst
integer :: totalnodenum,edgenodenum
integer :: GridNum
integer :: NumLon
integer :: NumLan
integer :: ig
integer :: id, idmax
integer :: im
integer :: k
integer :: ix,iy
character*70 :: datafile
character*70 :: GridFile
character*70 :: FdFile
character*70 :: OputFile
character*70 :: oputmodel
character*70 :: inpmodel
character*30 :: tempch


! Input parameters
!----------------------------------------------------

open(22,file='inverseprb.inp')
read(22,*)FdFile
read(22,*)datafile
read(22,*)OputFile
read(22,*)oputmodel
write(*,*)oputmodel
read(22,*)inpmodel
write(*,*)inpmodel
read(22,*)idmax
read(22,*)MinLon,MaxLon
read(22,*)MinLan,MaxLan
read(22,*)vconst
read(22,*)lambda
read(22,*)cdval,cmval
close(22)

write(*,*)inpmodel

open(32,file=inpmodel)
read(32,*)totalnodenum
read(32,*)edgenodenum
GridNum=totalnodenum-edgenodenum
do im=1,edgenodenum
    read(32,*)edgex(im),edgey(im),edgev(im)
end do

do im=1,GridNum
    read(32,*)nodex(im),nodey(im),m0(im)
    write(*,*)m0(im),im
    m(im)=1.0/m0(im)
    cminv(im)=1.0/cmval
end do
close(32)


open(25,file=FdFile)
! Caculate G.
!------------------------------------------------------------

do id=1,idmax
    do im=1,edgenodenum
        read(25,*)Gedge(id,im)
    end do
    do im=1,GridNum
        read(25,*)G(id,im)
    end do
end do


!------------------------------------------------------------

open(26,file=datafile)
do id=1,idmax
    read(26,*)d(id)
    do im=1,edgenodenum
        d(id)=d(id)-Gedge(id,im)/vconst
    end do
    cdinv(id)=1.0/cdval
    write(*,*)d(id)
end do



open(28,file='inverse2.log')
!------------------------------------------------------------



! Giving initial value s0, r0, p0
!-----------------------------------------------------------
!write(*,*)(cdinv(id),id=1,idmax)
!write(*,*)(cminv(im),im=1,GridNum)
!pause

do id=1,idmax
   s(id)=d(id)
   do im=1,GridNum
       s(id)=s(id)-G(id,im)*m(im)
   end do
end do

r=0
do im=1,GridNum
    do id=1,idmax
        r(im)=r(im)+G(id,im)*cdinv(id)*s(id)
    end do
end do

do im=1,GridNum
   rabs(im)=abs(r(im))
   p(im)=r(im)
end do

rmax=maxval(rabs(1:GridNum))
minm=minval(m(1:GridNum))
maxm=maxval(m(1:GridNum))

!-------------------------------------------------------



open(27,file=OputFile)
! Caculating the best m
!--------------------------------------------------------

k=0
    write(*,*)k
    write(27,*)k
    write(*,*)(1/m(im),im=1,GridNum)
    write(27,*)(1/m(im),im=1,GridNum)
    write(*,*)(r(im),im=1,GridNum)
    write(27,*)(r(im),im=1,GridNum)
    write(27,*)(rabs(im),im=1,GridNum)
    write(*,*)rmax
    write(27,*)rmax
do while( rmax>e .and. ((minm>0.1) .and. (maxm<1)))
!do while( rmax>e )

    write(*,*)k
    write(27,*)k
    !write(*,*)(1/m(im),im=1,GridNum)
    write(27,*)(1/m(im),im=1,GridNum)
    write(27,*)(r(im),im=1,GridNum)
    write(27,*)(rabs(im),im=1,GridNum)
    write(*,*)rmax
    write(27,*)rmax

    k=k+1

! Caculating Alpha(k)
    do id=1,idmax
        Gp(id)=0
        do im=1,GridNum
            Gp(id)=Gp(id)+G(id,im)*p(im)
        end do
    end do
    
    AlphaUp=0
    do im=1,GridNum
        AlphaUp=AlphaUp+r(im)*r(im)
    end do

    AlphaDown1=0
    do id=1,idmax
        AlphaDown1=AlphaDown1+Gp(id)*cdinv(id)*Gp(id)
    end do
    AlphaDown2=0
    do im=1,GridNum
        AlphaDown2=AlphaDown2+p(im)*cminv(im)*p(im)
    end do
    AlphaDown2=lambda*AlphaDown2
    

    Alpha=AlphaUp/(AlphaDown1+AlphaDown2)

!    write(*,*)"Alpha",Alpha
!    pause

! Caculating m(k+1), r(k+1), Beta(k+1), will use Gp upstairs.
    do im=1,GridNum
        m(im)=m(im)+Alpha*p(im)
    end do
    minm=minval(m(1:GridNum))
    maxm=maxval(m(1:GridNum))


    BetaDown=0
    do im=1,GridNum
        BetaDown=BetaDown+r(im)*r(im)
    end do

    do im=1,GridNum
        do id=1,idmax 
            r(im)=r(im)-Alpha*G(id,im)*cdinv(id)*Gp(id)
        end do
        r(im)=r(im)-Alpha*lambda*cminv(im)*p(im)
    end do


    do im=1,GridNum
       rabs(im)=abs(r(im))
    end do

    rmax=maxval(rabs(1:GridNum))

    BetaUp=0
    do im=1,GridNum
        BetaUp=BetaUp+r(im)*r(im)
    end do
    
    Beta=BetaUp/BetaDown
    write(*,*)"Beta",Beta

! Caculating p(k+1)
    do im=1,GridNum
        p(im)=r(im)+Beta*p(im)
    end do
end do

!-----------------------------------------------------------


close(25)
close(26)
close(27)
close(28)

open(33,file=oputmodel)
do im=1,GridNum
    write(33,*)1/m(im)
end do
close(33)

open(29,file="velresult.txt")
write(29,*)totalnodenum
write(29,*)edgenodenum
do iy=1,edgenodenum
    write(29,*)edgex(iy),edgey(iy),edgev(iy)
end do
do iy=1,GridNum
        write(29,*)nodex(iy),nodey(iy),1/m(iy)
end do
close(29)


end


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

