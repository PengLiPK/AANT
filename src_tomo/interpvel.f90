program interpvel

use strct
implicit none
type(vstrct) :: v(maxgrid)
type(tstrct), target :: vorig(maxgrid)
type(pstrct), pointer :: prv(:)
real(kind=8) :: dxv,dyv
real(kind=8) :: minx,maxx
real(kind=8) :: miny,maxy
real(kind=8) :: intx,inty
integer :: vxnum,vynum
integer :: ttlgnum
integer :: vorigxnum,vorigynum
integer :: vorignum
integer :: ix,iy
integer :: iv
character(len=70) :: velfile


open(22,file="interpvel.inp")
read(22,*)velfile
read(22,*)minx,maxx
read(22,*)miny,maxy
read(22,*)intx,inty
read(22,*)dxv,dyv

vorigxnum=nint((maxx-minx)/dxv)+1
vorigynum=nint((maxy-miny)/dyv)+1

vxnum=nint((maxx-minx)/intx)+1
vynum=nint((maxy-miny)/inty)+1
ttlgnum=vxnum*vynum
do iy=1,vynum
    do ix=1,vxnum
        v(ix+(iy-1)*vxnum)%x=minx+(ix-1)*intx
        v(ix+(iy-1)*vynum)%y=miny+(iy-1)*inty
    end do
end do

open(29,file=velfile)
read(29,*)vorignum
do iv=1,vorignum
    read(29,*)vorig(iv)%x,vorig(iv)%y,vorig(iv)%t
    vorig(iv)%dxx=0
    vorig(iv)%num=iv
    vorig(iv)%stat=0
end do
close(29)

open(30,file='velxy.txt')
write(30,*)ttlgnum
do iv=1,ttlgnum

    call localcood(v(iv)%x,v(iv)%y,vorig,prv,&
        &dxv,dyv,minx,miny,vorigxnum)

    v(iv)%vel=bilinear(prv(1)%p%t,prv(2)%p%t,prv(3)%p%t,prv(4)%p%t,&
          &prv(1)%p%x,prv(1)%p%y,prv(4)%p%x,prv(4)%p%y,&
          &v(iv)%x,v(iv)%y)
    write(30,*)v(iv)%x,v(iv)%y,v(iv)%vel,iv
end do
write(*,*)"line579"
close(30)

stop
end
