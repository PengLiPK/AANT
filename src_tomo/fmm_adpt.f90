! FMM with Delaunay triangles
!----------------------------------------------------------------

program fmm

! Use file fmm.inp, source file.

use strct
use linklist
implicit none
type(node_linklist),pointer :: node_head,node_tail
type(node_linklist),pointer :: node_p
type(srstrct) :: source(maxsource)
type(rcstrct) :: receiver(maxreceiver)
type(tstrct) :: v(maxgrid)
type(pstrct),pointer :: ptri(:)
type(gnstrct),target :: gn(maxvel)
type(pgnstrct),pointer :: pgn(:),pgntemp
type(srstrct) :: vtx(maxsource)
real(kind=8) :: g(maxdata,maxvel)
real(kind=8) :: gnmax,gnsum,gnavr
real(kind=8) :: dx,dy
real(kind=8) :: minx,maxx
real(kind=8) :: miny,maxy
real(kind=8) :: newndx(maxtri),newndy(maxtri)
real(kind=8) :: start,mid,finish
real(kind=8) :: rfxrg,rfyrg
real(kind=8) :: rfxint,rfyint
real(kind=8) :: xtemp,ytemp,veltemp,stemp
real(kind=8) :: vconst,voutside
real(kind=8) :: thrshd0,thrshd1,thrshd2
real(kind=8) :: fqthr,period
integer :: tri(maxtri,3)
integer :: trinum
integer :: nvtx
integer :: hwnodenum(maxvel)
integer :: hwnum,hwmax,newndnum
integer :: edgenode
integer :: ns,is
integer :: nr,ir
integer :: nt,it
integer :: vnum
integer :: gridtype
integer :: imethod
integer :: icoord
integer :: iv
integer :: id,im
integer :: idmax,immax
integer :: date(3)
integer :: time(3)
character(len=70) :: sourcefile
character(len=70) :: receiverfile
character(len=70), allocatable :: rfilename(:)
character(len=70) :: velfile
character(len=70) :: trifile
character(len=70) :: newndfile
character(len=70) :: vtxfile

call cpu_time(start)
call idate(date)
call itime(time)
write(*,2000)date,time
2000 format('Date: ',i2.2,'/',i2.2,'/',i4.4,'; Time: ',&
& i2.2,':',i2.2,':',i2.2)


! Read parameters.
!-----------------------------------------------------------------------
open(21,file='fmm_adpt.inp',status='old')
read(21,*)period
read(21,*)sourcefile
read(21,*)receiverfile
read(21,*)velfile
read(21,*)trifile
read(21,*)ns
read(21,*)idmax
read(21,*)dx,dy
read(21,*)rfxint,rfyint
read(21,*)rfxrg,rfyrg
read(21,*)minx,maxx
read(21,*)miny,maxy
read(21,*)thrshd0,thrshd1,thrshd2
read(21,*)imethod
read(21,*)icoord
read(21,*)gridtype
read(21,*)vconst
read(21,*)nvtx
read(21,*)vtxfile
read(21,*)voutside
close(21)

open(22,file=sourcefile,status='old')
do is=1,ns
    read(22,*)source(is)%x,source(is)%y
    if((source(is)%x .gt. maxx) .or. &
      &(source(is)%x .lt. minx) .or. &
      &(source(is)%y .gt. maxy) .or. &
      &(source(is)%y .lt. miny))then
        write(*,*)"Number ",is," source is outside study area!!"
        stop
    end if
end do
close(22)

allocate(rfilename(ns))

open(23,file=receiverfile,status='old')
do ir=1,ns
    read(23,*)rfilename(ir)
end do
close(23)

open(24,file=vtxfile,status='old')
do iv=1,nvtx
    read(24,*)vtx(iv)%x,vtx(iv)%y
end do



open(25,file=velfile,status='old')
read(25,*)vnum
read(25,*)edgenode
do iv=1,vnum
    read(25,*)v(iv)%x,v(iv)%y,v(iv)%t
    v(iv)%dxx=0
    v(iv)%num=iv
    v(iv)%stat=0
end do
close(25)

open(26,file=trifile,status='old')
read(26,*)trinum
do it=1,trinum
    read(26,*)tri(it,1),tri(it,2),tri(it,3)
end do
close(26)

!-----------------------------------------------------------------------


! Calcu traveltimes, ray paths and Frechet derivative
!----------------------------------------------------------------------
open(27,file='traveltime.txt',status='replace')
open(28,file="fd.txt",status='replace')
open(29,file='traveltimbyfd.txt',status='replace')

do is=1,ns
    open(30,file=rfilename(is),status='old')
    read(30,*)nr
    do ir=1,nr
        read(30,*)receiver(ir)%x,receiver(ir)%y
        if((receiver(ir)%x .gt. maxx) .or. &
          &(receiver(ir)%x .lt. minx) .or. &
          &(receiver(ir)%y .gt. maxy) .or. &
          &(receiver(ir)%y .lt. miny))then
            write(*,*)"Number ",ir," receiver is outside study area!!"
            stop
        end if
    end do
    close(30)


    call cacut(trinum,tri,is,source(is),receiver,v,vnum,nr,dx,dy,minx,&
    &maxx,miny,maxy,nt,rfxint,rfyint,rfxrg,rfyrg,imethod,icoord,gridtype,&
    &nvtx,vtx,voutside)
    do ir=1,nr
        write(27,*)receiver(ir)%t,receiver(ir)%x,receiver(ir)%y,&
        &   source(is)%x,source(is)%y
    end do
!    call system('cp source00001.out rename.out')
!    call rename('source00001.out','rename.out')
    write(*,*)"source",is,"has finished."
    write(*,*)nt,"grids has been caculated."
end do

close(27)
close(29)
!-------------------------------------------------------------------------


call cpu_time(mid)
write(*,2001)start,mid,mid-start


! Adjust nodes to increase the uniqueness of solution
!------------------------------------------------------------------------
! Varibles: idmax: max data
!           immax: number of velocity model nodes
!           g: 'G' in Gm=d
!           node_head: head of velocity nodes linklist, the head dosen't
!                      store velocity node
!           node_p: pointer point to the velocity nodes linklist
!           id,im,iv: Counters
!           gn: 1st norm of g
!           gnsum,gnavr,gnmax: sum, average and max of gn
!-------------------------------------------------------------------------
immax=vnum
! Read frechet derivative and grid node
rewind(28)
do id=1,idmax
    do im=1,immax
        read(28,*)g(id,im)
    end do
end do
close(28)


open(51,file=velfile,status='old')
read(51,*)vnum
read(51,*)edgenode
allocate(node_head)
allocate(node_tail)
node_head%next=>node_tail
node_tail%prev=>node_head
nullify(node_head%prev)
nullify(node_tail%next)

do iv=1,immax
    read(51,*)xtemp,ytemp,veltemp
    call list_ins(node_tail,xtemp,ytemp,veltemp,2)
end do
close(51)

allocate(pgn(immax))
! Ray density of each node from G
open(52,file='G_1st_norm.txt',status='replace')
gnsum=0
do im=1,immax
    gn(im)%val=0
    do id=1,idmax
        if(g(id,im) .gt. 0)then
            gn(im)%val=gn(im)%val+1
        end if
    end do
    gn(im)%num=im
    pgn(im)%p=>gn(im)
    write(*,*)gn(im)
    write(52,*)gn(im)
    if(im .gt. edgenode)then
        gnsum=gnsum+gn(im)%val
    end if
end do
gnavr=gnsum/(immax-edgenode)
write(*,*)"Mean: ",gnavr
write(52,*)"Mean: ",gnavr
close(52)

gnmax=maxval(gn(1:immax)%val)
write(*,*)"maxgn: ",gnmax
!call list_output(node_head)

! Delete nodes
open(52,file='edge_delgrid.txt',status='replace')
open(53,file='delgrid.txt',status='replace')
node_p=>node_head%next
do im=1,immax
    if(gnmax .eq. 0)then
        write(*,*)"All entries of matrix are 0!! &
        &Stop on line 169 of fmmrf2od_v7.f90"
        stop
    else if(gn(im)%val .lt. thrshd1)then
        write(*,*)"line170m,node,gn",im,gn(im)
        if((abs(node_p%x-minx) .lt. 1D-10) .or.&
          &(abs(node_p%x-maxx) .lt. 1D-10) .or.&
          &(abs(node_p%y-miny) .lt. 1D-10) .or.&
          &(abs(node_p%y-maxy) .lt. 1D-10))then
            write(52,*)node_p%x,node_p%y,vconst
        end if
        write(53,*)node_p%x,node_p%y,gn(im)
        call list_del(node_p)
    else
        node_p=>node_p%next
    end if
end do
close(52)
close(53)
newndfile='screen'
call list_output(node_head,newndfile)

write(*,*)"idmax,im,thrshd0",idmax,immax,thrshd0

! Add nodes when nodes number less than threshd0*(number of data)
if(immax .lt. nint(real(idmax)*thrshd0))then
    ! Sort pgn, which will be used in adding nodes
    allocate(pgntemp)
    do im=1,immax-1
        do id=im+1,immax
            if(pgn(id)%p%val .gt. pgn(im)%p%val)then
                pgntemp%p=>pgn(im)%p
                pgn(im)%p=>pgn(id)%p
                pgn(id)%p=>pgntemp%p
            end if
        end do
    end do
    deallocate(pgntemp)
    do im=1,immax
        write(*,*)"line254",pgn(im)%p
    end do
    do im=1,immax
        write(*,*)"line262",gn(im)
    end do

    ! Add nodes
    if((nint(real(idmax)*thrshd0)-immax) .gt. immax/2)then
        hwmax=immax/2
    else
        hwmax=idmax/4-immax
    end if
    hwnum=0
    do im=1,hwmax
        if(pgn(im)%p%val .gt. thrshd2)then
            hwnum=hwnum+1
            hwnodenum(hwnum)=pgn(im)%p%num
            write(*,*)hwnodenum(hwnum)
        end if
    end do

    fqthr=period*2.6
    call addnodes(newndx,newndy,newndnum,hwnodenum,hwnum,gn,v,tri,trinum,fqthr)

    deallocate(pgn)

    allocate(ptri(1:3))
    do im=1,newndnum
        write(*,*)newndx(im),newndy(im)
        call baryc_locat(newndx(im),newndy(im),v,tri,&
             &trinum,ptri)
        call baryc(stemp,1.0/ptri(1)%p%t,1.0/ptri(2)%p%t,1.0/ptri(3)%p%t,&
             &ptri(1)%p%x,ptri(1)%p%y,ptri(2)%p%x,ptri(2)%p%y,&
             &ptri(3)%p%x,ptri(3)%p%y,newndx(im),newndy(im))
        veltemp=1.0/stemp
        call list_ins(node_tail,newndx(im),newndy(im),veltemp,2)
    end do
    deallocate(ptri)
end if

newndfile='newnode.txt'
call list_output(node_head,newndfile)

!-------------------------------------------------------------------


call idate(date)
call itime(time)
write(*,2000)date,time

call cpu_time(finish)
write(*,2001)start,finish,finish-start
2001 format('start: ',f8.4,'; finish: ',f8.4,&
&'; Time consume: ',f8.4,' s.')

stop

end

!-----------------------------------------------------------------------



!This subroutine caculate travel time for each grid.
!--------------------------------------------------------------------------
subroutine cacut(trinum,tri,isr,source,receiver,vinp,vnum,nr,dx,dy,minx,&
  &maxx,miny,maxy,ttlgnum,rfxint,rfyint,rfxrg,rfyrg,imethod,icoodnt,gridtype,&
  &nvtx,vtx,voutside)
use strct
implicit none
type(srstrct) :: rcenter
type(srstrct) :: source
type(rcstrct) :: receiver(maxreceiver)
type(srstrct) :: path(maxpathnode)
type(tstrct) :: vinp(maxgrid)
type(tstrct), target :: vorig(maxgrid)
type(pstrct), pointer :: prv(:)
type(pstrct), pointer :: ptri(:)
type(tstrct), target :: travelt(maxgrid)
type(tstrct), target :: rtravelt(maxgrid)
type(pstrct), pointer :: ptravelt(:)
type(pstrct), pointer :: prtravelt(:)
type(pstrct), pointer :: ptemp
type(srstrct) :: vtx(maxsource)
real(kind=8) :: v(maxgrid)
real(kind=8) :: rv(maxgrid)
real(kind=8) :: voutside
real(kind=8) :: rvtemp
real(kind=8) :: fd(maxvel)
real(kind=8) :: dxx(maxpathnode)
real(kind=8) :: dyy
real(kind=8) :: rdxx(maxpathnode)
real(kind=8) :: rdyy
real(kind=8) :: dx,dy
real(kind=8) :: rdx,rdy
real(kind=8) :: minx,maxx
real(kind=8) :: miny,maxy
real(kind=8) :: rminx,rmaxx
real(kind=8) :: rminy,rmaxy
real(kind=8) :: rfxrg,rfyrg
real(kind=8) :: rfrg(4)
real(kind=8) :: rfxint,rfyint
real(kind=8) :: tempara1,tempara2
real(kind=8) :: dxv,dyv
real(kind=8) :: temptt
integer :: nvtx,rst
integer :: nb(maxnbnode),nbnode(maxnbnode),nbtail
integer :: tempnbnode(maxnbnode),tempnbnum
integer :: nbnum,inb
integer :: gridtype
integer :: tri(maxtri,3)
integer :: trinum
integer :: vxnum,vynum
integer :: vnum
integer :: imethod
integer :: icoodnt
integer :: ist,iist
integer :: minxl,maxxl
integer :: rminxl,rmaxxl
integer :: minyl,maxyl
integer :: rminyl,rmaxyl
integer :: xnum,ynum
integer :: rxnum,rynum
integer :: ttlgnum,rttlgnum
integer :: rcnum
integer :: gnum,rgnum
integer :: ig
integer :: iv
integer :: ip,iip,rip
integer :: ix,iy
integer :: isr
integer :: nr,ir
integer :: ipth,ippth
integer :: iipnum
integer :: i,j,k
character(len=10) :: sourcenum
character(len=70) :: tfilename
character(len=70) :: pathfname
character(len=70) :: errfname


vorig=vinp
write(*,*)"vnum",vnum
! initial grid are all far away
travelt%stat=-1
travelt%nbstat=0
rtravelt%stat=-1
travelt%nbstat=0

allocate(ptravelt(1:maxgrid))
allocate(prtravelt(1:maxgrid))
allocate(prv(1:4))
allocate(ptri(1:3))



minxl=1
maxxl=nint((maxx-minx)/dx)+1
minyl=1
maxyl=nint((maxy-miny)/dy)+1
xnum=maxxl-minxl+1
ynum=maxyl-minyl+1
ttlgnum=xnum*ynum



! Choose coodinate for outter grids.
!--------------------------------------------------------
if(icoodnt .eq. 1)then
    dxx=(pi*radii)*dx/180.0
    dyy=(pi*radii)*dy/180.0
else if(icoodnt .eq. 2)then
    dyy=(pi*radii)*dy/180.0
    do iy=1,ynum
        dxx(iy)=cos((miny+((real(iy)-0.5)*dy))*pi/180.0)&
        &   *pi*radii*dx/180.0
    end do
end if

do ig=1,ttlgnum
    travelt(ig)%x=minx+dx*(ig-((ig-1)/xnum)*xnum-1)
    travelt(ig)%y=miny+dy*((ig-1)/xnum)
    travelt(ig)%num=ig
end do



do iy=1,ynum
    do ig=(iy-1)*xnum+1,iy*xnum
        travelt(ig)%dxx=dxx(iy)
    end do
end do

!---------------------------------------------------------------------


! Determine the center and the range of refine area.
!---------------------------------------------------------------------
rdx=dx/rfxint
rdy=dy/rfyint
rcenter%x=(aint((source%x-minx)/rdx))*rdx+minx
rcenter%y=(aint((source%y-miny)/rdy))*rdy+miny
write(*,*)"Line344 center of refined area",rcenter%x,rcenter%y

if((rcenter%x-rfxrg) .lt. minx)then
    rminx=minx
else
    rminx=rcenter%x-rfxrg
end if
if((rcenter%x+rfxrg) .gt. maxx)then
    rmaxx=maxx
else
    rmaxx=rcenter%x+rfxrg
end if

if((rcenter%y-rfyrg) .lt. miny)then
    rminy=miny
else
    rminy=rcenter%y-rfyrg
end if
if((rcenter%y+rfyrg) .gt. maxy)then
    rmaxy=maxy
else
    rmaxy=rcenter%y+rfyrg
end if

rminxl=nint((rminx-rcenter%x)/rdx)
rmaxxl=nint((rmaxx-rcenter%x)/rdx)
rminyl=nint((rminy-rcenter%y)/rdy)
rmaxyl=nint((rmaxy-rcenter%y)/rdy)

rxnum=rmaxxl-rminxl+1
rynum=rmaxyl-rminyl+1
rttlgnum=rxnum*rynum
!---------------------------------------------------------------------

write(*,*)"line380",rminx,rmaxx,rminy,rmaxy,rminxl,rmaxxl,rminyl,rmaxyl

! Choose coordinate for refined grids
!--------------------------------------------------------------------
if(icoodnt .eq. 1)then
    rdxx=(pi*radii)*rdx/180.0
    rdyy=(pi*radii)*rdy/180.0
else if(icoodnt .eq. 2)then
    rdyy=(pi*radii)*rdy/180.0
    do iy=1,rynum
        rdxx(iy)=cos((rminy+((real(iy)-0.5)*rdy))*pi/180.0)&
        &   *pi*radii*rdx/180.0
    end do
end if


do ig=1,rttlgnum
    rtravelt(ig)%x=rminx+rdx*(ig-((ig-1)/rxnum)*rxnum-1)
    rtravelt(ig)%y=rminy+rdy*((ig-1)/rxnum)
    rtravelt(ig)%num=ig
end do

do iy=1,rynum
    do ig=(iy-1)*rxnum+1,iy*rxnum
        rtravelt(ig)%dxx=rdxx(iy)
    end do
end do
!---------------------------------------------------------------------



! Caculate velocity of refined grids with bilinear interpolation method.
!-----------------------------------------------------------------------
open(230,file='refinevel.txt',status='replace')
write(230,*)rttlgnum
do ig=1,rttlgnum
    call poly_locat(rtravelt(ig)%x,rtravelt(ig)%y,nvtx,vtx,rst)
    if(rst .eq. -1)then
        rv(ig)=voutside
    else if(rst .eq. 1)then
        call baryc_locat(rtravelt(ig)%x,rtravelt(ig)%y,vorig,tri,&
             &trinum,ptri)
        call baryc(rvtemp,1.0/ptri(1)%p%t,1.0/ptri(2)%p%t,1.0/ptri(3)%p%t,&
             &ptri(1)%p%x,ptri(1)%p%y,ptri(2)%p%x,ptri(2)%p%y,&
             &ptri(3)%p%x,ptri(3)%p%y,rtravelt(ig)%x,rtravelt(ig)%y)
        rv(ig)=1.0/rvtemp
    else
        write(*,*)"Error in poly_locat!! rst=",rst
        stop
    end if
    write(230,*)ig,rv(ig),rtravelt(ig)%x,rtravelt(ig)%y
end do
close(230)
!----------------------------------------------------------------------

write(*,*)"line426"

! Initial values around the source.
!----------------------------------------------------------------------
rip=0 ! Counter of living grid node 
ist=-1 ! Counter of updating living grid nodes with same travel time
rcnum=1-rminxl-rminyl*rxnum
write(*,*)"line438",rcnum,rttlgnum,rdx,rdy
tempara1=abs(source%x-rdx*aint(source%x/rdx))
tempara2=abs(source%y-rdy*aint(source%y/rdy))
if(tempara1 .le. 1D-4*rdx)then
    if(tempara2 .le. 1D-4*rdy)then
        rgnum=rcnum
        rtravelt(rgnum)%t=0
        rtravelt(rgnum)%stat=1
        rip=rip+1
        prtravelt(rip)%p=>rtravelt(rgnum)
        ist=ist+1
    else
        ! Calculate 8 points around source
        do iy=0,1
            do ix=-1,1
                rgnum=rcnum+iy*rxnum+ix
                ! Make sure rgnum in the refine number
                if((rgnum .ge. 1) .and. (rgnum .le. rttlgnum))then
                    ! Make sure rgnum and rcnum/rcnum+rxnum in the same row
                    if(abs(rtravelt(rgnum)%y-rtravelt(rcnum+iy*rxnum)%y) .lt. &
                    &rdy/100.0)then
                        rtravelt(rgnum)%t=sqrt((((source%x-rtravelt(rgnum)%x)*&
                        &   rtravelt(rgnum)%dxx/rdx)**2)+&
                        &   (((source%y-rtravelt(rgnum)%y)*rdyy/rdy)**2))/&
                        &   rv(rgnum)
                        rtravelt(rgnum)%stat=1
                        rip=rip+1
                        prtravelt(rip)%p=>rtravelt(rgnum)
                        ist=ist+1
                    end if
                end if
            end do
            rgnum=rcnum+(iy*3-1)*rxnum
            if((rgnum .ge. 1) .and. (rgnum .le. rttlgnum))then
                rtravelt(rgnum)%t=sqrt((((source%x-rtravelt(rgnum)%x)*&
                &   rtravelt(rgnum)%dxx/rdx)**2)+&
                &   (((source%y-rtravelt(rgnum)%y)*rdyy/rdy)**2))/&
                &   rv(rgnum)
                rtravelt(rgnum)%stat=1
                rip=rip+1
                prtravelt(rip)%p=>rtravelt(rgnum)
                ist=ist+1
            end if
        end do
    end if
else if(tempara2 .le. 1D-4*rdy)then
    do ix=0,1
        do iy=-1,1
            rgnum=rcnum+ix+iy*rxnum
            if((rgnum .ge. 1) .and. (rgnum .le. rttlgnum))then
                if(abs(rtravelt(rgnum)%y-rtravelt(rcnum+iy*rxnum)%y) .lt. &
                &rdy/100.0)then
                    rtravelt(rgnum)%t=sqrt((((source%x-rtravelt(rgnum)%x)*&
                    &   rtravelt(rgnum)%dxx/rdx)**2)+&
                    &   (((source%y-rtravelt(rgnum)%y)*rdyy/rdy)**2))/&
                    &   rv(rgnum)
                    rtravelt(rgnum)%stat=1
                    rip=rip+1
                    prtravelt(rip)%p=>rtravelt(rgnum)
                    ist=ist+1
                end if
            end if
        end do
        rgnum=rcnum+(ix*3-1)
        if ((rgnum .ge. 1) .and. (rgnum .le. rttlgnum))then
            if(abs(rtravelt(rgnum)%y-rtravelt(rcnum)%y) .lt. &
            &rdy/100.0)then
                rtravelt(rgnum)%t=sqrt((((source%x-rtravelt(rgnum)%x)*&
                &   rtravelt(rgnum)%dxx/rdx)**2)+&
                &   (((source%y-rtravelt(rgnum)%y)*rdyy/rdy)**2))/&
                &   rv(rgnum)
                rtravelt(rgnum)%stat=1
                rip=rip+1
                prtravelt(rip)%p=>rtravelt(rgnum)
                ist=ist+1
            end if
        end if
    end do
else
    do ix=0,1
        do iy=0,1
            do ig=0,1
                rgnum=rcnum+ix+(ix*2-1)*ig+iy*rxnum+(iy*2-1)*abs(ig-1)*rxnum
                if ((rgnum .ge. 1) .and. (rgnum .le. rttlgnum))then
                    if(abs(rtravelt(rgnum)%y-rtravelt(rcnum+&
                    &iy*rxnum+(iy*2-1)*abs(ig-1)*rxnum)%y) .lt. &
                    &rdy/100.0)then
                        rtravelt(rgnum)%t=sqrt((((source%x-rtravelt(rgnum)%x)*&
                        &   rtravelt(rgnum)%dxx/rdx)**2)+&
                        &   (((source%y-rtravelt(rgnum)%y)*rdyy/rdy)**2))/&
                        &   rv(rgnum)
                        rtravelt(rgnum)%stat=1
                        rip=rip+1
                        prtravelt(rip)%p=>rtravelt(rgnum)
                        ist=ist+1
                    end if
                end if
            end do
            rgnum=rcnum+ix+iy*rxnum
            if ((rgnum .ge. 1) .and. (rgnum .le. rttlgnum))then
                if(abs(rtravelt(rgnum)%y-rtravelt(rcnum+iy*rxnum)%y) .lt. &
                &rdy/100.0)then
                    rtravelt(rgnum)%t=sqrt((((source%x-rtravelt(rgnum)%x)*&
                    &   rtravelt(rgnum)%dxx/rdx)**2)+&
                    &   (((source%y-rtravelt(rgnum)%y)*rdyy/rdy)**2))/&
                    &   rv(rgnum)
                    rtravelt(rgnum)%stat=1
                    rip=rip+1
                    prtravelt(rip)%p=>rtravelt(rgnum)
                    ist=ist+1
                end if
            end if
        end do
    end do
end if
!----------------------------------------------------------------------

! Sort the initial values, with bubble sort method.
!----------------------------------------------------------------------
allocate(ptemp)
if(rip .gt. 1)then
    do ix=1,rip-1
        do iy=ix+1,rip
            if(prtravelt(iy)%p%t .lt. prtravelt(ix)%p%t)then
                ptemp%p=>prtravelt(ix)%p
                prtravelt(ix)%p=>prtravelt(iy)%p
                prtravelt(iy)%p=>ptemp%p
            end if
        end do
    end do
end if
deallocate(ptemp)
!---------------------------------------------------------------------        
write(*,*)"line568",rip
!do ix=1,rip
!    write(*,*)prtravelt(ix)%p
!end do


! Check if the coordinate of outter grid is as same as refined grid. If
! so, assign the value to outter grid.
!---------------------------------------------------------------------
ip=0
do ig=1,rip
    tempara1=prtravelt(ig)%p%x-dx*&
    &   aint((prtravelt(ig)%p%x+(0.01*rdx))/dx)
    tempara2=prtravelt(ig)%p%y-dy*&
    &   aint((prtravelt(ig)%p%y+(0.01*rdy))/dy)
    if((abs(tempara2) .lt. 0.000001) &
    & .and. (abs(tempara1) .lt. 0.000001))then
       gnum=nint((prtravelt(ig)%p%x-minx)/dx)+1+&
       &      nint((prtravelt(ig)%p%y-miny)/dy)*xnum
       travelt(gnum)%t=prtravelt(ig)%p%t
       travelt(gnum)%stat=1
       ip=ip+1
       ptravelt(ip)%p=>travelt(gnum)
    end if
end do
!----------------------------------------------------------------------



! Refined grids caculation begins.
!------------------------------------------------------------------------
rfrg(1)=rcenter%x+rfxrg-rdx/10.0
rfrg(2)=rcenter%x-rfxrg+rdx/10.0
rfrg(3)=rcenter%y+rfxrg-rdy/10.0
rfrg(4)=rcenter%y-rfxrg+rdy/10.0
i=1
j=1
k=1
nb=0
nbtail=0
do while( (prtravelt(rip)%p%x .lt. rfrg(1)) .and.& 
&         (prtravelt(rip)%p%x .gt. rfrg(2)) .and.&
&         (prtravelt(rip)%p%y .lt. rfrg(3)) .and.&
&         (prtravelt(rip)%p%y .gt. rfrg(4)) )

    i=i+1
    ! Update narrow band grid.
    !----------------------------------------------------------------------
    tempnbnum=0
    do iip=rip-ist,rip
        iipnum=prtravelt(iip)%p%num
        if(imethod .eq. 1)then
            call march1(iipnum,rtravelt,rv,rxnum,rttlgnum,rdyy)
        else if(imethod .eq. 2)then
            j=j+1
            call march2(iipnum,rtravelt,rv,rxnum,rttlgnum,rdyy,nbnode,nbnum)
        else
            write(*,*)"Methods parameter is neither 1 or 2!!"
        end if
        
        ! Add narrow band nodes to tempnbnode
        tempnbnode(tempnbnum+1:tempnbnum+nbnum)=nbnode(1:nbnum)
        tempnbnum=tempnbnum+nbnum

    end do
    
    !----------------------------------------------------------------------

    do inb=1,tempnbnum
    if(rtravelt(tempnbnode(inb))%nbstat .eq. 0)then
        nbtail=nbtail+1
        nb(nbtail)=tempnbnode(inb)
        rtravelt(nb(nbtail))%nbstat=nbtail
        rtravelt(nb(nbtail))%stat=0
        call upheap2d(rtravelt,nb,nbtail)
    else if(rtravelt(tempnbnode(inb))%nbstat .gt. 0)then
        rtravelt(tempnbnode(inb))%stat=0
        call updateheap2d(rtravelt,nb,rtravelt(tempnbnode(inb))%nbstat,nbtail)
        end if
    end do

    ! Find alive grid, delete the root of heap and its children with same
    ! values as the root.
    !----------------------------------------------------------------------
    ist=0
    prtravelt(rip+1)%p=>rtravelt(nb(1))
    prtravelt(rip+1)%p%stat=1
    prtravelt(rip+1)%p%nbstat=0
    nb(1)=nb(nbtail)
    rtravelt(nb(1))%nbstat=1
    nbtail=nbtail-1
    call downheap2d(rtravelt,nb,1,nbtail)
    

    ! Insert refined grid value into outter grid.
    tempara1=prtravelt(rip+1)%p%x-dx*&
    &   anint(prtravelt(rip+1)%p%x/dx)
    tempara2=prtravelt(rip+1)%p%y-dy*&
    &   anint(prtravelt(rip+1)%p%y/dy)
    if((abs(tempara1) .lt. 1d-3*rdx) .and.& 
    &  (abs(tempara2) .lt. 1d-3*rdy))then
       gnum=nint((prtravelt(rip+1)%p%x-minx)/dx)+1+&
       &    nint((prtravelt(rip+1)%p%y-miny)/dy)*xnum
       travelt(gnum)%t=prtravelt(rip+1)%p%t
       travelt(gnum)%stat=1
       ip=ip+1
       ptravelt(ip)%p=>travelt(gnum)
    end if

    rip=rip+1

    do while(rtravelt(nb(1))%t .eq. prtravelt(rip)%p%t)
        ist=ist+1
        prtravelt(rip+1)%p=>rtravelt(nb(1))
        prtravelt(rip+1)%p%stat=1
        prtravelt(rip+1)%p%nbstat=0
        nb(1)=nb(nbtail)
        rtravelt(nb(1))%nbstat=1
        nbtail=nbtail-1
        call downheap2d(rtravelt,nb,1,nbtail)
        rip=rip+1
    end do


    if(ist .gt. 0)then
        do iist=1,ist
            iip=rip-iist+1
            prtravelt(iip)%p%stat=1
            tempara1=prtravelt(iip)%p%x-dx*&
            &   anint(prtravelt(iip)%p%x/dx)
            tempara2=prtravelt(iip)%p%y-dy*&
            &   anint(prtravelt(iip)%p%y/dy)
            if((abs(tempara1) .lt. 1d-3*rdx) .and.& 
            &  (abs(tempara2) .lt. 1d-3*rdy))then
                gnum=nint((prtravelt(iip)%p%x-minx)/dx)+1+&
                &    nint((prtravelt(iip)%p%y-miny)/dy)*xnum
                travelt(gnum)%t=prtravelt(iip)%p%t
                travelt(gnum)%stat=1
                ip=ip+1
                ptravelt(ip)%p=>travelt(gnum)
            end if
        end do
    end if

    !----------------------------------------------------------------

end do


! "refine.out" records the information of refined grids which are caculated
! in this step
open(120,file='refine.out',status='replace')
do iist=1,rip
    write(120,*)prtravelt(iist)%p
end do
close(120)

!-------------------------------------------------------------------------


write(*,*)"Refine part ends!"



!####Caculate non-refined ereas.
!------------------------------------------------------------------------

! Velocity of each node
ist=ip-1
open(232,file='vfile.txt',status='replace')
write(232,*)ttlgnum
do iv=1,ttlgnum

    call poly_locat(travelt(iv)%x,travelt(iv)%y,nvtx,vtx,rst)
    if(rst .eq. -1)then
        v(iv)=voutside
    else if(rst .eq. 1)then
        call baryc_locat(travelt(iv)%x,travelt(iv)%y,vorig,tri,&
             &trinum,ptri)
        call baryc(rvtemp,1.0/ptri(1)%p%t,1.0/ptri(2)%p%t,1.0/ptri(3)%p%t,&
             &ptri(1)%p%x,ptri(1)%p%y,ptri(2)%p%x,ptri(2)%p%y,&
             &ptri(3)%p%x,ptri(3)%p%y,travelt(iv)%x,travelt(iv)%y)
        v(iv)=1.0/rvtemp
     else
         write(*,*)"error in poly_locat!! rst=",rst
         stop
     end if

     write(232,*)travelt(iv)%x,travelt(iv)%y,v(iv),iv

end do
close(232)


i=1
j=1
k=1
nb=0
nbtail=0
! Update live, narrow band grid.
!-------------------------------------------------------------------------
do while(ip .lt. ttlgnum)
    i=i+1
    ! Update narrow band grid.
    !-------------------------------------------------------------------------
    tempnbnum=0
    do iip=ip-ist,ip
        iipnum=ptravelt(iip)%p%num
        if(imethod .eq. 1)then
            call march1(iipnum,travelt,v,xnum,ttlgnum,dyy)
        else if(imethod .eq. 2)then
            j=j+1
            call march2(iipnum,travelt,v,xnum,ttlgnum,dyy,nbnode,nbnum)
        else
            write(*,*)"Method parameter is neither 1 or 2!!!"
        end if

        ! Add narrow band nodes to tempnbnode
        tempnbnode(tempnbnum+1:tempnbnum+nbnum)=nbnode(1:nbnum)
        tempnbnum=tempnbnum+nbnum
    end do

    !----------------------------------------------------------------------

    ! Add narrow band nodes to nb
    do inb=1,tempnbnum
        if(travelt(tempnbnode(inb))%nbstat .eq. 0)then
            nbtail=nbtail+1
            nb(nbtail)=tempnbnode(inb)
            travelt(nb(nbtail))%nbstat=nbtail
            travelt(nb(nbtail))%stat=0
            call upheap2d(travelt,nb,nbtail)
        else if(travelt(tempnbnode(inb))%nbstat .gt. 0)then
            travelt(tempnbnode(inb))%stat=0
            call updateheap2d(travelt,nb,travelt(tempnbnode(inb))%nbstat,nbtail)
        end if
    end do
    
    ! Find alive grid
    !----------------------------------------------------------------------
    ist=0
    ptravelt(ip+1)%p=>travelt(nb(1))
    ptravelt(ip+1)%p%stat=1
    ptravelt(ip+1)%p%nbstat=0
    nb(1)=nb(nbtail)
    travelt(nb(1))%nbstat=1
    nbtail=nbtail-1
    call downheap2d(travelt,nb,1,nbtail)
    ip=ip+1
    do while((travelt(nb(1))%t .eq. ptravelt(ip)%p%t) .and. &
            &(nbtail .ge. 1))
        ist=ist+1
        ptravelt(ip+1)%p=>travelt(nb(1))
        ptravelt(ip+1)%p%stat=1
        ptravelt(ip+1)%p%nbstat=0
        nb(1)=nb(nbtail)
        travelt(nb(1))%nbstat=1
        nbtail=nbtail-1
        call downheap2d(travelt,nb,1,nbtail)
        ip=ip+1
    end do


end do
close(233)
!-------------------------------------------------------------------------



! Find the ray path and frechet derivative
!--------------------------------------------------------------------------

write(sourcenum,1000)isr
pathfname="path"//trim(sourcenum)//".txt"
open(101,file=pathfname,status='replace')
do ir=1,nr
    ! Calculate travel time for one receiver
    call localcood(receiver(ir)%x,receiver(ir)%y,travelt,prv,dx,dy,&
    &minx,maxx,miny,maxy,xnum,ynum)
    receiver(ir)%t=bilinear(prv(1)%p%t,prv(2)%p%t,prv(3)%p%t,prv(4)%p%t,&
                  &prv(1)%p%x,prv(1)%p%y,prv(4)%p%x,prv(4)%p%y,&
                  &receiver(ir)%x,receiver(ir)%y)
    
    ! Calculate ray path for one receiver-source. Travel times grid is
    ! regular grid. Bilinear interpolation will be used.
    call raypath(travelt,receiver(ir),source,path,&
    & ipth,dx,dy,minx,maxx,miny,maxy,xnum,ynum)



    ! Velocity grid is irregular grid. Barycentric interpolation will be
    ! used.
    if(gridtype .eq. 1)then
        dxv=vorig(2)%x-vorig(1)%x
        vxnum=nint((maxx-minx)/dxv)+1
        dyv=vorig(vxnum+1)%y-vorig(1)%y
        vynum=nint((maxy-miny)/dyv)+1
        call frech_regular(vorig,path,ipth,fd,dxv,dyv,minx,&
        &   maxx,miny,maxy,vxnum,vynum)
    else if(gridtype .eq. 2)then
        call frech_tri(vorig,path,ipth,fd,tri,trinum)
    else
        write(*,*)"Grid type is neither 1 or 2, please check line 14 in&
        & 'fmm.inp'."
    end if

    
    ! Verify the travel time by frechet derivative
    temptt=0
    do ippth=1,vnum
        write(28,*)fd(ippth),ippth
        temptt=temptt+fd(ippth)/vorig(ippth)%t
    end do
    write(*,*)"temptt by fd",temptt
    write(29,*)temptt,receiver(ir)%x,receiver(ir)%y,source%x,source%y

    write(101,1101)ipth
1101 format('N ',i6)
    do ippth=1,ipth
        write(101,*)path(ippth)
    end do
end do

close(101)

!------------------------------------------------------------------------


! Write travel time in files.
!------------------------------------------------------------------------
write(sourcenum,1000)isr
1000 format(i5.5)
tfilename="source"//trim(sourcenum)//".out"
write(*,*)tfilename
errfname="err"//trim(sourcenum)//".txt"

open(102,file=tfilename,status='replace')
open(103,file='t_refine.out',status='replace')
open(108,file=errfname,status='replace')
do ist=1,ip-1
    if(ptravelt(ist)%p%t .gt. ptravelt(ist+1)%p%t)then
        write(108,*)ist+1,ptravelt(ist+1)%p,&
        &   ptravelt(ist+1)%p%t-ptravelt(ist)%p%t
    end if
end do
do ist=1,ip
    write(102,1001)ptravelt(ist)%p
    write(103,1001)travelt(ist)
1001 format(1x,f16.12,1x,f16.12,1x,f16.12,1x,f16.12,1x,i6.6,1x,i1,1x,i1)
end do
close(102)    
close(103)
close(108)
!-------------------------------------------------------------------------


deallocate(ptravelt)
deallocate(prtravelt)
deallocate(prv)
deallocate(ptri)
end subroutine cacut

!------------------------------------------------------------------------



! This subroutine uses 1st order upwind method to update the travel times
! of narrow band grids.
!-----------------------------------------------------------------------
subroutine march1(num,travelt,v,xnum,ttlgnum,dyy)

use strct
implicit none
type(tstrct) :: travelt(maxgrid)
real(kind=8) :: v(maxgrid)
real(kind=8) :: dyy
real(kind=8) :: temptx,tempty
real(kind=8) :: a,b,c
integer :: num
integer :: xnum
integer :: ttlgnum
integer :: gnum
integer :: iudlr1,iudlr2
integer :: xsol,ysol

xsol=0
ysol=0
do iudlr1=-1,1,2
  do iudlr2=-1,1,2
  ! Update right, left, down and up grid
     gnum=num+(iudlr1-1)*iudlr2/2+&
     &    xnum*((iudlr1+1)*iudlr2/2)
     ! Check if these grids exist.
     if((gnum .gt. 0) .and. (gnum .le. ttlgnum))then
        if(((real(iudlr2)*travelt(gnum)%x) .lt. &
         &  (real(iudlr2)*travelt(num)%x)) .or. &
         & abs(travelt(gnum)%x-travelt(num)%x) .lt.&
         & 0.0001)then
            
            ! Make sure the  grid is not alive
            if(travelt(gnum)%stat .lt. 1)then

                ! Check upwind direction of x coodinate
                if((travelt(gnum+1)%x .gt. travelt(gnum)%x) .and.&
                & (travelt(gnum+1)%stat .eq. 1))then
                    if((travelt(gnum-1)%x .lt. travelt(gnum)%x) .and.&
                    & (travelt(gnum-1)%stat .eq. 1))then
                        temptx=min(travelt(gnum-1)%t,travelt(gnum+1)%t)
                        xsol=1
                    else
                        temptx=travelt(gnum+1)%t
                        xsol=1
                    end if
                else if((travelt(gnum-1)%x .lt. travelt(gnum)%x)&
                & .and. (travelt(gnum-1)%stat .eq. 1))then
                    temptx=travelt(gnum-1)%t
                    xsol=1
                else
                    temptx=0
                    xsol=0
                end if

                ! Check upwind direction of y coodinate
                if(((gnum-xnum) .gt. 0) .and. &
                & (travelt(gnum-xnum)%stat .eq. 1))then
                    if(((gnum+xnum) .le. ttlgnum) .and. &
                    & (travelt(gnum+xnum)%stat .eq. 1))then
                        tempty=min(travelt(gnum-xnum)%t,&
                        &      travelt(gnum+xnum)%t)
                        ysol=1
                    else
                        tempty=travelt(gnum-xnum)%t
                        ysol=1
                    end if
                else if(((gnum+xnum) .le. ttlgnum) .and. &
                &  (travelt(gnum+xnum)%stat .eq. 1))then
                    tempty=travelt(gnum+xnum)%t
                    ysol=1
                else
                    tempty=0
                    ysol=0
                end if

                a=xsol/(travelt(gnum)%dxx**2)+ysol/(dyy**2)
                b=-2.0*((xsol*temptx)/(travelt(gnum)%dxx**2)+&
                & (ysol*tempty)/(dyy**2))
                c=(xsol*((temptx/travelt(gnum)%dxx)**2))+ &
                & (ysol*((tempty/dyy)**2))-1.0/(v(gnum)**2)
                travelt(gnum)%t=(-b+sqrt(b**2-4.0*a*c))/(2*a)
                travelt(gnum)%stat=2
            end if
        end if
     end if
  end do
end do

end subroutine march1
!-----------------------------------------------------------------------


! This subroutine uses 2nd order upwind method to update the travel times
! of narrow band grids.
!-----------------------------------------------------------------------
subroutine march2(num,travelt,v,xnum,ttlgnum,dyy,nbnode,nbnum)

use strct
implicit none
type(tstrct) :: travelt(maxgrid)
real(kind=8) :: v(maxgrid)
real(kind=8) :: dyy
real(kind=8) :: temptx,tempty
real(kind=8) :: temptx1,temptx2
real(kind=8) :: tempty1,tempty2
real(kind=8) :: a,b,c
real(kind=8) :: vx,vy
real(kind=8) :: vavr
integer :: nbnode(maxnbnode)
integer :: nbnum
integer :: num
integer :: xnum
integer :: ttlgnum
integer :: gnum
integer :: iudlr1,iudlr2
integer :: xsol,ysol

nbnum=0
xsol=0
ysol=0
do iudlr1=-1,1,2
  do iudlr2=-1,1,2
  ! Update right, left, down and up grid
     gnum=num+(iudlr1-1)*iudlr2/2+&
     &    xnum*((iudlr1+1)*iudlr2/2)
     ! Check if these grids exist.
     if((gnum .gt. 0) .and. (gnum .le. ttlgnum))then
        if(((real(iudlr2)*travelt(gnum)%x) .lt. &
         &  (real(iudlr2)*travelt(num)%x)) .or. &
         & abs(travelt(gnum)%x-travelt(num)%x) .lt.&
         & 0.0001)then
            
            ! Make sure the  grid is not alive
            if(travelt(gnum)%stat .lt. 1)then

                ! Check upwind direction of x coodinate
                if((travelt(gnum+1)%x .gt. travelt(gnum)%x) .and.&
                & (travelt(gnum+1)%stat .eq. 1))then
                    if((travelt(gnum-1)%x .lt. travelt(gnum)%x) .and.&
                    & (travelt(gnum-1)%stat .eq. 1))then
                        !Check second order
                        if((travelt(gnum+2)%x .gt. travelt(gnum+1)%x) .and.&
                        & (travelt(gnum+2)%stat .eq. 1))then
                            if((travelt(gnum-2)%x .lt. travelt(gnum-1)%x) .and.&
                            & (travelt(gnum-2)%stat .eq. 1))then
                                temptx1=(2.0*travelt(gnum-1)%t-&
                                & 0.5*travelt(gnum-2)%t)
                                temptx2=(2.0*travelt(gnum+1)%t-&
                                & 0.5*travelt(gnum+2)%t)
                                if(temptx1 .le. temptx2)then
                                    temptx=temptx1
                                    vx=(3.0*v(gnum)+4.0*v(gnum-1)+v(gnum-2))/8.0
                                else
                                    temptx=temptx2
                                    vx=(3.0*v(gnum)+4.0*v(gnum+1)+v(gnum+2))/8.0
                                end if
                                xsol=2
                            else
                                temptx1=travelt(gnum-1)%t
                                temptx2=travelt(gnum+1)%t
                                if(temptx1 .le. temptx2)then
                                    temptx=temptx1
                                    vx=(v(gnum)+v(gnum-1))/2.0
                                else
                                    temptx=temptx2
                                    vx=(v(gnum)+v(gnum-1))/2.0
                                end if
                                xsol=1
                            end if
                        else
                            temptx1=travelt(gnum-1)%t
                            temptx2=travelt(gnum+1)%t
                            if(temptx1 .le. temptx2)then
                                temptx=temptx1
                                vx=(v(gnum)+v(gnum-1))/2.0
                            else
                                temptx=temptx2
                                vx=(v(gnum)+v(gnum+1))/2.0
                            end if
                            xsol=1
                        end if
                    else if((travelt(gnum+2)%x .gt. travelt(gnum+1)%x) .and.&
                    & (travelt(gnum+2)%stat .eq. 1))then
                        temptx=2.0*travelt(gnum+1)%t-0.5*travelt(gnum+2)%t
                        vx=(3.0*v(gnum)+4.0*v(gnum+1)+v(gnum+2))/8.0
                        xsol=2
                    else
                        temptx=travelt(gnum+1)%t
                        vx=(v(gnum)+v(gnum+1))/2.0
                        xsol=1
                    end if
                else if((travelt(gnum-1)%x .lt. travelt(gnum)%x)&
                & .and. (travelt(gnum-1)%stat .eq. 1))then
                    if((travelt(gnum-2)%x .lt. travelt(gnum-1)%x) .and.&
                    & (travelt(gnum-2)%stat .eq. 1))then
                        temptx=2.0*travelt(gnum-1)%t-0.5*travelt(gnum-2)%t
                        vx=(3.0*v(gnum)+4.0*v(gnum-1)+v(gnum-2))/8.0
                        xsol=2
                    else
                        temptx=travelt(gnum-1)%t
                        vx=(v(gnum)+v(gnum-1))/2.0
                        xsol=1
                    end if
                else
                    temptx=0
                    vx=v(gnum)
                    xsol=0
                end if

                ! Check upwind direction of y coodinate
                if(((gnum-xnum) .gt. 0) .and. &
                & (travelt(gnum-xnum)%stat .eq. 1))then
                    if(((gnum+xnum) .le. ttlgnum) .and. &
                    & (travelt(gnum+xnum)%stat .eq. 1))then
                        ! Check 2nd order
                        if(((gnum-2*xnum) .gt. 0) .and. &
                        & (travelt(gnum-2*xnum)%stat .eq. 1))then
                            if(((gnum+2*xnum) .le. ttlgnum) .and. &
                            & (travelt(gnum+2*xnum)%stat .eq. 1))then
                                tempty1=(2.0*travelt(gnum-xnum)%t-&
                                & 0.5*travelt(gnum-2*xnum)%t)
                                tempty2=(2.0*travelt(gnum+xnum)%t-&
                                & 0.5*travelt(gnum+2*xnum)%t)
                                if(tempty1 .le. tempty2)then
                                    tempty=tempty1
                                    vy=(3.0*v(gnum)+4.0*v(gnum-xnum)&
                                    &   +v(gnum-2*xnum))/8.0
                                else
                                    tempty=tempty2
                                    vy=(3.0*v(gnum)+4.0*v(gnum+xnum)&
                                    &   +v(gnum+2*xnum))/8.0
                                end if
                                ysol=2
                            else
                                tempty1=travelt(gnum-xnum)%t
                                tempty2=travelt(gnum+xnum)%t
                                if(tempty1 .le. tempty2)then
                                    tempty=tempty1
                                    vy=(v(gnum)+v(gnum-xnum))/2.0
                                else
                                    tempty=tempty2
                                    vy=(v(gnum)+v(gnum+xnum))/2.0
                                end if
                                ysol=1
                            end if
                        else
                            tempty1=travelt(gnum-xnum)%t
                            tempty2=travelt(gnum+xnum)%t
                            if(tempty1 .le. tempty2)then
                                tempty=tempty1
                                vy=(v(gnum)+v(gnum-xnum))/2.0
                            else
                                tempty=tempty2
                                vy=(v(gnum)+v(gnum+xnum))/2.0
                            end if
                            ysol=1
                        end if
                    else if(((gnum-2*xnum) .gt. 0) .and. &
                    & (travelt(gnum-2*xnum)%stat .eq. 1))then
                        tempty=2.0*travelt(gnum-xnum)%t-&
                        & 0.5*travelt(gnum-2*xnum)%t
                        vy=(3.0*v(gnum)+4.0*v(gnum-xnum)+v(gnum-2*xnum))/8.0
                        ysol=2
                    else
                        tempty=travelt(gnum-xnum)%t
                        vy=(v(gnum)+v(gnum-xnum))/2
                        ysol=1
                    end if
                else if(((gnum+xnum) .le. ttlgnum) .and. &
                &  (travelt(gnum+xnum)%stat .eq. 1))then
                    if(((gnum+2*xnum) .le. ttlgnum) .and. &
                    & (travelt(gnum+2*xnum)%stat .eq. 1))then
                        tempty=2.0*travelt(gnum+xnum)%t-&
                        & 0.5*travelt(gnum+2*xnum)%t
                        vy=(3.0*v(gnum)+4.0*v(gnum+xnum)+v(gnum+2*xnum))/8.0
                        ysol=2
                    else
                        tempty=travelt(gnum+xnum)%t
                        vy=(v(gnum)+v(gnum+xnum))/2.0
                        ysol=1
                    end if
                else
                    tempty=0
                    vy=v(gnum)
                    ysol=0
                end if

                if(xsol .gt. 0)then
                    if(ysol .gt. 0)then
                        vavr=(vx+vy)/2.0
                    else
                        vavr=vx
                    end if
                else 
                    if(ysol .gt. 0)then
                        vavr=vy
                    else
                        write(*,*)"xsol,ysol:",xsol,ysol,"error in march2!"
                    end if
                end if

                if(xsol .lt. 2)then
                    if(ysol .lt. 2)then
                        a=xsol/(travelt(gnum)%dxx**2)+ysol/(dyy**2)
                        b=-2.0*((xsol*temptx)/(travelt(gnum)%dxx**2)+(ysol*tempty)/(dyy**2))
                        c=xsol*((temptx/travelt(gnum)%dxx)**2)+ &
                        & ysol*((tempty/dyy)**2)-1.0/(vavr**2)
                    else if(ysol .eq. 2)then
                        a=xsol/(travelt(gnum)%dxx**2)+9.0/(4.0*(dyy**2))
                        b=-2.0*(xsol*temptx)/(travelt(gnum)%dxx**2)-3.0*tempty/(dyy**2)
                        c=xsol*((temptx/travelt(gnum)%dxx)**2)+(tempty/dyy)**2 &
                        & -1.0/(vavr**2)
                    else
                        write(*,*)"ysol:",ysol,"error in march2!"
                    end if
                else if(xsol .eq. 2)then 
                    if(ysol .lt. 2)then
                        a=ysol/(dyy**2)+9.0/(4.0*(travelt(gnum)%dxx**2))
                        b=-2.0*(ysol*tempty)/(dyy**2)-3.0*temptx/(travelt(gnum)%dxx**2)
                        c=ysol*((tempty/dyy)**2)+(temptx/travelt(gnum)%dxx)**2 &
                        & -1.0/(vavr**2)
                    else if(ysol .eq. 2)then
                        a=9.0/(4.0*(dyy**2))+9.0/(4.0*(travelt(gnum)%dxx**2))
                        b=-3.0*(temptx/(travelt(gnum)%dxx**2)+tempty/(dyy**2))
                        c=(temptx/travelt(gnum)%dxx)**2+(tempty/dyy)**2 &
                        & -1.0/(vavr**2)
                    else
                        write(*,*)"ysol:",ysol,"error in march2!"
                    end if
                else
                    write(*,*)"ysol:",ysol,"error in march2!"
                end if
                travelt(gnum)%t=(-b+sqrt(b**2-4.0*a*c))/(2*a)
                travelt(gnum)%stat=2
                nbnum=nbnum+1
                nbnode(nbnum)=gnum
            end if
        end if
     end if
  end do
end do

end subroutine march2
!-----------------------------------------------------------------------


