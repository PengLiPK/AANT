!function: bilinear (bilinear interpolation)
!subroutine: bil_gradient (gradient of bilinear interpolation)
!            localcood (locate the center of line segment in
!                       regular grids)
!            raypath (find raypath when giving source and receiver)
!            frech (calculate frechet derivative of regular grids)
!            baryc (barycentric interplation in 2D (triangle))
!            baryc_grad (gradient of barycentric interpolation)
!            baryc_locat (locate the triangle for barycentric
!                         interpolation)
!            poly_locat (determine a point in or out of an polygon)
!            addnodes (add nodes in Delaunay triangles where ray 
!                       weights are high)
!------------------------------------------------------------------
module strct
    implicit none
    integer,parameter :: maxgrid=250000
    integer,parameter :: maxpathnode=10000
    integer,parameter :: maxtri=5000
    integer,parameter :: maxsource=500
    integer,parameter :: maxreceiver=500
    integer,parameter :: maxdata=25000
    integer,parameter :: maxvel=1000
    integer,parameter :: maxnbnode=5000
    real,parameter :: radii=6371.0
    real,parameter :: pi=3.14159
! Source data structure
    type srstrct
        real(kind=8) :: x
        real(kind=8) :: y
    end type

! stat=1,0,-1 means alive, narrowband, far away.
! stat=2 is temporary, to avoid duplication.
    type tstrct
        real(kind=8) :: x
        real(kind=8) :: y
        real(kind=8) :: t
        real(kind=8) :: dxx
        integer :: num
        integer :: stat
    end type
       
    type pstrct
        type(tstrct),pointer :: p
    end type

    type nb_linklist
        type(tstrct),pointer :: p
        type(nb_linklist),pointer :: prev
        type(nb_linklist),pointer :: next
    end type
! Velocity structure
    type vstrct
        real(kind=8) :: x
        real(kind=8) :: y
        real(kind=8) :: vel
    end type
! Receiver structure
    type rcstrct
        real(kind=8) :: x
        real(kind=8) :: y
        real(kind=8) :: t
    end type

! G 1st norm structure
    type gnstrct
        real(kind=8) :: val
        integer :: num
    end type

    type pgnstrct
        type(gnstrct),pointer :: p
    end type
contains

! Delete element in linklist
subroutine nb_del(item)
implicit none
type(nb_linklist),pointer :: item
type(nb_linklist),pointer :: prev,next

prev=>item%prev
next=>item%next
deallocate(item)
if(associated(prev))prev%next=>next
if(associated(next))next%prev=>prev
item=>next

return

end subroutine nb_del

subroutine nb_output(list)
implicit none
type(nb_linklist),pointer :: list,p

if(associated(list%prev) .eqv. .false.)then
    p=>list%next
else
    p=>list
end if

do while(associated(p))
    write(*,*)p%p
    p=>p%next
end do


return

end subroutine nb_output

! Bilinear interpolation for f(x,y)
!---------------------------------------------------------------
real(kind=8) function bilinear(f11,f21,f12,f22,x1,y1,x2,y2,x,y)

implicit none
real(kind=8) :: f11,f21,f12,f22
real(kind=8) :: x1,x2,y1,y2
real(kind=8) :: x,y

if( ((x2-x1)*(y2-y1)) .le. 0 )then
    write(*,*)"error in function bilinear(strct.f90)"&
            &,((x2-x1)*(y2-y1))
end if

bilinear=(f11*(x2-x)*(y2-y) + f21*(x-x1)*(y2-y)&
   &+f12*(x2-x)*(y-y1) + f22*(x-x1)*(y-y1))&
   &/((x2-x1)*(y2-y1))

return
end function bilinear
!-------------------------------------------------------------------


! Gradient of f(x,y) with bilinear interpolation
!-------------------------------------------------------------------
subroutine bil_gradient(gradx,grady,f11,f21,f12,f22,x1,y1,x2,y2,x,y)

implicit none
real(kind=8) :: f11,f21,f12,f22
real(kind=8) :: gradx,grady
real(kind=8) :: x1,x2,y1,y2
real(kind=8) :: x,y

if( ((x2-x1)*(y2-y1)) .le. 0 )then
    write(*,*)"error in subroutine bil_gradient(strct.f90)"&
            &,((x2-x1)*(y2-y1))
end if

gradx=((f21-f11)*(y2-y) + (f22-f12)*(y-y1))&
     &/((x2-x1)*(y2-y1))

grady=((f12-f11)*(x2-x) + (f22-f21)*(x-x1))&
     &/((x2-x1)*(y2-y1))

return

end subroutine bil_gradient
!----------------------------------------------------------------

! Barycentric interpolation in 2D (triangle)
!--------------------------------------------------------------------
subroutine baryc(f,f1,f2,f3,x1,y1,x2,y2,x3,y3,x,y)

implicit none
real(kind=8) :: f
real(kind=8) :: x,y
real(kind=8) :: f1,f2,f3
real(kind=8) :: x1,x2,x3
real(kind=8) :: y1,y2,y3
real(kind=8) :: area
real(kind=8) :: area1,area2,area3

area=abs(x1*y2-x2*y1+x2*y3-x3*y2+x3*y1-x1*y3)

if( area .le. 0 )then
    write(*,*)"area of triangle is 0! error in subroutine &
            &baryc(strct.f90)",area
end if

area1=abs(x*y2-x2*y+x2*y3-x3*y2+x3*y-x*y3)
area2=abs(x1*y-x*y1+x*y3-x3*y+x3*y1-x1*y3)
area3=area-area1-area2

f=(f1*area1+f2*area2+f3*area3)/area

return

end subroutine baryc
!---------------------------------------------------------------------


! Gradient of barycentric interpolation in 2D (triangle)
!--------------------------------------------------------------------
subroutine baryc_grad(gradx,grady,f1,f2,f3,x1,y1,x2,y2,x3,y3,x,y)

implicit none
real(kind=8) :: gradx,grady
real(kind=8) :: x,y
real(kind=8) :: f1,f2,f3
real(kind=8) :: x1,x2,x3
real(kind=8) :: y1,y2,y3
real(kind=8) :: area
real(kind=8) :: area1,area2,area3
real(kind=8) :: gx1,gx2,gx3
real(kind=8) :: gy1,gy2,gy3

area=abs(x1*y2-x2*y1+x2*y3-x3*y2+x3*y1-x1*y3)

if( area .le. 0 )then
    write(*,*)"area of triangle is 0! error in subroutine &
            &baryc_grad(strct.f90)",area
end if

area1=x*y2-x2*y+x2*y3-x3*y2+x3*y-x*y3
if(area1 .ge. 0)then
    gx1=y2-x2*y+x2*y3-x3*y2+x3*y-y3
    gy1=x*y2-x2+x2*y3-x3*y2+x3-x*y3
else
    gx1=-(y2-x2*y+x2*y3-x3*y2+x3*y-y3)
    gy1=-(x*y2-x2+x2*y3-x3*y2+x3-x*y3)
end if

area2=x1*y-x*y1+x*y3-x3*y+x3*y1-x1*y3
if(area1 .ge. 0)then
    gx2=x1*y-y1+y3-x3*y+x3*y1-x1*y3
    gy2=x1-x*y1+x*y3-x3+x3*y1-x1*y3
else
    gx2=-(x1*y-y1+y3-x3*y+x3*y1-x1*y3)
    gy2=-(x1-x*y1+x*y3-x3+x3*y1-x1*y3)
end if

area3=x1*y2-x2*y1+x2*y-x*y2+x*y1-x1*y
if(area1 .ge. 0)then
    gx3=x1*y2-x2*y1+x2*y-y2+y1-x1*y
    gy3=x1*y2-x2*y1+x2-x*y2+x*y1-x1
else
    gx3=-(x1*y2-x2*y1+x2*y-y2+y1-x1*y)
    gy3=-(x1*y2-x2*y1+x2-x*y2+x*y1-x1)
end if

gradx=(f1*gx1+f2*gx2+f3*gx3)/area
grady=(f1*gy1+f2*gy2+f3*gy3)/area

return

end subroutine baryc_grad
!---------------------------------------------------------------------



! This subroutine determine the local coordinate of point (x,y).
!----------------------------------------------------------------------
subroutine localcood(x,y,travelt,prv,dx,dy,&
            &minx,maxx,miny,maxy,xnum,ynum)

implicit none
type(tstrct), target :: travelt(maxgrid)
type(pstrct), pointer :: prv(:)
real(kind=8) :: x,y
real(kind=8) :: minx,maxx
real(kind=8) :: miny,maxy
real(kind=8) :: dx,dy
integer :: pnum
integer :: xnum,ynum
integer :: tmplnumx,tmplnumy

allocate(prv(1:4))
if(x .le. minx)then
    if(y .le. miny)then
        pnum=1
    else if(y .lt. maxy)then
        tmplnumy=int((y-miny)/dy)
        pnum=tmplnumy*xnum+1
    else
        pnum=(ynum-2)*xnum+1
    end if
else if(x .lt. maxx)then
    if(y .le. miny)then
        tmplnumx=int((x-minx)/dx)
        pnum=tmplnumx+1
    else if(y .lt. maxy)then
        tmplnumx=int((x-minx)/dx)
        tmplnumy=int((y-miny)/dy)
        pnum=tmplnumy*xnum+tmplnumx+1
    else
        tmplnumx=int((x-minx)/dx)
        pnum=(ynum-2)*xnum+tmplnumx+1
    end if
else
    if(y .le. miny)then
        pnum=xnum-1
    else if(y .lt. maxy)then
        tmplnumy=int((y-miny)/dy)
        pnum=tmplnumy*xnum+xnum-1
    else
        pnum=(ynum-2)*xnum+xnum-1
    end if
end if
prv(1)%p=>travelt(pnum)
prv(2)%p=>travelt(pnum+1)
prv(3)%p=>travelt(pnum+xnum)
prv(4)%p=>travelt(pnum+xnum+1)

return
end subroutine
!----------------------------------------------------------------------

! This subroutine locates the triangle for barycentric interpolation
!----------------------------------------------------------------------
subroutine baryc_locat(x,y,velnode,tri,trinum,ptri)

implicit none
type(tstrct), target :: velnode(maxgrid)
type(pstrct), pointer :: ptri(:)
real(kind=8) :: x,y
real(kind=8) :: xt(3),yt(3)
real(kind=8) :: a,b,c
integer :: tri(maxtri,3)
integer :: trinum
integer :: it,i

allocate(ptri(1:3))
do i=1,3
    nullify(ptri(i)%p)
end do


do it=1,trinum
    do i=1,3
        xt(i)=velnode(tri(it,i)+1)%x-x
        yt(i)=velnode(tri(it,i)+1)%y-y
    end do
    a=xt(1)*yt(2)-xt(2)*yt(1)
    b=xt(2)*yt(3)-xt(3)*yt(2)
    if( a*b .ge. -1D-10 )then
        c=xt(3)*yt(1)-xt(1)*yt(3)
        if( (b*c .ge. -1D-10) .and. (c*a .ge. -1D-10) )then
            ptri(1)%p=>velnode(tri(it,1)+1)
            ptri(2)%p=>velnode(tri(it,2)+1)
            ptri(3)%p=>velnode(tri(it,3)+1)
            exit
        end if
    end if
end do

if(associated(ptri(1)%p) .eqv. .false.)then
    write(*,*)"point",x,y,"is not in triangles"
    stop
end if

return
end subroutine baryc_locat
!----------------------------------------------------------------------



! This subroutine test a point in or out of an convex polygon
! rst=1, inside; =-1, outside.
!----------------------------------------------------------------------
subroutine poly_locat(x,y,nvtx,vtx,rst)

implicit none
type(srstrct) :: vtx(maxsource)
real(kind=8) :: x,y
real(kind=8), allocatable :: xp(:),yp(:)
real(kind=8), allocatable :: a(:)
real(kind=8) :: b
integer :: nvtx,rst
integer :: it,i

allocate(xp(1:nvtx))
allocate(yp(1:nvtx))
allocate(a(1:nvtx))

do i=1,nvtx
    xp(i)=vtx(i)%x-x
    yp(i)=vtx(i)%y-y
end do

do i=1,nvtx-1
    a(i)=xp(i)*yp(i+1)-xp(i+1)*yp(i)
end do
a(nvtx)=xp(nvtx)*yp(1)-xp(1)*yp(nvtx)

b=a(1)
do i=2,nvtx
    b=b*a(i)
    if(b .le. 0)exit
end do

if(b .le. 0)then
    rst=-1
else
    rst=1
end if

return
end subroutine poly_locat
!----------------------------------------------------------------------



! This subroutine finds ray path with knowing travel time of each grid.
!-----------------------------------------------------------------------
subroutine raypath(travelt,receiver,source,path,&
           &    ipth,dx,dy,minx,maxx,miny,maxy,xnum,ynum)

implicit none
type(srstrct) :: path(maxpathnode)
type(rcstrct) :: receiver
type(srstrct) :: source
type(tstrct), target :: travelt(maxgrid)
type(pstrct), pointer :: prv(:)
real(kind=8) :: dtx,dty
real(kind=8) :: dx,dy
real(kind=8) :: minx,maxx
real(kind=8) :: miny,maxy
integer :: ipth
integer :: ith
integer :: xnum,ynum

allocate(prv(1:4))
ipth=1
path(ipth)%x=receiver%x
path(ipth)%y=receiver%y

do while((abs(path(ipth)%x-source%x) .gt. dx) .or.&
&       (abs(path(ipth)%y-source%y) .gt. dy))
        
! Find the local rectangle for next ray path point caculation.
!-------------------------------------------------------------------
    call localcood(path(ipth)%x,path(ipth)%y,travelt,prv,dx,dy,&
    &minx,maxx,miny,maxy,xnum,ynum)
!-------------------------------------------------------------------
! Local rectangle finding ends.

! Caculate next point of ray path
!-------------------------------------------------------------------

    call bil_gradient(dtx,dty,prv(1)%p%t,prv(2)%p%t,prv(3)%p%t,&
    &prv(4)%p%t,prv(1)%p%x,prv(1)%p%y,prv(4)%p%x,prv(4)%p%y,&
    &path(ipth)%x,path(ipth)%y)

    path(ipth+1)%x=path(ipth)%x-(((dx+dy)/8.0)*dtx/sqrt(dtx**2+dty**2))
    path(ipth+1)%y=path(ipth)%y-(((dx+dy)/8.0)*dty/sqrt(dtx**2+dty**2))
    
    ith=0
    if(path(ipth+1)%x .lt. minx)then
        path(ipth+1)%x=minx
        ith=ith+1
    else if(path(ipth+1)%x .gt. maxx)then
        path(ipth+1)%x=maxx
        ith=ith+1
    else if(path(ipth+1)%y .lt. miny)then
        path(ipth+1)%y=miny
        ith=ith+1
    else if(path(ipth+1)%y .gt. maxy)then
        path(ipth+1)%y=maxy
        ith=ith+1
    end if


    ipth=ipth+1
!------------------------------------------------------------------

end do

path(ipth)%x=source%x
path(ipth)%y=source%y
deallocate(prv)

return

end subroutine raypath
!----------------------------------------------------------------------


! This subroutine finds Frechet Derivative for inverse with regular 
! velocity grid.
!-----------------------------------------------------------------------
subroutine frech_regular(vel,path,ipath,fd,dx,dy,minx,maxx,miny,maxy,&
            &xnum,ynum)

implicit none
type(srstrct) :: path(maxpathnode)
type(tstrct), target :: vel(maxgrid)
type(pstrct), pointer :: prv(:)
real(kind=8) :: fd(maxvel)
real(kind=8) :: dx,dy
real(kind=8) :: minx,maxx
real(kind=8) :: miny,maxy
real(kind=8) :: x,y
real(kind=8) :: length
real(kind=8) :: a
integer :: gridtype
integer :: ipath
integer :: xnum,ynum
integer :: i

allocate(prv(1:4))

fd=0

do i=1,ipath-1
    x=(path(i)%x+path(i+1)%x)/2
    y=(path(i)%y+path(i+1)%y)/2
    length=(sqrt((path(i)%x-path(i+1)%x)**2&
          &+(path(i)%y-path(i+1)%y)**2))*radii&
          &*pi/180

    call localcood(x,y,vel,prv,dx,dy,minx,maxx,miny,maxy,xnum,ynum)

    a=(prv(4)%p%x-prv(1)%p%x)*(prv(4)%p%y-prv(1)%p%y)
    
    fd(prv(1)%p%num)=length*(prv(4)%p%x-x)*(prv(4)%p%y-y)/a&
                    &+fd(prv(1)%p%num)
    fd(prv(2)%p%num)=length*(x-prv(1)%p%x)*(prv(4)%p%y-y)/a&
                    &+fd(prv(2)%p%num)
    fd(prv(3)%p%num)=length*(prv(4)%p%x-x)*(y-prv(1)%p%y)/a&
                    &+fd(prv(3)%p%num)
    fd(prv(4)%p%num)=length*(x-prv(1)%p%x)*(y-prv(1)%p%y)/a&
                    &+fd(prv(4)%p%num)

end do


return
end subroutine frech_regular
!-----------------------------------------------------------------------


! This subroutine finds Frechet Derivative for inverse with Delaunay
! triangle velocity grid.
!-----------------------------------------------------------------------
subroutine frech_tri(vel,path,ipath,fd,tri,trinum)

implicit none
type(srstrct) :: path(maxpathnode)
type(tstrct), target :: vel(maxgrid)
type(pstrct), pointer :: prv(:)
real(kind=8) :: fd(maxvel)
real(kind=8) :: x,y
real(kind=8) :: length
real(kind=8) :: area
real(kind=8) :: area1,area2,area3
integer :: tri(maxtri,3)
integer :: trinum
integer :: ipath
integer :: i

allocate(prv(1:3))

fd=0

do i=1,ipath-1
    x=(path(i)%x+path(i+1)%x)/2
    y=(path(i)%y+path(i+1)%y)/2
    length=(sqrt((path(i)%x-path(i+1)%x)**2&
          &+(path(i)%y-path(i+1)%y)**2))*radii&
          &*pi/180

    call baryc_locat(x,y,vel,tri,trinum,prv)
    area=abs(prv(1)%p%x*prv(2)%p%y - prv(2)%p%x*prv(1)%p%y&
            &+ prv(2)%p%x*prv(3)%p%y - prv(3)%p%x*prv(2)%p%y&
            &+ prv(3)%p%x*prv(1)%p%y - prv(1)%p%x*prv(3)%p%y)
    
    area1=abs(x*prv(2)%p%y - prv(2)%p%x*y&
            &+ prv(2)%p%x*prv(3)%p%y - prv(3)%p%x*prv(2)%p%y&
            &+ prv(3)%p%x*y - x*prv(3)%p%y)

    area2=abs(prv(1)%p%x*y - x*prv(1)%p%y&
            &+ x*prv(3)%p%y - prv(3)%p%x*y&
            &+ prv(3)%p%x*prv(1)%p%y - prv(1)%p%x*prv(3)%p%y)
    
    area3=area-area1-area2
    fd(prv(1)%p%num)=length*area1/area+fd(prv(1)%p%num)
    fd(prv(2)%p%num)=length*area2/area+fd(prv(2)%p%num)
    fd(prv(3)%p%num)=length*area3/area+fd(prv(3)%p%num)

end do


return
end subroutine frech_tri
!-----------------------------------------------------------------


! This subroutine add nodes in Delaunay triangles with high ray weights.
! (With triangles which have at least two high ray weights vertexes.)
! Input: nodenum (nodes's numbers which ray weights are high)
!        num (the total num of 'nodenum')
!        weight (ray weight of velocity nodes)
!        tri (Delaunay triangles), 
!        trinum (numbers of D tri), 
!        velnode (velocity nodes),
! Output: x (x coordinates of new nodes which will be added)
!         y (y coordinates of new nodes which will be added)
!         i (the total number of new nodes which will be added)
!-----------------------------------------------------------------------
subroutine addnodes(x,y,i,nodenum,num,weight,velnode,tri,trinum,fqthr)

implicit none
type(tstrct), target :: velnode(maxgrid)
type(gnstrct) :: weight(maxvel)
real(kind=8) :: x(maxtri),y(maxtri)
real(kind=8) :: weightsum
real(kind=8) :: dist(3)
real(kind=8) :: distmin
real(kind=8) :: fqthr
real(kind=8) :: xtemp,ytemp
integer :: tri(maxtri,3)
integer :: nodenum(maxvel)
integer :: it(3)
integer :: itri,trinum
integer :: ind,num
integer :: i,j,ids


i=0
do itri=1,trinum
    j=0
    it=0
    do ind=1,num
        if(nodenum(ind) .eq. (tri(itri,1)+1))then
            j=j+1
            it(1)=it(1)+1
            if(it(1) .ge. 2)then
                write(*,*)"Error in addnode, nodes's numbers &
                &of high ray weight have duplicate"
                stop
            end if
            if(j .ge. 2)exit
        else if(nodenum(ind) .eq. (tri(itri,2)+1))then
            j=j+1
            it(2)=it(2)+1
            if(it(2) .ge. 2)then
                write(*,*)"Error in addnode, nodes's numbers &
                &of high ray weight have duplicate"
                stop
            end if
            if(j .ge. 2)exit
        else if(nodenum(ind) .eq. (tri(itri,3)+1))then
            j=j+1
            it(3)=it(3)+1
            if(it(3) .ge. 2)then
                write(*,*)"Error in addnode, nodes's numbers &
                &of high ray weight have duplicate"
                stop
            end if
            if(j .ge. 2)exit
        end if
    end do
    if(j .ge. 2)then
        weightsum=weight(tri(itri,1)+1)%val+&
                 &weight(tri(itri,2)+1)%val+&
                 &weight(tri(itri,3)+1)%val
        xtemp=(weight(tri(itri,1)+1)%val*velnode(tri(itri,1)+1)%x&
             &+weight(tri(itri,2)+1)%val*velnode(tri(itri,2)+1)%x&
             &+weight(tri(itri,3)+1)%val*velnode(tri(itri,3)+1)%x&
             &)/weightsum
        ytemp=(weight(tri(itri,1)+1)%val*velnode(tri(itri,1)+1)%y&
             &+weight(tri(itri,2)+1)%val*velnode(tri(itri,2)+1)%y&
             &+weight(tri(itri,3)+1)%val*velnode(tri(itri,3)+1)%y&
             &)/weightsum
        do ids=1,3
            call gcdist(xtemp,ytemp,velnode(tri(itri,ids)+1)%x,&
                       &velnode(tri(itri,ids)+1)%y,dist(ids))
        end do
        distmin=minval(dist(1:3))
        if(distmin .gt. fqthr)then
            i=i+1
            x(i)=xtemp
            y(i)=ytemp
        else
            write(*,*)"distmin,fqthr,xy,tri",distmin,fqthr,xtemp,ytemp,&
                    &velnode(tri(itri,1)+1)%x,velnode(tri(itri,1)+1)%y,&
                    &velnode(tri(itri,2)+1)%x,velnode(tri(itri,2)+1)%y,&
                    &velnode(tri(itri,3)+1)%x,velnode(tri(itri,3)+1)%y
        end if
    end if
end do


return
end subroutine addnodes
!-----------------------------------------------------------------


! Calculate the great-circle distance.
!--------------------------------------------------------------------
subroutine gcdist(lon1,lan1,lon2,lan2,dist)

implicit none
real(kind=8) :: lan1,lan2
real(kind=8) :: lon1,lon2
real(kind=8) :: tmp1,tmp2
real(kind=8) :: theta
real(kind=8) :: dist

tmp1=(sin((lan1-lan2)*pi/360.0))**2
tmp2=((sin((lon1-lon2)*pi/360.0))**2)*cos(lan1*pi/180.0)*cos(lan2*pi/180.0)

theta=2*asin(sqrt(tmp1+tmp2))

dist=radii*theta

end subroutine gcdist
!--------------------------------------------------------------------------


end module strct
