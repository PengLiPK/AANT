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
