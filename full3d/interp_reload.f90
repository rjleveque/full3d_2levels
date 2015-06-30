subroutine interp_reload(lst, lend, nvar, x, y, q)

    use amr_reload_module
    implicit none
    integer, intent(in) :: lst, lend, nvar
    real(kind=8), intent(in) :: x, y
    real(kind=8), intent(inout) :: q(nvar)
    integer :: nx,ny,loc,locaux,ii,jj,level,mptr,mitot,mjtot
    integer :: ivar,iaux,i,j, iadd
    real(kind=8) :: xlow, ylow, xhi, yhi, dx, dy,wx,wy, xm,ym
    real(kind=8) :: xp, yp

    iadd(ivar,i,j)  = loc + ivar - 1 + nvar*((j-1)*mitot+i-1)

    level = lend
    do while (level .ge. lst)
        mptr = lstart(level)
        do while (mptr .ne. 0)
            nx      = node(ndihi,mptr) - node(ndilo,mptr) + 1
            ny      = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
            xlow = rnode(cornxlo,mptr)
            ylow = rnode(cornylo,mptr)
            dx = hxposs(level)
            dy = hyposs(level)
            loc     = node(store1, mptr)
            locaux  = node(storeaux,mptr)
            xhi = xlow + nx*dx
            yhi = ylow + ny*dy
            mitot   = nx + 2*nghost
            mjtot   = ny + 2*nghost
            xm = xlow - (nghost+0.5d0)*dx
            xp = xhi  + (nghost+0.5d0)*dx
            ym = ylow - (nghost+0.5d0)*dy
            yp = yhi  + (nghost+0.5d0)*dy
            if ((x .ge. xm) .and. (x .le. xp) .and. (y .ge. ym) .and.  &
                (y .le. yp)) then
                ii = int((x-xm)/dx)
                jj = int((y-ym)/dy)
                wx = x - (xm + ii*dx)
                wy = y - (ym + jj*dy)
                do ivar=1,nvar
                    q(ivar) =   &
                        (1.d0-wx)*(1.d0-wy)*alloc(iadd(ivar,ii,jj)) &
                      + (1.d0-wx)*wy*alloc(iadd(ivar,ii,jj+1))  &
                      + wx*(1.d0-wy)*alloc(iadd(ivar,ii+1,jj))  &
                      + wx*wy*alloc(iadd(ivar,ii+1,jj+1)) 
                end do

                return
            end if
            mptr = node(levelptr, mptr)
        end do

        level = level - 1
    end do


    write(6,*) 'Data not found in reload grids'
    write(6,*) '(x,y)=', x, y
    write(6,*) '(xlower,xupper)', xlower, xupper
    write(6,*) '(ylower,yupper)', ylower, yupper
    stop
      
end subroutine interp_reload
