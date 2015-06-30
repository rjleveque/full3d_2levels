! -------------------------------------------------------------------
      subroutine flag2refine(mx,my,mz,mbc,meqn,maux,xlower,ylower, &
                        zlower,dx,dy,dz,t, &
                        level,tolsp,q,aux,amrflags,DONTFLAG,DOFLAG)
! -------------------------------------------------------------------

! ::::::::::::::::::::: flag2refine ::::::::::::::::::::::::::::::::::
!
! User routine to control flagging of points for refinement.
!
! Default version computes spatial difference dq in each direction and
! for each component of q and flags any point where this is greater than
! the tolerance tolsp.  This is consistent with what the routine errsp did in
! earlier versions of amrclaw (4.2 and before).
!
! This routine can be copied to an application directory and modified to
! implement some other desired refinement criterion.
!
! The logical function allowflag(x,y,t,level) is called to check whether
! further refinement at this level is allowed at this particular location
! and time.  The default library version of this routine returns .true.
! for all arguments.  Copy that routine to the application directory and
! modify it if needed to restrict the region where refinement is allowed.
!
! Points may also be flagged for refining based on a Richardson estimate
! of the error, obtained by comparing solutions on the current grid and a
! coarsened grid.  Points are flagged if the estimated error is larger than
! the parameter tol in amr2ez.data, provided tol>0.  If tol<=0 then
! the coarsening and Richardson estimation is not performed!
! This is a change from previous versions (4.2 and before) of amrclaw.
! Note: in previous versions, the routine errf1 used a function
! allowed(x,y,level) that has been replaced by the allowflag.  This new
! function is also used in Richardson estimation if that is invoked.
!
!
!    q   = grid values including ghost cells (bndry vals at specified
!          time have already been set, so can use ghost cell values too)
!
!  aux   = aux array on this grid patch
!
! amrflags  = array to be flagged with either the value
!             DONTFLAG (no refinement needed)  or
!             DOFLAG   (refinement desired)
!
! tolsp = tolerance specified by user in input file amr2ez.data, used in default
!         version of this routine as a tolerance for spatial differences.

! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      implicit none


      real (kind=8) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc)
      real (kind=8) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc)
      real (kind=8) :: amrflags(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc)
      logical     allowflag
      external    allowflag
      real (kind=8) :: DOFLAG, DONTFLAG, xlower, ylower, zlower, dx, dy, dz, tolsp, t
      integer :: mx, my, mz, mbc, meqn, maux

      integer :: ccount, level, nx, i, j, k, l
      real (kind=8) :: x, y, z, tv, wavelength, cpipe, cinner

      REAL (kind=8) :: rho, lambda, mu, rho2, lambda2, mu2, rho3, lambda3, mu3, pipe_inner, pipe_outer
      common /cparam/  rho, lambda, mu, rho2, lambda2, mu2, rho3, lambda3, mu3, pipe_inner, pipe_outer


!     # loop over interior points on this grid:

      cpipe = dsqrt((lambda2 + 2.d0*mu2)/rho2)
      cinner = dsqrt((lambda + 2.d0*mu)/rho)
      wavelength = 2.d0*(pipe_outer - pipe_inner)/cpipe*cinner
      if (wavelength/dx .ge. 8.d0) then
        nx = (0.125d0*wavelength)/dx
    else
        nx = 1
    end if

      do k = 1,mz
      z = zlower + (k-0.5d0)*dz

      do j = 1,my
        y = ylower + (j-0.5d0)*dy
        do i = 1,mx
          x = xlower + (i-0.5d0)*dx

          amrflags(i,j,k) = DONTFLAG

          if (allowflag(x,y,z,t,level)) then
!            # check to see if we should flag this point for refinement.
!            # Here the default test is taken from errsp.f in previous
!            # versions of amrclaw -- flag this point if dq > tolsp:
                 ccount = 0
                 tv = 0.d0
                 do l = i-nx+1,i+nx
                     if (l .ge. 1 .and. l .le. mx) then
                         ccount = ccount + 1

                         tv = tv + dabs(q(7,l,j,k) - q(7,l-1,j,k))
                     end if
                 end do
                 if (tv .ge. ccount/(2.d0*nx)*tolsp) then
                     amrflags(i,j,k) = DOFLAG
                 end if
          end if

        end do
    end do
end do

      return
      end
