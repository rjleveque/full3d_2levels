
!     =========================================
      logical function allowflag(x,y,z,t,level)
!     =========================================
!     # Indicate whether the grid point at (x,y,z,t) at this refinement level
!     # is allowed to be flagged for further refinement.
!
!     # This is useful if you wish to zoom in on some structure in a 
!     # known location but don't want the same level of refinement elsewhere.  
!     # Points are flagged only if one of the errors is greater than the 
!     # corresponding tolerance.
!
!     # For example, to allow refinement of Level 1 grids everywhere but
!     # of finer grids only for  y >= 0.4:
!     # allowed(x,y,z,t,level) = (level.le.1 .or. y.ge.0.4d0) 
!
!     # This routine is called from routine flag2refine.
!     # If Richardson error estimates are used (if tol>0) then this routine
!     # is also called from errf1.

      implicit none

      REAL (kind=8) :: t0wall, amplitude, pulse_span, trans_halfdepth, trans_halfwidth
      common /combc/ t0wall, amplitude, pulse_span, trans_halfdepth, trans_halfwidth

      real (kind=8) :: x, y, z, t
      integer level

      real (kind=8) xp, yp, zp

      call mapc2p(x,y,z,xp,yp,zp)

      !if ((dsqrt(yp*yp + zp*zp) .le. trans_halfwidth) .or. (t.le.50.d0)) then
      if ((dsqrt(yp*yp + zp*zp) .le. 0.015d0) .or. (t.le.48.5d0)) then
            allowflag = .true.
      else
            allowflag = .false.
      end if

      return
      end
