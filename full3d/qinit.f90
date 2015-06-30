subroutine qinit(meqn,mbc,mx,my,mz,xlower,ylower,zlower,dx,dy,dz,q,maux,aux)

    ! Set initial conditions for the q array.
    ! Interpolate from reloaded solution if available
    ! Otherwise, start with trivial ICs

    use amr_reload_module, only: rstfile, lfine

    implicit none

    integer, intent(in) :: mbc,mx,my,mz,maux,meqn
    real(kind=8), intent(in) :: xlower,ylower,zlower,dx,dy,dz
    real(kind=8), intent(in) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc)
    real(kind=8), intent(inout) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc)
    real(kind=8) :: xcell, ycell, zcell, xp, yp, zp, r, s
    real(kind=8) :: q_interp(3)

    REAL (kind=8) :: rho, lambda, mu, rho2, lambda2, mu2, rho3, lambda3, mu3, pipe_inner, pipe_outer
    common /cparam/  rho, lambda, mu, rho2, lambda2, mu2, rho3, lambda3, mu3, pipe_inner, pipe_outer

    REAL (kind=8) :: t0wall, amplitude, pulse_span, trans_halfdepth, trans_halfwidth
    common /combc/ t0wall, amplitude, pulse_span, trans_halfdepth, trans_halfwidth



    integer :: i,j,k,m


!   In order to debug, set initial conditions from analytic expression instead
!   of interpolating 2d results
    if (1 .eq. 2) then
       do k=1,mz
          zcell = zlower + (k-0.5d0)*dz
          do j=1,my
             ycell = ylower + (j-0.5d0)*dy
             do i=1,mx
                xcell = xlower + (i-0.5d0)*dx
                call mapc2p(xcell, ycell, zcell, xp, yp, zp)
                r = dsqrt(yp*yp + zp*zp)
                call interp_reload(1, lfine, 3, xp, r, q_interp)
                q(1,i,j,k) = -q_interp(1)
                q(2,i,j,k) = -q_interp(1)
                q(3,i,j,k) = -q_interp(1)
                q(4,i,j,k) = 0.d0
                q(5,i,j,k) = 0.d0
                q(6,i,j,k) = 0.d0
                q(7,i,j,k) = q_interp(2)
                q(8,i,j,k) = yp / (r + 1.e-10) * q_interp(3)
                q(9,i,j,k) = zp / (r + 1.e-10) * q_interp(3)
             end do
          end do
       end do

    else
       s = dsqrt((lambda + 2.d0*mu)/rho)
       do k=1,mz
          do j=1,my
             do i=1,mx
                xcell = xlower + (i-0.5d0)*dx
                q(1,i,j,k) = -amplitude*dexp(-((48.d0 - (xcell-trans_halfdepth)/s - 3.d0*t0wall)/t0wall)**2)
                q(2,i,j,k) = q(1,i,j,k)
                q(3,i,j,k) = q(1,i,j,k)
                q(4,i,j,k) = 0.d0
                q(5,i,j,k) = 0.d0
                q(6,i,j,k) = 0.d0
                q(7,i,j,k) = -q(1,i,j,k)/rho/s
                q(8,i,j,k) = 0.d0
                q(9,i,j,k) = 0.d0
             end do
          end do
       end do
   end if

       return
       end
