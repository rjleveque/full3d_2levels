subroutine mapc2p(xc, yc, zc, xp, yp, zp)
    implicit none

    real(kind=8), intent(in) :: xc, yc, zc
    real(kind=8), intent(out) :: xp, yp, zp
    real(kind=8) :: w

      REAL (kind=8) :: rho, lambda, mu, rho2, lambda2, mu2, rho3, lambda3, mu3, pipe_inner, pipe_outer
      common /cparam/  rho, lambda, mu, rho2, lambda2, mu2, rho3, lambda3, mu3, pipe_inner, pipe_outer


      REAL (kind=8) :: t0wall, amplitude, pulse_span, trans_halfdepth, trans_halfwidth
      common /combc/ t0wall, amplitude, pulse_span, trans_halfdepth, trans_halfwidth

    xp = xc
    yp = yc
    zp = zc

    return
end
