subroutine setaux(mbc,mx,my,mz,xlower,ylower,zlower,dx,dy,dz,maux,aux)

    ! Called at start of computation before calling qinit, and
    ! when AMR is used, also called every time a new grid patch is created.
    ! Use to set auxiliary arrays
    !   aux(1:maux, 1-mbc:mx+mbc, 1-mbc:my+mbc, 1-mbc:mz+mbc).
    ! Note that ghost cell values may need to be set if the aux arrays
    ! are used by the Riemann solver(s).

    implicit none
    integer, intent(in) :: mbc,mx,my,mz,maux
    real(kind=8), intent(in) :: xlower,ylower,zlower,dx,dy,dz
    real(kind=8), intent(out) ::  aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc)

    real(kind=8) :: xcell, ycell, cp, cs, rho_cell, lambda_cell, mu_cell
    real(kind=8) :: xpcorn(4),ypcorn(4),mag
    integer :: i,j,k
    
    REAL (kind=8) :: rho, lambda, mu, rho2, lambda2, mu2, rho3, lambda3, mu3, pipe_inner, pipe_outer
    common /cparam/  rho, lambda, mu, rho2, lambda2, mu2, rho3, lambda3, mu3, pipe_inner, pipe_outer

    ! Loop over all cells
    do k=1-mbc,mz + mbc
      do j=1-mbc,my + mbc 
        ycell = ylower + (j-0.5d0)*dy
        do i=1-mbc,mx + mbc 
          xcell = xlower + (i-0.5d0)*dx
          
          if (xcell .le. pipe_inner) then
            rho_cell = rho
            lambda_cell = lambda
            mu_cell = mu
          else if (xcell .le. pipe_outer) then
            rho_cell = rho2
            lambda_cell = lambda2
            mu_cell = mu2
          else
            rho_cell = rho3
            lambda_cell = lambda3
            mu_cell = mu3
          end if
            
          aux(1,i,j,k) = lambda_cell
          aux(2,i,j,k) = mu_cell

          ! Calculate p and s wave speeds
          cp = dsqrt((lambda_cell + 2.d0*mu_cell)/rho_cell)
          cs = dsqrt(mu_cell/rho_cell)
          aux(3,i,j,k) = cp
          aux(4,i,j,k) = cs

          ! Calculate normal components / area ratios at lower walls

          ! Note that the z-component is ignored here as it is identically mapped to itself
          call mapc2p(xcell - 0.5d0*dx, ycell - 0.5d0*dy, 0.d0, xpcorn(1), ypcorn(1), mag)
          call mapc2p(xcell - 0.5d0*dx, ycell + 0.5d0*dy, 0.d0, xpcorn(2), ypcorn(2), mag)
          call mapc2p(xcell + 0.5d0*dx, ycell + 0.5d0*dy, 0.d0, xpcorn(3), ypcorn(3), mag)
          call mapc2p(xcell + 0.5d0*dx, ycell - 0.5d0*dy, 0.d0, xpcorn(4), ypcorn(4), mag)

          ! lower wall in x-direction
          mag = dsqrt((xpcorn(2) - xpcorn(1))*(xpcorn(2) - xpcorn(1)) + (ypcorn(2) - ypcorn(1))*(ypcorn(2) - ypcorn(1)))
          aux(5,i,j,k) = (ypcorn(2) - ypcorn(1))/mag
          aux(6,i,j,k) = (xpcorn(1) - xpcorn(2))/mag
          aux(7,i,j,k) = mag/dy

          ! lower wall in y-direction
          mag = dsqrt((xpcorn(4) - xpcorn(1))*(xpcorn(4) - xpcorn(1)) + (ypcorn(4) - ypcorn(1))*(ypcorn(4) - ypcorn(1)))
          aux(8,i,j,k) = (ypcorn(1) - ypcorn(4))/mag
          aux(9,i,j,k) = (xpcorn(4) - xpcorn(1))/mag
          aux(10,i,j,k) = mag/dx

          ! lower wall in z-direction (cross-product of diagonals)... also capacity value
          mag = (xpcorn(3) - xpcorn(1))*(ypcorn(2) - ypcorn(4)) - (xpcorn(2) - xpcorn(4))*(ypcorn(3) - ypcorn(1))
          aux(11,i,j,k) = 0.5d0*dabs(mag)/(dx*dy)
          aux(12,i,j,k) = aux(11,i,j,k)

        end do
      end do
    end do

end subroutine setaux
