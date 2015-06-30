subroutine setprob

      implicit none
      REAL (kind=8) :: rho, lambda, mu, rho2, lambda2, mu2, rho3, lambda3, mu3, pipe_inner, pipe_outer
      common /cparam/  rho, lambda, mu, rho2, lambda2, mu2, rho3, lambda3, mu3, pipe_inner, pipe_outer
      
      REAL (kind=8) :: t0wall, amplitude, pulse_span, trans_halfdepth, trans_halfwidth
      common /combc/ t0wall, amplitude, pulse_span, trans_halfdepth, trans_halfwidth

      character(len=200) :: restart_directory, checkpoint_file
!
     ! local variables:
     character*12 :: fname
     integer :: iunit, i, nrst
!
      iunit = 7
      fname = 'setprob.data'
!     # open the unit with new routine from Clawpack 4.4 to skip over
!     # comment lines starting with #:
      call opendatafile(iunit, fname)

!
       read(7,*) restart_directory
       read(7,*) checkpoint_file
       read(7,*) rho
       read(7,*) lambda
       read(7,*) mu
       read(7,*) rho2
       read(7,*) lambda2
       read(7,*) mu2
       read(7,*) rho3
       read(7,*) lambda3
       read(7,*) mu3
       read(7,*) t0wall
       read(7,*) amplitude
       read(7,*) pulse_span
       read(7,*) trans_halfdepth
       read(7,*) trans_halfwidth
       read(7,*) pipe_inner
       read(7,*) pipe_outer

       ! If this run is initialized with 2D axisymmetric data, reload that solution
       if (restart_directory .ne. '') then
            ! need to count whitespace in filenames for truncation in reload
            nrst = 0
            do i=1,200
                if (restart_directory(i:i) .ne. ' ') then
                    nrst = nrst + 1
                else
                    exit
                end if
            end do
            call amr2_reload(restart_directory, checkpoint_file, nrst)
       end if

      return
      end
