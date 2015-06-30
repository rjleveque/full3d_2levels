subroutine amr2_reload(restart_directory, checkpoint_file, nrst)

    use amr_reload_module
    implicit none

    integer :: nsteps,nvar,naux,nrst
    real(kind=8) :: time
    ! Local variables
    integer :: i, iaux, mw, level
    integer :: ndim, mcapa1, mindim
    integer :: nstart, nv1, nx, ny, lentotsave, num_gauge_SAVE
    integer :: omp_get_max_threads, maxthreads
    real(kind=8) :: ratmet, cut, dtinit, dt_max
    logical :: vtime, rest, output_t0 
    real(kind=8) :: dxmin, dymin
    real(kind=8) :: x,y
    real(kind=8), allocatable :: q(:)
    integer :: k

    character(len=200) :: restart_directory, checkpoint_file, rstfile_ignore
    character(len=364) :: format_string, clawfile, amrfile, outfile, dbugfile, matfile, parmfile
    character(len=nrst) :: truncated_restart_directory

    truncated_restart_directory = restart_directory(1:nrst)
    clawfile = truncated_restart_directory // 'claw.data'
    amrfile = truncated_restart_directory // 'amr.data'
    outfile = truncated_restart_directory // 'fort.amr'
    dbugfile = truncated_restart_directory // 'fort.debug'
    matfile = truncated_restart_directory // 'fort.nplot'
    parmfile = truncated_restart_directory // 'fort.parameters'
    rstfile = truncated_restart_directory // checkpoint_file

    ! Read in claw.data and amr.data files to get parameters
    ! that were used in the axisymmetric run.  Many of these
    ! parameters are irrelevant to the restart.

    ! Open AMRClaw primary parameter file
    call opendatafile(inunit,clawfile)
   ! Number of space dimensions, not really a parameter but we read it in and
    ! check to make sure everyone is on the same page. 
    read(inunit,"(i1)") ndim  
    if (ndim /= 2) then
        print *,'Error ***   ndim = 2 is required,  ndim = ',ndim
        print *,'*** Are you sure input has been converted'
        print *,'*** to Clawpack 5.x form?'
        stop
    endif
          
    ! Domain variables
    read(inunit,*) xlower, ylower
    read(inunit,*) xupper, yupper
    read(inunit,*) nx, ny
    read(inunit,*) nvar    ! meqn
    read(inunit,*) mwaves
    read(inunit,*) naux
    read(inunit,*) t0

    ! ==========================================================================
    ! Output Options
    ! Output style
    read(inunit,*) output_style
    if (output_style == 1) then
        read(inunit,*) nout
        read(inunit,*) tfinal
        read(inunit,*) output_t0

        iout = 0
    else if (output_style == 2) then
        read(inunit,*) nout
        allocate(tout(nout))
        read(inunit,*) (tout(i), i=1,nout)
        output_t0 = (tout(1) == t0)
        ! Move output times down one index
        if (output_t0) then
            nout = nout - 1
            do i=1,nout
                tout(i) = tout(i+1)
            enddo
        endif
        iout = 0
        tfinal = tout(nout)
    else if (output_style == 3) then
        read(inunit,*) iout
        read(inunit,*) nstop
        read(inunit,*) output_t0
        nout = 0
        tfinal = rinfinity
    else
        stop "Error ***   Invalid output style."
    endif

    ! Error checking
    if ((output_style == 1) .and. (nout > 0)) then
        allocate(tout(nout))
        do i=1,nout
            tout(i) = t0 + i * (tfinal - t0) / real(nout,kind=8)
        enddo
    endif

    ! What and how to output
    read(inunit,*) output_format
    allocate(output_q_components(nvar))
    read(inunit,*) (output_q_components(i),i=1,nvar)
    if (naux > 0) then
        allocate(output_aux_components(naux))
        read(inunit,*) (output_aux_components(i),i=1,naux)
        read(inunit,*) output_aux_onlyonce
    endif
    ! ==========================================================================

    ! ==========================================================================
    !  Algorithm parameters

    read(inunit,*) possk(1)   ! dt_initial
    read(inunit,*) dt_max     ! largest allowable dt
    read(inunit,*) cflv1      ! cfl_max
    read(inunit,*) cfl        ! clf_desired
    read(inunit,*) nv1        ! steps_max
      
    if (output_style /= 3) then
        nstop = nv1
    endif

    read(inunit,*) vtime      ! dt_variable
    if (vtime) then
        method(1) = 2
    else
        method(1) = 1
    endif

    read(inunit,*) method(2)  ! order
    iorder = method(2)
    read(inunit,*) method(3)  ! order_trans

    read(inunit,*) dimensional_split
       if (dimensional_split > 1) then
           print *, '*** ERROR ***  dimensional_split = ', dimensional_split
           print *, ' Strang splitting not supported in amrclaw'
           stop
       endif

    read(inunit,*) method(4)   ! verbosity
    read(inunit,*) method(5)   ! src_split
    read(inunit,*) mcapa1
    
    read(inunit,*) use_fwaves
    allocate(mthlim(mwaves))
    read(inunit,*) (mthlim(mw), mw=1,mwaves)

    ! Boundary conditions
    read(inunit,*) nghost
    read(inunit,*) mthbc(1),mthbc(3)
    read(inunit,*) mthbc(2),mthbc(4)

    ! 1 = left, 2 = right 3 = bottom 4 = top boundary
    xperdom = (mthbc(1) == 2 .and. mthbc(2) == 2)
    yperdom =  (mthbc(3) == 2 .and. mthbc(4) == 2)
    spheredom =  (mthbc(3) == 5 .and. mthbc(4) == 5)

    if ((mthbc(1).eq.2 .and. mthbc(2).ne.2) .or. &
        (mthbc(2).eq.2 .and. mthbc(1).ne.2)) then
        
        print *, '*** ERROR ***  periodic boundary conditions: '
        print *, '  mthbc(1) and mthbc(2) must BOTH be set to 2'
        stop
    endif

    if ((mthbc(3).eq.2 .and. mthbc(4).ne.2) .or. &
        (mthbc(4).eq.2 .and. mthbc(3).ne.2)) then

        print *, '*** ERROR ***  periodic boundary conditions: '
        print *, '  mthbc(3) and mthbc(4) must BOTH be set to 2'
        stop
    endif

    if ((mthbc(3).eq.5 .and. mthbc(4).ne.5) .or. &
        (mthbc(4).eq.5 .and. mthbc(3).ne.5)) then
    
        print *, '*** ERROR ***  sphere bcs at top and bottom: '
        print *, '  mthbc(3) and mthbc(4) must BOTH be set to 5'
        stop
    endif

    if (spheredom .and. .not. xperdom) then

        print *,'*** ERROR ***  sphere bcs at top and bottom: '
        print *,'must go with periodic bcs at left and right  '
        stop
    endif

    ! ==========================================================================
    !  Restart and Checkpointing

    read(inunit,*) rest
    read(inunit,*) rstfile_ignore  ! since we set rstfile above


    read(inunit,*) checkpt_style
    if (checkpt_style == 0) then
        ! Never checkpoint:
        checkpt_interval = iinfinity

    else if (checkpt_style == 2) then
        read(inunit,*) nchkpt
        allocate(tchk(nchkpt))
        read(inunit,*) (tchk(i), i=1,nchkpt)

    else if (checkpt_style == 3) then
        ! Checkpoint every checkpt_interval steps on coarse grid
        read(inunit,*) checkpt_interval
    endif

    close(inunit)

    ! ==========================================================================
    !  Refinement Control
    call opendatafile(inunit, amrfile)

    read(inunit,*) mxnest
    if (mxnest <= 0) then
        stop 'Error ***   mxnest (amrlevels_max) <= 0 not allowed'
    endif
          
    if (mxnest > maxlv) then
        stop 'Error ***   mxnest > max. allowable levels (maxlv) in common'
    endif
      
    ! Anisotropic refinement always allowed in 5.x:
    read(inunit,*) (intratx(i),i=1,max(1,mxnest-1))
    read(inunit,*) (intraty(i),i=1,max(1,mxnest-1))
    read(inunit,*) (kratio(i), i=1,max(1,mxnest-1))
    read(inunit,*)

    do i=1,mxnest-1
        if ((intratx(i) > max1d) .or. (intraty(i) > max1d)) then 
            print *, ""
            format_string = "(' *** Error: Refinement ratios must be no " // &
                            "larger than max1d = ',i5,/,'     (set max1d" // &
                            " in amr_module.f90)')"
            print format_string, max1d
            stop
        endif
    enddo

    if (naux > 0) then
        allocate(auxtype(naux))
        read(inunit,*) (auxtype(iaux), iaux=1,naux)
    endif
    read(inunit,*)
              
    read(inunit,*) flag_richardson
    read(inunit,*) tol            ! for richardson
    read(inunit,*) flag_gradient
    read(inunit,*) tolsp          ! for gradient
    read(inunit,*) kcheck
    read(inunit,*) ibuff
    read(inunit,*) cut
    read(inunit,*) verbosity_regrid

    ! read verbose/debugging flags
    read(inunit,*) dprint
    read(inunit,*) eprint
    read(inunit,*) edebug
    read(inunit,*) gprint
    read(inunit,*) nprint
    read(inunit,*) pprint
    read(inunit,*) rprint
    read(inunit,*) sprint
    read(inunit,*) tprint
    read(inunit,*) uprint

    close(inunit)
    ! Finished with reading in parameters

    call reload(nsteps,time,nvar,naux)

    write(6,*) 'lfine = ',lfine
    write(6,*) 'time = ',time
    allocate(output_aux_components(naux))

    output_format = 1
    output_aux_components = 0

end subroutine amr2_reload
