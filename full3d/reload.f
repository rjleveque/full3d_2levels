
c
c ---------------------------------------------------------
c
      subroutine reload(nsteps,time,nvar,naux)

c
      use amr_reload_module
      implicit double precision (a-h,o-z)
 
 
      logical foundFile
      dimension intrtx(maxlv),intrty(maxlv),intrtt(maxlv)
c
c :::::::::::::::::::::::::::: RELOAD ::::::::::::::::::::::::::::::::
c read back in the check point files written by subr. check.
c
c Modified version of restrt that only loads old data
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
c     !! Now allow user-specified file name !!
c     rstfile  = 'restart.data'

      write(6,*) 'Attempting to reload data '
      write(6,*) '  checkpoint file: ',trim(rstfile)
      inquire(file=trim(rstfile),exist=foundFile)
      if (.not. foundFile) then
        write(*,*)" Did not find checkpoint file!"
        stop
      endif
      open(rstunit,file=trim(rstfile),status='old',form='unformatted')
      rewind rstunit

      read(rstunit) lenmax,lendim,isize

c     # need to allocate for dynamic memory:
c     call restrt_alloc(isize)
      allocate(alloc(isize))
      memsize = isize

      read(rstunit) (alloc(i),i=1,lendim)
      read(rstunit) hxposs,hyposs,possk,icheck
      read(rstunit) lfree,lenf
      read(rstunit) rnode,node,lstart,newstl,listsp,tl,
     1       ibuf,mstart,ndfree,lfine,iorder,mxnest,
     2       intratx,intraty,kratio,iregsz,jregsz,
     2       iregst,jregst,iregend,jregend,
     3       numgrids,kcheck1,nsteps,time,
     3       matlabu
      read(rstunit) avenumgrids, iregridcount,
     1              evol,rvol,rvoll,lentot,tmass0,cflmax

      close(rstunit) 

      write(outunit,100) nsteps,time
      write(6,100) nsteps,time
 100  format(/,' Data comes from calculating over ',i5,' steps',
     1        /,'  (time = ',e15.7,')')
c
c     skip error checking that refinement ratios have not changed

      do i = 1, min(mxnold-1,mxnest-1)
        if ( (intratx(i) .ne. intrtx(i)) .or.
     .       (intraty(i) .ne. intrty(i)) ) then
c    .       (kratio(i) .ne.  intrtt(i) .and. .not. varRefTime) ) then
        write(outunit,*) 
     .  " not allowed to change existing refinement ratios on Restart"
        write(outunit,*)" Old ratios:"
        write(*,*)      " Old ratios:"
        write(outunit,903)(intrtx(j),j=1,mxnold-1)
        write(*,903)      (intrtx(j),j=1,mxnold-1)
        write(outunit,903)(intrty(j),j=1,mxnold-1)
        write(*,903)      (intrty(j),j=1,mxnold-1)
c       write(outunit,903)(intrtt(j),j=1,mxnold-1)
c       write(*,903)      (intrtt(j),j=1,mxnold-1)
 903    format(6i3)
        stop
       endif
      end do

c     if (varRefTime) then  ! reset intrat to previously saved ratios, not input ratios
      if (.true.) then  ! reset intrat to previously saved ratios, not input ratios
        do i = 1, mxnold-1
            kratio(i) = intrtt(i)
        end do
      endif

c
c adjust free list of storage in case size has changed.
c
      idif = memsize - isize
      print *,'+++', isize,memsize
      if (idif .gt. 0) then
          lfree(lenf,1) = isize + 2
          call reclam(isize+1,idif)
      else if (idif .lt. 0) then
            write(outunit,900) isize, memsize
            write(*,900)       isize, memsize
 900        format(' size of alloc not allowed to shrink with ',/,
     .             ' restart old size ',i7,' current size  ',i7)
            stop
      endif
c
      return
      end
