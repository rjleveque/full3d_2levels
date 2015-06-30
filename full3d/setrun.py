""" 
Module to set up run time parameters for Clawpack -- classic code.

The values set in the function setrun are then written out to data files
that will be read in by the Fortran code.
    
""" 

import os
import numpy as np

# Used for obtain Z_target when code is called by run_imp_sweep
from clawpack.clawutil.data import ClawData
#imp_sweep_data = ClawData()
#imp_sweep_data.read('imp_sweep.data', force=True)
#Z_target = imp_sweep_data.Z_target
#Z_target = 8.0e6

# Flag to use only one octant with symmetry
use_symmetry = 1

# Physical parameters used at some point
trans_halfwidth = 0.025
domain_radius = 0.215
pipe_outer_radius = 0.1916 
pipe_inner_radius = 0.1866
trans_halfdepth = 0.1162 #how far the transducer face is from the center of the pipe

# Adjust the Concrete material parameters to obtain target impedence
#rho0 = 1970
#lambda0 = 11075734000.0
#mu0 = 7801397000.0
#c0 = ((lambda0 + 2.0*mu0)/rho0)**(0.5)
#Z0 = c0*rho0
#scaling_factor = Z_target/Z0
#rho3 = scaling_factor**(0.5)*rho0
#lambda3 = scaling_factor**(1.5)*lambda0
#mu3 = scaling_factor**(1.5)*mu0
 
#------------------------------
def setrun(claw_pkg='amrclaw'):
#------------------------------
    
    """ 
    Define the parameters used for running Clawpack.

    INPUT:
        claw_pkg expected to be "amrclaw" for this setrun.

    OUTPUT:
        rundata - object of class ClawRunData 
    
    """ 
    
    from clawpack.clawutil import data 
    
    
    assert claw_pkg.lower() == 'amrclaw',  "Expected claw_pkg = 'amrclaw'"

    num_dim = 3
    rundata = data.ClawRunData(claw_pkg, num_dim)

    #------------------------------------------------------------------
    # Problem-specific parameters to be written to setprob.data:
    #------------------------------------------------------------------
    # Sample setup to write one line to setprob.data ...
    probdata = rundata.new_UserData(name='probdata',fname='setprob.data')
    probdata.add_param('restart_directory',  '', 'directory of restart files (empty string for NONE)')
#    probdata.add_param('restart_directory',  '../../acoustics/_output/', 'directory of restart files (empty string for NONE)')
    probdata.add_param('checkpoint_file', 'fort.chk00072', 'checkpoint file name (emptry string for NONE)')
    # Water inside pipe
    probdata.add_param('rho',        1000.0,  'density') #kg/m^3 ==> Speeds in m/us
    probdata.add_param('lambda',     2.202256e-3,  'Lame parameter lambda')
    probdata.add_param('mu',         0.0,  'Lame parameter mu')
    # Steel pipe
    probdata.add_param('rho2',       7850.0,  'density') #kg/m^3 ==> Speeds in m/us
    probdata.add_param('lambda2',    1.145405275e-1,  'Lame parameter lambda')
    probdata.add_param('mu2',        8.215201625e-2,  'Lame parameter mu')
    # Water outside pipe
    probdata.add_param('rho3',       1000.0,  'density') #kg/m^3 ==> Speeds in m/s
    probdata.add_param('lambda3',    2.202256e-3,  'Lame parameter lambda')
    probdata.add_param('mu3',        0.0,  'Lame parameter mu')
    # Uncured concrete outside pipe 
#    probdata.add_param('rho3',         rho3,    'scaled density') # Scaled Concrete
#    probdata.add_param('lambda3',      lambda3,   'scaled Lame parameter lambda')
#    probdata.add_param('mu3',          mu3,    'scaled Lame parameter mu')
    probdata.add_param('t0wall',      1.0,  'time duration of force') # For transducer pulse at boundary
    probdata.add_param('amplitude',   1.0e-5,  'max amplitude of force') # For transducer pulse at boundary
#    probdata.add_param('pulse_span',  2.0e-5, 'time span of initial pulse')
    probdata.add_param('pulse_span',  0.0, 'no source')
    probdata.add_param('trans_halfdepth', trans_halfdepth, 'half-depth of the transducer')
    probdata.add_param('trans_halfwidth', trans_halfwidth, 'half-width of the transducer')
    probdata.add_param('pipe_inner_radius', pipe_inner_radius, 'inner radius of pipe')
    probdata.add_param('pipe_outer_radius', pipe_outer_radius, 'outer radius of pipe')
        
    #------------------------------------------------------------------
    # Standard Clawpack parameters to be written to claw.data:
    #------------------------------------------------------------------

    clawdata = rundata.clawdata  # initialized when rundata instantiated


#    # ---------------
#    # Spatial domain (focus on resolving transducer):
#    # ---------------
#
#    # Number of space dimensions:
#    clawdata.num_dim = num_dim
#    
#    # Lower and upper edge of computational domain:
#    clawdata.lower[0] = trans_halfdepth                            # xlower
#    clawdata.upper[0] = domain_radius                          # xupper
#    clawdata.lower[1] = -0.5*(1.0 - trans_halfdepth/domain_radius) # ylower
#    clawdata.upper[1] = 0.5*(1.0 - trans_halfdepth/domain_radius)  # yupper
#    clawdata.lower[2] = -0.5*(clawdata.upper[0]-clawdata.lower[0]) # zlower
#    clawdata.upper[2] = 0.5*(clawdata.upper[0]-clawdata.lower[0])  # zupper
#   
#    # Number of grid cells:
#    num_cells_transducer = 11.0
#    clawdata.num_cells[2] = int(np.floor(num_cells_transducer/trans_halfwidth*(clawdata.upper[2]-clawdata.lower[2])))      # mz
#    clawdata.num_cells[0] = clawdata.num_cells[2]      # mx    
#    clawdata.num_cells[1] = clawdata.num_cells[2]      # mz
    
    # ---------------
    # Spatial domain (focus on resolving pipe):
    # ---------------

    # Number of space dimensions:
    clawdata.num_dim = num_dim
    
    if (use_symmetry == 1):
        num_cells_pipe = 1.0
        dx = (pipe_outer_radius-pipe_inner_radius)/num_cells_pipe
#        num_cells_interior = np.ceil((pipe_inner_radius-trans_halfdepth)/dx)
        num_cells_interior = np.ceil((pipe_inner_radius-0.16)/dx)
        num_cells_exterior = np.ceil((domain_radius-pipe_outer_radius)/dx)
        # Lower and upper edge of computational domain:
        clawdata.lower[0] = pipe_inner_radius - num_cells_interior*dx # xlower
        clawdata.upper[0] = pipe_outer_radius + num_cells_exterior*dx   # xupper
        clawdata.lower[1] = -(clawdata.upper[0]-clawdata.lower[0])  # ylower
        clawdata.upper[1] = 0.0                                            # yupper
        clawdata.lower[2] = 0.0                                            # zlower
        clawdata.upper[2] = (clawdata.upper[0]-clawdata.lower[0])       # zupper
        
        # Number of grid cells:
        clawdata.num_cells[0] = int(num_cells_interior + num_cells_pipe + num_cells_exterior)      # mx
        clawdata.num_cells[1] = int(clawdata.num_cells[0])      # my    
        clawdata.num_cells[2] = int(clawdata.num_cells[0])      # mz

    else:
        num_cells_pipe = 6.0
        dx = (pipe_outer_radius-pipe_inner_radius)/num_cells_pipe
        num_cells_interior = np.ceil((pipe_inner_radius-trans_halfdepth)/dx)
        num_cells_exterior = np.ceil((domain_radius-pipe_outer_radius)/dx)
        # Lower and upper edge of computational domain:
        clawdata.lower[0] = pipe_inner_radius - num_cells_interior*dx # xlower
        clawdata.upper[0] = pipe_outer_radius + num_cells_exterior*dx   # xupper
        clawdata.lower[1] = -0.5*(clawdata.upper[0] - clawdata.lower[0]) # ylower
        clawdata.upper[1] = 0.5*(clawdata.upper[0] - clawdata.lower[0])  # yupper
        clawdata.lower[2] = -0.5*(clawdata.upper[0] - clawdata.lower[0])      # zlower
        clawdata.upper[2] = 0.5*(clawdata.upper[0] - clawdata.lower[0])      # zupper
        
        # Number of grid cells:
        clawdata.num_cells[0] = int(num_cells_interior + num_cells_pipe + num_cells_exterior)      # mx
        clawdata.num_cells[1] = clawdata.num_cells[0]      # my    
        clawdata.num_cells[2] = clawdata.num_cells[0]      # mz

    # ---------------
    # Size of system:
    # ---------------

    # Number of equations in the system:
    clawdata.num_eqn = 9

    # Number of auxiliary variables in the aux array (initialized in setaux)
    clawdata.num_aux = 12
    
    # Index of aux array corresponding to capacity function, if there is one:
    clawdata.capa_index = 12
    
    
    # -------------
    # Initial time:
    # -------------

    # should be set to final time in reload file if reloading
    #clawdata.t0 = 0.000000
    #clawdata.t0 = 20.0
    clawdata.t0 = 48.0
    

    # Restart from checkpoint file of a previous run?
    # Note: If restarting, you must also change the Makefile to set:
    #    RESTART = True
    # If restarting, t0 above should be from original run, and the
    # restart_file 'fort.chkNNNNN' specified below should be in 
    # the OUTDIR indicated in Makefile.

    clawdata.restart = False               # True to restart from prior results
    clawdata.restart_file = 'fort.chk00006'  # File to use for restart data
    
    
    # -------------
    # Output times:
    #--------------

    # Specify at what times the results should be written to fort.q files.
    # Note that the time integration stops after the final output time.
 
    clawdata.output_style = 1
 
    if clawdata.output_style==1:
        # Output ntimes frames at equally spaced times up to tfinal:
        # Can specify num_output_times = 0 for no output
        clawdata.num_output_times = 20 #80
        clawdata.tfinal = 52.0 #6e-5 #0.0001800000 #0.00003
        clawdata.output_t0 = True  # output at initial (or restart) time?
        
    elif clawdata.output_style == 2:
        # Specify a list or numpy array of output times:
        # Include t0 if you want output at the initial time.
        clawdata.output_times =  [48.0, 65.0, 90.0, 135.0, 180.0]
 
    elif clawdata.output_style == 3:
        # Output every step_interval timesteps over total_steps timesteps:
        clawdata.output_step_interval = 1
        clawdata.total_steps = 10
        clawdata.output_t0 = True  # output at initial (or restart) time?
        

    clawdata.output_format = 'ascii'      # 'ascii' or 'binary' 

    clawdata.output_q_components = 'all'   # could be list such as [True,True]
    clawdata.output_aux_components = 'none'  # could be list
    clawdata.output_aux_onlyonce = True    # output aux arrays only at t0
    

    # ---------------------------------------------------
    # Verbosity of messages to screen during integration:  
    # ---------------------------------------------------

    # The current t, dt, and cfl will be printed every time step
    # at AMR levels <= verbosity.  Set verbosity = 0 for no printing.
    clawdata.verbosity = 1
    
    

    # --------------
    # Time stepping:
    # --------------

    # if dt_variable==True:  variable time steps used based on cfl_desired,
    # if dt_variable==False: fixed time steps dt = dt_initial always used.
    clawdata.dt_variable = True
    
    # Initial time step for variable dt.  
    # (If dt_variable==0 then dt=dt_initial for all steps)
    clawdata.dt_initial = 1.00000e-05
    
    # Max time step to be allowed if variable dt used:
    clawdata.dt_max = 1.000000e+99
    
    # Desired Courant number if variable dt used 
    clawdata.cfl_desired = 0.500000
    # max Courant number to allow without retaking step with a smaller dt:
    clawdata.cfl_max = 1.000000
    
    # Maximum number of time steps to allow between output times:
    clawdata.steps_max = 500000


    # ------------------
    # Method to be used:
    # ------------------

    # Order of accuracy:  1 => Godunov,  2 => Lax-Wendroff plus limiters
    clawdata.order = 1
    
    # Use dimensional splitting? 
    clawdata.dimensional_split = False
    
    # For unsplit method, transverse_waves can be 
    #  0 or 'none'      ==> donor cell (only normal solver used)
    #  1 or 'increment' ==> corner transport of waves
    #  2 or 'all'       ==> corner transport of 2nd order corrections too
    clawdata.transverse_waves = 0
    
    
    # Number of waves in the Riemann solution:
    clawdata.num_waves = 6
    
    # List of limiters to use for each wave family:  
    # Required:  len(limiter) == num_waves
    # Some options:
    #   0 or 'none'     ==> no limiter (Lax-Wendroff)
    #   1 or 'minmod'   ==> minmod
    #   2 or 'superbee' ==> superbee
    #   3 or 'vanleer'  ==> van Leer
    #   4 or 'mc'       ==> MC limiter
    clawdata.limiter = ['vanleer', 'vanleer', 'vanleer', 'vanleer', 'vanleer','vanleer']
    
    clawdata.use_fwaves = False    # True ==> use f-wave version of algorithms
    
    # Source terms splitting:
    #   src_split == 0 or 'none'    ==> no source term (src routine never called)
    #   src_split == 1 or 'godunov' ==> Godunov (1st order) splitting used, 
    #   src_split == 2 or 'strang'  ==> Strang (2nd order) splitting used,  not recommended.
    clawdata.source_split = 0  # use axi-symmetric source initially
    
    
    # --------------------
    # Boundary conditions:
    # --------------------

    # Number of ghost cells (usually 2)
    clawdata.num_ghost = 2
    
    # Choice of BCs at xlower and xupper:
    #   0 or 'user'     => user specified (must modify bcNamr.f to use this option)
    #   1 or 'extrap'   => extrapolation (non-reflecting outflow)
    #   2 or 'periodic' => periodic (must specify this at both boundaries)
    #   3 or 'wall'     => solid wall for systems where q(2) is normal velocity
    
    if (use_symmetry == 1):
        clawdata.bc_lower[0] = 'user'   # at xlower
        clawdata.bc_upper[0] = 'extrap'   # at xupper
    
        clawdata.bc_lower[1] = 'extrap'   # at ylower   
        clawdata.bc_upper[1] = 'wall'   # at yupper -- for half of symmetric domain
        
        clawdata.bc_lower[2] = 'wall'   # at zlower   -- for half of symmetric domain
        clawdata.bc_upper[2] = 'extrap'   # at zupper
    
    else:
        
        clawdata.bc_lower[0] = 'user'   # at xlower
        clawdata.bc_upper[0] = 'extrap'   # at xupper
    
        clawdata.bc_lower[1] = 'extrap'   # at ylower
        clawdata.bc_upper[1] = 'extrap'   # at yupper
        
        clawdata.bc_lower[2] = 'extrap'   # at zlower
        clawdata.bc_upper[2] = 'extrap'   # at zupper
        
                  
                  
    # --------------
    # Checkpointing:
    # --------------

    # Specify when checkpoint files should be created that can be
    # used to restart a computation.

    clawdata.checkpt_style = 2

    if clawdata.checkpt_style == 0:
      # Do not checkpoint at all
      pass

    elif clawdata.checkpt_style == 1:
      # Checkpoint only at tfinal.
      pass

    elif clawdata.checkpt_style == 2:
      # Specify a list of checkpoint times.  
      clawdata.checkpt_times = [53.]

    elif clawdata.checkpt_style == 3:
      # Checkpoint every checkpt_interval timesteps (on Level 1)
      # and at the final time.
      clawdata.checkpt_interval = 5

       
    # ---------------
    # Gauges:
    # ---------------
    gauges = rundata.gaugedata.gauges
    
    # Gauges across transducer face
    if (use_symmetry == 1):
        ngauges1d = 6
        gcount = 0
        h = trans_halfwidth/(ngauges1d - 1.0)
        x = trans_halfdepth
        for k,x in enumerate(np.hstack((np.linspace(0.16,0.18,5),np.linspace(0.19,0.215,4)))):
            for i in range(ngauges1d):
                y = -trans_halfwidth + i*h
                for j in range(ngauges1d):
                    z = j*h
                    if (np.sqrt(y*y+z*z) <= trans_halfwidth):
                        gauges.append([gcount, x, y, z, 0.0, 1e9])
                        gcount = gcount + 1

    else:
        ngauges1d = 12
        gcount = 0
        h = 2.0*trans_halfwidth/(ngauges1d - 1.0)
        x = trans_halfdepth
        for k,x in enumerate(np.hstack((np.linspace(trans_halfdepth,0.18,5),np.linspace(0.19,0.22,4)))):
            for i in range(ngauges1d):
                y = -trans_halfwidth + i*h
                for j in range(ngauges1d):
                    z = -trans_halfwidth + j*h
                    if (np.sqrt(y*y+z*z) <= trans_halfwidth):
                        gauges.append([gcount, x, y, z, 0.0, 1e9])
                        gcount = gcount + 1
                    
    # -----------------
    # Slice parameters:
    # -----------------

    rundata.slicedata.slices_xy = [0.0]
    rundata.slicedata.slices_xz = [0.0]    
    rundata.slicedata.slices_yz = [0.19]    
        
    # ---------------
    # AMR parameters:
    # ---------------
    amrdata = rundata.amrdata

    # max number of refinement levels:
    amrdata.amr_levels_max = 2
    
    # List of refinement ratios at each level (length at least amr_level_max-1)
    amrdata.refinement_ratios_x = [8]
    amrdata.refinement_ratios_y = [8]
    amrdata.refinement_ratios_z = [8]
    amrdata.refinement_ratios_t = [8]

    # Specify type of each aux variable in amrdata.auxtype.
    # This must be a list of length num_aux, each element of which is one
    # of:
    #   'center',  'capacity', 'xleft', or 'yleft'  (see documentation).
    amrdata.aux_type = ['center', 'center', 'center', 'center', 'xleft', 'xleft', 'xleft', 'yleft', 'yleft', 'yleft', 'zleft','capacity']

    # Flag for refinement based on Richardson error estimater:
    amrdata.flag_richardson = False    # use Richardson?
    amrdata.flag_richardson_tol = 1.000000e+00  # Richardson tolerance
    
    # Flag for refinement using routine flag2refine:
    amrdata.flag2refine = True      # use this?
    amrdata.flag2refine_tol = 5.0e-8  # tolerance used in this routine
    # Flags regions to be refined where tol is roughly the smallest amplitude of normal velocity to capture

    # steps to take on each level L between regriddings of level L+1:
    amrdata.regrid_interval = 2       

    # width of buffer zone around flagged points:
    # (typically the same as regrid_interval so waves don't escape):
    amrdata.regrid_buffer_width  = 2

    # clustering alg. cutoff for (# flagged pts) / (total # of cells
    # refined)
    # (closer to 1.0 => more small grids may be needed to cover flagged
    # cells)
    amrdata.clustering_cutoff = 0.700000

    # print info about each regridding up to this level:
    amrdata.verbosity_regrid = 4


    # -------------------
    # Refinement Regions:
    # -------------------
    #regions = rundata.regiondata.regions 
    # to specify regions of refinement append lines of the form
    #  [minlevel,maxlevel,t1,t2,x1,x2,y1,y2]
    #regions.append([1,3,0,1e9,0.,1.,-1,1])  # whole domain
    #regions.append([2,5,0,1e9,0.,1.,-0.5, 0.5])  # near centerline
    #regions.append([2,5,0,1e9,0.,.02,-yt,yt])  # near transducer


    #  ----- For developers ----- 
    # Toggle debugging print statements:
    amrdata.dprint = False      # print domain flags
    amrdata.eprint = False      # print err est flags
    amrdata.edebug = False      # even more err est flags
    amrdata.gprint = False      # grid bisection/clustering
    amrdata.nprint = False      # proper nesting output
    amrdata.pprint = False      # proj. of tagged points
    amrdata.rprint = False      # print regridding summary
    amrdata.sprint = False      # space/memory output
    amrdata.tprint = False      # time step reporting each level
    amrdata.uprint = False      # update/upbnd reporting
    
    
    return rundata

    # end of function setrun
    # ----------------------


if __name__ == '__main__':
    # Set up run-time parameters and write all data files.
    import sys
    rundata = setrun(*sys.argv[1:])
    rundata.write()
    
