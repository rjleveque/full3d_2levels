
""" 
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.
    
""" 
from mapc2p import mapc2p
import numpy, os, shutil
from pylab import linspace

from clawpack.clawutil.data import ClawData

# Use 1 for xy, 2 for xz
slice_type = 2


#--------------------------
def setplot(plotdata):
#--------------------------
    
    """ 
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of clawpack.visclaw.data.ClawPlotData.
    Output: a modified version of plotdata.
    
    """ 
    from clawpack.visclaw import colormaps

    plotdata.outdir = './_output'
    plotdata.clearfigures()  # clear any old figures,axes,items data

#   copy slice files to fort files
    if (slice_type == 1):
        prefix = "slice_xy1"
    elif (slice_type == 2):
        prefix = "slice_xz1"

    os.chdir(plotdata.outdir)
    for filename in os.listdir("."):
        if filename.startswith(prefix):
            shutil.copyfile(filename,filename.replace(prefix,"fort",1))
    os.chdir("..")

    probdata = ClawData()
    probdata.read('setprob.data',force=True)
    griddata = ClawData()
    griddata.read('claw.data',force=True)
        
    #Plot outline of interface
    def interface(current_data):
        from pylab import linspace,plot
        from numpy import cos,sin,ones

        if (slice_type == 1):
	     y = linspace(griddata.lower[1],griddata.upper[1],1000)
        elif (slice_type == 2):
	     y = linspace(griddata.lower[2],griddata.upper[2],1000)
             
        x1 = probdata.pipe_inner_radius*ones((1000,1))
        x2 = probdata.pipe_outer_radius*ones((1000,1))
        plot(x1,y,'g',linewidth=2.0)
        plot(x2,y,'g',linewidth=2.0)
          
    # Outputs the trace of the stress tensor
    def sigmatr(current_data):
        q = current_data.q
        return -(q[0,:,:] + q[1,:,:] + q[2,:,:])/3.0
      
    # Figure for stress trace - pcolor
    plotfigure = plotdata.new_plotfigure(name='Minus stress trace', figno=1)
    #plotfigure.show = False

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = 'auto' #[-0.02,0.04] #'auto'
    plotaxes.ylimits = 'auto' #[-0.00,0.1]
    plotaxes.title = 'Stress trace'
    plotaxes.afteraxes = "pylab.axis('scaled')"
    plotaxes.scaled = True      # so aspect ratio is 1
    plotaxes.afteraxes = interface

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.show = True 
    plotitem.plot_var = sigmatr
    plotitem.pcolor_cmap = colormaps.blue_white_red
    plotitem.pcolor_cmin = -9e-6
    plotitem.pcolor_cmax = 9e-6
    plotitem.add_colorbar = True
    #plotitem.amr_patchedges_show = [1]
    #plotitem.amr_celledges_show = [0]
    if (slice_type == 1):
        plotitem.MappedGrid = True
        plotitem.mapc2p = mapc2p
    # Figure for stress trace - pcolor
    plotfigure = plotdata.new_plotfigure(name='Minus stress trace (scale)', figno=2)
    plotfigure.show = False

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = 'auto' #[-0.02,0.04] #'auto'
    plotaxes.ylimits = 'auto' #[-0.00,0.1]
    plotaxes.title = 'Stress trace'
    plotaxes.afteraxes = "pylab.axis('scaled')"
    plotaxes.scaled = True      # so aspect ratio is 1
    plotaxes.afteraxes = interface

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.show = True 
    plotitem.plot_var = sigmatr
    plotitem.pcolor_cmap = colormaps.blue_white_red
    plotitem.pcolor_cmin = -3e-7
    plotitem.pcolor_cmax = 3e-7
    plotitem.add_colorbar = True
    #plotitem.amr_patchedges_show = [1]
    #plotitem.amr_celledges_show = [0]
    if (slice_type == 1):
        plotitem.MappedGrid = True
        plotitem.mapc2p = mapc2p

    ## Figure for grid cells
    plotfigure = plotdata.new_plotfigure(name='cells', figno=3)
    #plotfigure.show = False

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = 'Grid patches'
    plotaxes.scaled = True
    plotaxes.afteraxes = interface
    plotaxes.xlimits = 'auto' #[.06,.09]
    plotaxes.ylimits = 'auto' #[0.0,0.1]

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_patch')
    #plotitem.amr_patch_bgcolor = ['#ffeeee', '#eeeeff', '#eeffee']
    plotitem.amr_celledges_show = [1,0]
    plotitem.amr_patchedges_show = [1]
    if (slice_type == 1):
        plotitem.MappedGrid = True
        plotitem.mapc2p = mapc2p

    # Parameters used only when creating html and/or latex hardcopy
    # e.g., via clawpack.visclaw.frametools.printframes:

    plotdata.printfigs = True                # print figures
    plotdata.print_format = 'png'            # file format
    plotdata.print_framenos = 'all' #range(37)          # list of frames to print
    plotdata.print_fignos = 'all'            # list of figures to print
    plotdata.html = True                     # create html files of plots?
    plotdata.html_homelink = '../README.html'   # pointer for top of index
    plotdata.html_movie = 'JSAnimation'      # new style, or "4.x" for old style
    plotdata.latex = True                    # create latex file of plots?
    plotdata.latex_figsperline = 2           # layout of plots
    plotdata.latex_framesperline = 1         # layout of plots
    plotdata.latex_makepdf = False           # also run pdflatex?

    return plotdata

    
