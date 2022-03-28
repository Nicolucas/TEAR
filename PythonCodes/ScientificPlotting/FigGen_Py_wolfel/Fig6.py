import numpy as np
import matplotlib.pyplot as plt
plt.style.use('science')

import os, sys, time
from mpl_toolkits.axes_grid1.inset_locator import (inset_axes, TransformedBbox,
                                                   BboxPatch, BboxConnector)


sys.path.insert(0,"/import/freenas-m-03-geodynamics/jhayek/petsc-3.12.5/lib/petsc/bin/")
sys.path.insert(0,"/import/freenas-m-03-geodynamics/jhayek/TEAR/se2wave/utils/python")
sys.path.insert(0,"/import/freenas-m-03-geodynamics/jhayek/TEAR/processing/TEAR/PythonCodes/")
from se2waveload import *
from Lib_GeneralFunctions import *
from GeneratePaperFigs import *
from ModelIllustration import *


SMALL_SIZE = 14
MEDIUM_SIZE = 16
BIGGER_SIZE = 20
FontSizeControlFreak(SMALL_SIZE,MEDIUM_SIZE,BIGGER_SIZE)


from palettable.colorbrewer.diverging import PuOr_11_r as FieldColor

cmap = FieldColor.mpl_colormap

from matplotlib.colors import ListedColormap
import matplotlib.lines as mlines
from palettable.cartocolors.qualitative import Safe_5 as LineColor

cmapProf = ListedColormap(LineColor.mpl_colors[:])


def my_mark_inset(parent_axes, inset_axes, loc1a=1, loc1b=1, loc2a=2, loc2b=2, **kwargs):
    rect = TransformedBbox(inset_axes.viewLim, parent_axes.transData)
    pp = BboxPatch(rect, fill=False, **kwargs)
    parent_axes.add_patch(pp)
    p1 = BboxConnector(inset_axes.bbox, rect, loc1=loc1a, loc2=loc1b, **kwargs)
    inset_axes.add_patch(p1)
    p1.set_clip_on(False)
    p2 = BboxConnector(inset_axes.bbox, rect, loc1=loc2a, loc2=loc2b, **kwargs)
    inset_axes.add_patch(p2)
    p2.set_clip_on(False)

    return pp, p1, p2

def scientific(x, pos):
            # x:  tick value - ie. what you currently see in yticks
            # pos: a position - ie. the index of the tick (from 0 to 9 in this example)
            return '%.2E' % x

def PlotVelSweep(fname,path,Frames,PathPost,ZoomIn, VMaxMin=None,**kwargs):  

    fig = plt.figure(figsize = (9, 19),dpi=350) #constrained_layout=True
    gs = fig.add_gridspec(3, 1, wspace=0,hspace=0.12)

    ConnectAxis = []
    Ax4xLabel = []
    for i,gsSet in reversed(list(enumerate(gs))):
        gsSet = gsSet.subgridspec(2, 3, wspace=0,hspace=0.)

        ax0_00 = fig.add_subplot(gsSet[0,2])
        ax0_01 = fig.add_subplot(gsSet[1,2])
        ax0_10 = fig.add_subplot(gsSet[0:2, :2])
        ax0_10.yaxis.tick_left()
        ax0_00.yaxis.tick_right()
        ax0_01.yaxis.tick_right()

        ax0_00.set_xticklabels([])
        
        FormatAxNormal(ax0_00)
        FormatAxNormal(ax0_01)
        FormatAxNormal(ax0_10)
        

        
        ######## Load the necessary files for Ax0_10
        MeshFilename = os.path.join(path.format(PathPost[1]), "default_mesh_coor.pbin")
        se2_coor = se2wave_load_coordinates(MeshFilename)
        FieldFilename = os.path.join(path.format(PathPost[1]),fname.format(timestep=Frames[i]))
        TimeStepVal, CoorX, CoorY, FieldX, FieldY, FieldvelX, FieldvelY =  ExtractFields(FieldFilename, se2_coor)
        MeshFilename = os.path.join(path, "default_mesh_coor.pbin")
        
        # Get the values to normalize the rest of the subset of plots
        ############# Get max and min values of the velocity field
        if VMaxMin is None:
            Vmax  = FieldvelX.max()
            Vmin  = FieldvelX.min()
        else:
            Vmax  =  VMaxMin[i]
            Vmin  = -VMaxMin[i]
        ############# 
        
        
        img3 = ax0_10.pcolormesh(CoorX, CoorY, FieldvelX, vmax = Vmax, vmin = Vmin,**kwargs)
        
        # Colorbar
        cbaxes = inset_axes(ax0_10, width="40%",height="4%",loc=3, borderpad=2)
        plt.colorbar(img3,cax=cbaxes,orientation="horizontal", label=r"$v_{x}$ [m/s]")
        cbaxes.xaxis.set_label_position('top')
        
        # Timestamp
        ax0_10.annotate(text="time [s]: {0:.2f}".format(TimeStepVal),xy=[0.05,0.9], xycoords= "axes fraction")
        
        # Stamp for fault zone thickness
        ax0_10.annotate(text="$\delta={}\ h$".format(PathPost[1]),xy=[0.75,0.05], xycoords= "axes fraction")
        
        
        
        ######## Load the necessary files for Ax0_00
        MeshFilename = os.path.join(path.format(PathPost[0]), "default_mesh_coor.pbin")
        se2_coor = se2wave_load_coordinates(MeshFilename)
        FieldFilename = os.path.join(path.format(PathPost[0]),fname.format(timestep=Frames[i]))
        TimeStepVal, CoorX, CoorY, FieldX, FieldY, FieldvelX, FieldvelY =  ExtractFields(FieldFilename, se2_coor)
        MeshFilename = os.path.join(path, "default_mesh_coor.pbin")

        img1 = ax0_00.pcolormesh(CoorX, CoorY, FieldvelX, vmax = Vmax, vmin = Vmin, **kwargs)
        
        # Stamp for fault zone thickness
        ax0_00.annotate(text="$\delta={}\ h$".format(PathPost[0]),xy=[0.55,0.1], xycoords= "axes fraction")

        ######## Load the necessary files for Ax0_01
        MeshFilename = os.path.join(path.format(PathPost[2]), "default_mesh_coor.pbin")
        se2_coor = se2wave_load_coordinates(MeshFilename)
        FieldFilename = os.path.join(path.format(PathPost[2]),fname.format(timestep=Frames[i]))
        TimeStepVal, CoorX, CoorY, FieldX, FieldY, FieldvelX, FieldvelY =  ExtractFields(FieldFilename, se2_coor)
        MeshFilename = os.path.join(path, "default_mesh_coor.pbin")

        img2 = ax0_01.pcolormesh(CoorX, CoorY, FieldvelX, vmax = Vmax, vmin = Vmin,**kwargs)
        
        # Stamp for fault zone thickness
        ax0_01.annotate(text="$\delta={}\ h$".format(PathPost[2]),xy=[0.55,0.1], xycoords= "axes fraction")

        
        
        
        ax0_10.set_xlim([-ZoomIn[i],ZoomIn[i]])
        ax0_00.set_xlim([-ZoomIn[i],ZoomIn[i]])
        ax0_01.set_xlim([-ZoomIn[i],ZoomIn[i]])
        ax0_10.set_ylim([-ZoomIn[i],ZoomIn[i]])
        ax0_00.set_ylim([-ZoomIn[i],ZoomIn[i]])
        ax0_01.set_ylim([-ZoomIn[i],ZoomIn[i]])
        
        ConnectAxis.append(ax0_10)
        Ax4xLabel.append(ax0_01)
        
        ax0_10.set_ylabel("$y$ [m]")
        
    
    my_mark_inset(ConnectAxis[0],ConnectAxis[1],loc1a=3, loc1b=2, loc2a=4, loc2b=1, alpha=0.3,ls="--")
    my_mark_inset(ConnectAxis[1],ConnectAxis[2],loc1a=3, loc1b=2, loc2a=4, loc2b=1, alpha=0.3,ls="--")
    
    ConnectAxis[0].set_xlabel("$x$ [m]")
    Ax4xLabel[0].set_xlabel("$x$ [m]")

    return fig
        
def PlotVelSweepX2(fname,path,Frames,PathPost,ZoomIn,VMaxMin=None,**kwargs):  

    fig = plt.figure(figsize = (9, 15),dpi=300) #constrained_layout=True
    gs = fig.add_gridspec(3, 1, wspace=0,hspace=0.12)
    
    for i,gsSet in enumerate(gs):
        gsSet = gsSet.subgridspec(1, 2, wspace=0,hspace=0.)

        ax0_00 = fig.add_subplot(gsSet[0,0])
        ax0_10 = fig.add_subplot(gsSet[0,1])
        ax0_10.yaxis.tick_right()

        ax0_00.set_xticklabels([])
        ax0_00.set_yticklabels([])
        
        FormatAxNormal(ax0_00)
        FormatAxNormal(ax0_10)
        
        
        ######## Load the necessary files for Ax0_10
        MeshFilename = os.path.join(path.format(PathPost[1]), "default_mesh_coor.pbin")
        se2_coor = se2wave_load_coordinates(MeshFilename)
        FieldFilename = os.path.join(path.format(PathPost[1]),fname.format(timestep=Frames[i]))
        TimeStepVal, CoorX, CoorY, FieldX, FieldY, FieldvelX, FieldvelY =  ExtractFields(FieldFilename, se2_coor)
        MeshFilename = os.path.join(path, "default_mesh_coor.pbin")
        
        ############# Get max and min values of the velocity field
        if VMaxMin is None:
            Vmax  = FieldvelX.max()
            Vmin  = FieldvelX.min()
        else:
            Vmax  =  VMaxMin[i]
            Vmin  = -VMaxMin[i]
        ############# 
        
        img3 = ax0_10.pcolormesh(CoorX, CoorY, FieldvelX, vmax = Vmax, vmin = Vmin, **kwargs)
        
        
        # Stamp for fault zone thickness
        ax0_10.annotate(text="$\delta={}\ h$".format(PathPost[1]),xy=[0.7,0.05], xycoords= "axes fraction")
        
        
        
        ######## Load the necessary files for Ax0_00
        MeshFilename = os.path.join(path.format(PathPost[0]), "default_mesh_coor.pbin")
        se2_coor = se2wave_load_coordinates(MeshFilename)
        FieldFilename = os.path.join(path.format(PathPost[0]),fname.format(timestep=Frames[i]))
        TimeStepVal, CoorX, CoorY, FieldX, FieldY, FieldvelX, FieldvelY =  ExtractFields(FieldFilename, se2_coor)
        MeshFilename = os.path.join(path, "default_mesh_coor.pbin")

        img1 = ax0_00.pcolormesh(CoorX, CoorY, FieldvelX, vmax = Vmax, vmin = Vmin, **kwargs)
        
        # Stamp for fault zone thickness
        ax0_00.annotate(text="$\delta={}\ h$".format(PathPost[0]),xy=[0.7,0.05], xycoords= "axes fraction")

        # Colorbar
        cbaxes = inset_axes(ax0_00, width="40%",height="4%",loc=3, borderpad=2)
        plt.colorbar(img3,cax=cbaxes,orientation="horizontal", label=r"$v_{x}$ [m/s]")
        cbaxes.xaxis.set_label_position('top')
        
        # Timestamp
        ax0_00.annotate(text="time [s]: {0:.2f}".format(TimeStepVal),xy=[0.05,0.9], xycoords= "axes fraction")
        
        
        
        ax0_10.set_xlim([-ZoomIn[i],ZoomIn[i]])
        ax0_00.set_xlim([-ZoomIn[i],ZoomIn[i]])
        ax0_10.set_ylim([-ZoomIn[i],ZoomIn[i]])
        ax0_00.set_ylim([-ZoomIn[i],ZoomIn[i]])

    return fig


start_time = time.time()
fname = "step-{timestep:04}_wavefield.pbin"
path = "/import/freenas-m-03-geodynamics/jhayek/TEAR/Results/T2/Runs/TEAR35_TPV_T20_P3_025x025_A12phi65_Delta{}_3s/"

PathPost = ['1.0','1.43','4.00'];Frames = [ 50, 600, 2300];ZoomIn = [2000,3400,7000];

F1 = PlotVelSweep(fname,path,Frames,PathPost,ZoomIn,VMaxMin=[0.1,1.5,2],cmap=cmap)

OutFile = "/import/freenas-m-03-geodynamics/jhayek/SharedWolfel/Works/se2dr_Paper/Illustrations/FinalFigures/F{}.png"
F1.savefig(OutFile.format("6"))