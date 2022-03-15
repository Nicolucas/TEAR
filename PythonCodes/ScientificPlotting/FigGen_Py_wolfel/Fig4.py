import numpy as np
import matplotlib.pyplot as plt
plt.style.use('science')

import os, sys, time


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


###################################################################
###################### Reference solution
###################################################################
path = "/home/nico/Documents/TEAR/Codes_TEAR/ProfilePicking/Output/"
# Reference saved into a list of objects
RefList = [SSCreference(path + "Reference/sem2dpack/sem2d-{}-0.txt", "0km"),
           SSCreference(path + "Reference/sem2dpack/sem2d-{}-1.txt", "2km"),
           SSCreference(path + "Reference/sem2dpack/sem2d-{}-2.txt", "4km"),
           SSCreference(path + "Reference/sem2dpack/sem2d-{}-3.txt", "6km"),
           SSCreference(path + "Reference/sem2dpack/sem2d-{}-4.txt", "8km"),
          ]

pathTPV = "/home/nico/Documents/TEAR/Codes_TEAR/ProfilePicking/[TPV3]Results/"

# Reference saved into a list of objects

RefListTPV =  [TPV3reference(pathTPV + "Reference/sem2dpack/[TPV3]sem2dpack-{}-receiver-0.0e+00.txt", "0km"),
               TPV3reference(pathTPV + "Reference/sem2dpack/[TPV3]sem2dpack-{}-receiver-2.0e+03.txt", "2km"),
               TPV3reference(pathTPV + "Reference/sem2dpack/[TPV3]sem2dpack-{}-receiver-4.0e+03.txt", "4km"),
               TPV3reference(pathTPV + "Reference/sem2dpack/[TPV3]sem2dpack-{}-receiver-6.0e+03.txt", "6km"),
               TPV3reference(pathTPV + "Reference/sem2dpack/[TPV3]sem2dpack-{}-receiver-8.0e+03.txt", "8km"),
              ]
###################################################################
###################### Reference solution
###################################################################

# Figure 4

start_time = time.time()
fname = "step-{timestep:04}_wavefield.pbin"
path = "/home/nico/Documents/Documents/SharedWolfel/PaperData/220120FieldData/TEAR18_Kos_Sig_P3_025x025_A12phi65_Delta2.501/"


i=8180
FieldFilename = os.path.join(path,fname.format(timestep=i))

MeshFilename = os.path.join(path, "default_mesh_coor.pbin")
se2_coor = se2wave_load_coordinates(MeshFilename)



FileList = glob(os.path.join(path,"step-{timestep}_wavefield.pbin".format(timestep="*")))
l = [i.replace(os.path.join(path,'step-'),'').replace('_wavefield.pbin','') for i in FileList]

TimeStepVal, LCoorX, LCoorY, LFieldX, LFieldY, LFieldvelX, LFieldvelY =  ExtractFields(FieldFilename, se2_coor)

FolderTiltedPath = "/home/nico/Documents/Documents/SharedWolfel/PaperData/ConvTPV/Sigmoid-20220124/"
DataProfile = LoadPickleFile(Filename = "TEAR18_Kos_Sig_P3_025x025_A12phi65_Delta2.501-25x25-P3-62.525",FolderPath = FolderTiltedPath)
StressFromPickle = LoadPickleFile("/home/nico/Documents/Documents/SharedWolfel/PaperData/220120FieldData/TEAR18_Kos_Sig_P3_025x025_A12phi65_Delta2.501/Out/", "StressInAPickle")

DataProfile.reverse()

x0,y0 = 6600,3610
InsetAxis = [x0-200,x0+200,y0-200,y0+200]

F1, ax = PlotF4Setup(LCoorX, LCoorY, LFieldvelX, StressFromPickle,  
           ["X-Component Displacement ", "X-Component Displacement [m]"],
           TimeStepVal,InsetAxis,cmap=cmap)
del x0,y0,InsetAxis
# Tilted case plotting
iidx = 0
for iidx,Test1 in enumerate(DataProfile):
    ax[0].plot(Test1.Time, Test1.Slip, color= cmapProf.colors[iidx], linewidth=1.5, zorder=iidx)
    ax[1].plot(Test1.Time, Test1.SlipRate,  color= cmapProf.colors[iidx], linewidth=1.5, zorder=iidx) 


ax[0].set_xlabel("time [s]")

# F1.suptitle("Sigmoid Kostrov simulation")
[item.PlotReference(ax[0], "Slip", filtering=False) for item in RefList]
[item.PlotReference(ax[1], "SlipRate", filtering=False) for item in RefList]



Format_LabelsOnFig_formatAxis(F1, ax[:2],inverted=True, ncols = 3)

LabelizeAxisList(ax,Pos=[0.9, 0.9],fontsize=BIGGER_SIZE)

OutFile = "/home/nico/Documents/Documents/SharedWolfel/Works/se2dr_Paper/Illustrations/FinalFigures/F{}.pdf"
F1.savefig(OutFile.format("4"))