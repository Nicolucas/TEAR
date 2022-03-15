#!/usr/bin/env python
# coding: utf-8

# # se2dr Fig 2 and 5 plot compound

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
plt.style.use('science')

import os, sys, time


sys.path.insert(0,"/home/nico/Tools/petsc-3.12.5/lib/petsc/bin/")
sys.path.insert(0,"/home/nico/Documents/TEAR/Codes_TEAR/se2dr/se2wave/utils/python/")
sys.path.insert(0,"/home/nico/Documents/TEAR/Codes_TEAR/PythonCodes/LibFolder")
from se2waveload import *
from Lib_GeneralFunctions import *
from GeneratePaperFigs import *
from ModelIllustration import *


SMALL_SIZE = 12
MEDIUM_SIZE = 14
BIGGER_SIZE = 16
FontSizeControlFreak(SMALL_SIZE,MEDIUM_SIZE,BIGGER_SIZE)



# # Colormap selection
from palettable.colorbrewer.diverging import PuOr_11_r as FieldColor

cmap = FieldColor.mpl_colormap

from matplotlib.colors import ListedColormap
import matplotlib.lines as mlines
from palettable.cartocolors.qualitative import Safe_5 as LineColor

cmapProf = ListedColormap(LineColor.mpl_colors[:])


# # Extract the information for the profiles
# 
# ## First the reference

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


# # Now select the time snapshot of interest
start_time = time.time()
fname = "step-{timestep:04}_wavefield.pbin"
path = "/home/nico/Documents/Documents/SharedWolfel/PaperData/220120FieldData/TEAR14_Kos_T0_P3_025x025_A12phi65_5Sec/"


i=8180
FieldFilename = os.path.join(path,fname.format(timestep=i))

MeshFilename = os.path.join(path, "default_mesh_coor.pbin")
se2_coor = se2wave_load_coordinates(MeshFilename)

FileList = glob(os.path.join(path,"step-{timestep}_wavefield.pbin".format(timestep="*")))
l = [i.replace(os.path.join(path,'step-'),'').replace('_wavefield.pbin','') for i in FileList]



# # Extract the fields for velocity and displacement in each component

TimeStepVal, LCoorX, LCoorY, LFieldX, LFieldY, LFieldvelX, LFieldvelY =  ExtractFields(FieldFilename, se2_coor)

FolderTiltedPath = "/home/nico/Documents/Documents/SharedWolfel/PaperData/ConvTPV/20220120/"
ZeroTiltp3h25 = LoadPickleFile(Filename = "TEAR14_Kos_T0_P3_025x025_A12phi65_5Sec-Tilt0.0-P3-TPList_t9990_d25.025.pickle",FolderPath = FolderTiltedPath)
StressFromPickle = LoadPickleFile("/home/nico/Documents/Documents/SharedWolfel/PaperData/220120FieldData/TEAR14_Kos_T0_P3_025x025_A12phi65_5Sec/Out/", "StressInAPickle")

F1, ax = PlotFullSetup(LCoorX, LCoorY, LFieldX, LFieldvelX, StressFromPickle, 
           ["X-Component Displacement ", "X-Component Displacement [m]"],
           TimeStepVal,[8000-200,8000+200,-200,200],
            cmap=cmap)

# Tilted case plotting
iidx = 0
for iidx,Test1 in enumerate(ZeroTiltp3h25):
    ax[0].plot(Test1.Time, Test1.DispX, color= cmapProf.colors[iidx], linewidth=1.5, zorder=iidx)
    ax[1].plot(Test1.Time, Test1.VelX,  color= cmapProf.colors[iidx], linewidth=1.5, zorder=iidx) 

LabelsPerColor= ["25x25 - P3 - $\delta$50."]

#F1.suptitle("Mesh-aligned Kostrov simulation")
[item.PlotReference(ax[0], "Slip", filtering=False) for item in RefList]
[item.PlotReference(ax[1], "SlipRate", filtering=False) for item in RefList]

Format_LabelsOnFig_formatAxis(F1, ax[:2],inverted=True,ncolSim=2,AxLabelLocRef=ax[2])

############### Ax 2
GenKostrovCase(ax[2])

# receivers locations
for i in [0,1,2,3,4]:
    loc=[0.2*i,0.04]

    LocateSlipReceivers_MeshAligned(ax[2],loc,np.array(cmapProf.colors[-i-1]))
###########################################
ax[2].set_xlabel("$x$")
ax[2].set_ylabel("$y$")

LabelizeAxisList(ax,Pos=[0.9, 0.9],fontsize=BIGGER_SIZE)

OutFile = "/home/nico/Documents/Documents/SharedWolfel/Works/se2dr_Paper/Illustrations/FinalFigures/F{}.pdf"
F1.savefig(OutFile.format("2"))
