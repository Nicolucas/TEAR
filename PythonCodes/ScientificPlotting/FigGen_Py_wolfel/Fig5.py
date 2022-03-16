#!/usr/bin/env python
# coding: utf-8

# # se2dr Fig 2 and 5 plot compound

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


## Extract the information for the profiles

###################################################################
###################### Reference solution
###################################################################
pathRef = "/import/freenas-m-03-geodynamics/jhayek/SharedWolfel/PaperData/References/feik/"
# Reference saved into a list of objects
RefList = [SSCreference(pathRef + "Kos_sem2dpack-{}-receiver-0.txt", "0km"),
           SSCreference(pathRef + "Kos_sem2dpack-{}-receiver-1.txt", "2km"),
           SSCreference(pathRef + "Kos_sem2dpack-{}-receiver-2.txt", "4km"),
           SSCreference(pathRef + "Kos_sem2dpack-{}-receiver-3.txt", "6km"),
           SSCreference(pathRef + "Kos_sem2dpack-{}-receiver-4.txt", "8km"),
          ]
# Reference saved into a list of objects
RefListTPV =  [TPV3reference(pathRef + "[TPV3]sem2dpack-{}-receiver-0.0e+00.txt", "0km"),
               TPV3reference(pathRef + "[TPV3]sem2dpack-{}-receiver-2.0e+03.txt", "2km"),
               TPV3reference(pathRef + "[TPV3]sem2dpack-{}-receiver-4.0e+03.txt", "4km"),
               TPV3reference(pathRef + "[TPV3]sem2dpack-{}-receiver-6.0e+03.txt", "6km"),
               TPV3reference(pathRef + "[TPV3]sem2dpack-{}-receiver-8.0e+03.txt", "8km"),
              ]
###################################################################
###################### Reference solution
###################################################################

## Now select the time snapshot of interest

start_time = time.time()
fname = "step-{timestep:04}_wavefield.pbin"
path = "/import/freenas-m-03-geodynamics/jhayek/TEAR/Results/T2/Runs/TEAR35_TPV_T0_P3_025x025_A12phi65_Delta1.001_3s/"


i=4630
FieldFilename = os.path.join(path,fname.format(timestep=i))

MeshFilename = os.path.join(path, "default_mesh_coor.pbin")
se2_coor = se2wave_load_coordinates(MeshFilename)

FileList = glob(os.path.join(path,"step-{timestep}_wavefield.pbin".format(timestep="*")))
l = [i.replace(os.path.join(path,'step-'),'').replace('_wavefield.pbin','') for i in FileList]


## Extract the fields for velocity and displacement in each component
TimeStepVal, LCoorX, LCoorY, LFieldX, LFieldY, LFieldvelX, LFieldvelY =  ExtractFields(FieldFilename, se2_coor)


FolderProfilesPath = "/import/freenas-m-03-geodynamics/jhayek/SharedWolfel/PaperData/CorrectedSimulations/20220120/"
Profile2Plot = LoadPickleFile(Filename = "TEAR14_Kos_T0_P3_025x025_A12phi65_5Sec-Tilt0.0-P3-TPList_t9990_d25.025.pickle",FolderPath = FolderProfilesPath)
StressFromPickle = LoadPickleFile(path+"/Out/", "StressInAPickle")

F1, ax = PlotFullSetup(LCoorX, LCoorY, LFieldX, LFieldvelX, StressFromPickle, 
           ["X-Component Displacement ", "X-Component Displacement [m]"],
           TimeStepVal,[8000-200,8000+200,-200,200],
            cmap=cmap)

# Tilted case plotting
iidx = 0
for iidx,Test1 in enumerate(Profile2Plot):
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

# OutFile = "/home/nico/Documents/Documents/SharedWolfel/Works/se2dr_Paper/Illustrations/FinalFigures/F{}.pdf"
# F1.savefig(OutFile.format("5"))
OutFile = "/import/freenas-m-03-geodynamics/jhayek/SharedWolfel/Works/se2dr_Paper/Illustrations/FinalFigures/F{}.png"
F1.savefig(OutFile.format("5"))