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


# Figure 3

start_time = time.time()
fname = "step-{timestep:04}_wavefield.pbin"
path = "/import/freenas-m-03-geodynamics/jhayek/TEAR/Results/T2/Runs/TEAR35_Kos_T20_P3_025x025_A12phi65_Delta2_4s/"


i=4630
FieldFilename = os.path.join(path,fname.format(timestep=i))

MeshFilename = os.path.join(path, "default_mesh_coor.pbin")
se2_coor = se2wave_load_coordinates(MeshFilename)

FileList = glob(os.path.join(path,"step-{timestep}_wavefield.pbin".format(timestep="*")))
l = [i.replace(os.path.join(path,'step-'),'').replace('_wavefield.pbin','') for i in FileList]


TimeStepVal, LCoorX, LCoorY, LFieldX, LFieldY, LFieldvelX, LFieldvelY =  ExtractFields(FieldFilename, se2_coor)


FolderProfilesPath = "/import/freenas-m-03-geodynamics/jhayek/SharedWolfel/PaperData/CorrectedSimulations/20220315/"
DataProfile = LoadPickleFile(Filename = "TEAR35_Kos_T20_P3_025x025_A12phi65_Delta2_4s-Tilt20.0-P3-TPList_t4630_d50.0.pickle",FolderPath = FolderProfilesPath)

x0,y0 = 7350,2675
InsetAxis = [x0-200,x0+200,y0-200,y0+200]
F1, ax = Plot4KomaSetup(LCoorX, LCoorY, LFieldX, LFieldvelX, 
           ["X-Component Displacement ", "X-Component Displacement [m]"],
           TimeStepVal,InsetAxis,
            cmap=cmap)
del x0,y0,InsetAxis


# Tilted case plotting
iidx = 0
for iidx,Test1 in enumerate(DataProfile):
    ax[0].plot(Test1.Time, Test1.DispX, color= cmapProf.colors[iidx], linewidth=1.5, zorder=iidx)
    ax[1].plot(Test1.Time, Test1.VelX,  color= cmapProf.colors[iidx], linewidth=1.5, zorder=iidx) 


ax[0].set_xlabel("time [s]")

#F1.suptitle("Tilting (20deg) Kostrov simulation")
[item.PlotReference(ax[0], "Slip", filtering=False) for item in RefList]
[item.PlotReference(ax[1], "SlipRate", filtering=False) for item in RefList]


Format_LabelsOnFig_formatAxis(F1, ax[:2],inverted=True, ncols = 3, HeightBbox=1.2)


LabelizeAxisList(ax,Pos=[0.9, 0.9],fontsize=BIGGER_SIZE)


# OutFile = "/home/nico/Documents/Documents/SharedWolfel/Works/se2dr_Paper/Illustrations/FinalFigures/F{}.pdf"
# F1.savefig(OutFile.format("3"))
OutFile = "/import/freenas-m-03-geodynamics/jhayek/SharedWolfel/Works/se2dr_Paper/Illustrations/FinalFigures/F{}.png"
F1.savefig(OutFile.format("3"))