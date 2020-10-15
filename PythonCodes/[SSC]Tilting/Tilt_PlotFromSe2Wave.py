import os, sys
from glob import glob
from matplotlib import pyplot as plt
from scipy.interpolate import RectBivariateSpline

sys.path.insert(0,"/home/nico/Documents/TEAR/Codes_TEAR/PythonCodes/LibFolder")
from Lib_GeneralFunctions import *
from Lib_ProfilePlotting import *
from Lib_ProfileProcessing import *

from se2waveload import *


NumPoints = 1000
delta = 50.005 #not used
TiltAngle = 45.0
LocIni,LocEnd = [0, delta], [8000, delta]
LocIni,LocEnd = list(ApplyTilting(TiltAngle,LocIni[0],LocIni[1])),list(ApplyTilting(TiltAngle,LocEnd[0],LocEnd[1]))

path = "/home/nico/Documents/TEAR/Codes_TEAR/se2dr/se2wave/"
#path = "/home/nico/Documents/TEAR/Codes_TEAR/plot-utils_se2wave/se2wave/"
#path = "/media/nico/Elements/Simulations/20200728/SSCdeg2/"
#path = "/home/nico/Documents/TEAR/Codes_TEAR/ProfilePicking/Output/20200729/TPV3-P1-Default/"
filename = os.path.join(path,"default_mesh_coor.pbin")

OutputFolder="/home/nico/Documents/TEAR/Codes_TEAR/ProfilePicking/Output/"+GetTodayDate()+"/"
CreateFolder(OutputFolder)

se2_coor = se2wave_load_coordinates(filename)

# Change between specific timestep(file) or just the last one
LastTimeStep=True
if (LastTimeStep):
    files = glob(os.path.join(path,"step-*_wavefield.pbin"))
    w_filename= sorted(files)[-1]
else:
    w_filename = os.path.join(path,"step-2000_wavefield.pbin")

# Load wavefield file
se2_field = se2wave_load_wavefield(w_filename,True,True)

# Separate field components into matrices
LCoorX, LCoorY = SeparateList(se2_coor['coor'], se2_coor['nx'].item(), se2_coor['ny'].item())
LFieldX, LFieldY = Tilt_SeparateList(se2_field['vel'], se2_field['nx'].item(), se2_field['ny'].item(),-TiltAngle)

# Create the SPline function in a specific Field
SplineFunction = [RectBivariateSpline(LCoorX[:,0], LCoorY[0,:], LFieldX), 
                  RectBivariateSpline(LCoorX[:,0], LCoorY[0,:], LFieldY)]

# Get a profile between two coordinates using the SPline function 
ArrayDist, CompX, CompY = GetProfileData(LocIni,LocEnd,NumPoints, SplineFunction)

TimeTxt = "t = {}s".format(round(se2_field["time"].item(),5)) # Timestamp label

#BuildAndSaveDomainFig(LCoorX,LCoorY,LFieldX, LocIni, LocEnd, TimeTxt,
#                      "Displacement field Parallel-component [m]", 
#                      OutputFolder+"DispField_XComp.pdf")

#PlotProfileInter(ArrayDist, CompX, "Displacement field X-component [m]", 
#                "/home/nico/Documents/TEAR/Codes_TEAR/ProfilePicking/Output/Profile_XComp.pdf")


BuildAndSaveDomainFig(LCoorX,LCoorY,LFieldX, LocIni, LocEnd, TimeTxt,
                      "Velocity field Parallel-component [m/s]    ", 
                      OutputFolder+"VelField_XComp.pdf")

#PlotProfileInter(ArrayDist, CompX, "Velocity field Parallel-component [m]", 
#                 OutputFolder+ "Profile_XComp.pdf", delta)