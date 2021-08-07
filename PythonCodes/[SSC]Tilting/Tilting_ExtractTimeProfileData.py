import os, sys, time
import numpy as np
from scipy.interpolate import RectBivariateSpline

#sys.path.insert(0,"/import/freenas-m-03-geodynamics/jhayek/TEAR/processing/TEAR/PythonCodes/LibFolder")

sys.path.insert(0,"/home/nico/Tools/petsc-3.12.5/lib/petsc/bin/")
sys.path.insert(0,"/home/nico/Documents/TEAR/Codes_TEAR/se2dr/se2wave/utils/python/")
sys.path.insert(0,"/home/nico/Documents/TEAR/Codes_TEAR/PythonCodes/LibFolder")

from Lib_GeneralFunctions import *
from Lib_ProfilePlotting import *
from Lib_ProfileProcessing import *

from se2waveload import *

def ExtractFieldsPerTS_tilting(ListTimeProfileObj, w_filename, se2_coor,TiltAngle):
    se2_field = se2wave_load_wavefield(w_filename,True,True)
    TimeStep = se2_field["time"].item()

    LCoorX, LCoorY = SeparateList(se2_coor['coor'], se2_coor['nx'].item(), se2_coor['ny'].item())
    LFieldX, LFieldY = Tilt_SeparateList(se2_field['displ'], se2_field['nx'].item(), se2_field['ny'].item(), TiltAngle)
    LFieldvelX, LFieldvelY = Tilt_SeparateList(se2_field['vel'], se2_field['nx'].item(), se2_field['ny'].item(), TiltAngle)

    SplineDispl = [RectBivariateSpline(LCoorX[:,0], LCoorY[0,:], LFieldX, kx=1, ky=1), 
                    RectBivariateSpline(LCoorX[:,0], LCoorY[0,:], LFieldY, kx=1, ky=1)]
    SplineVel = [RectBivariateSpline(LCoorX[:,0], LCoorY[0,:], LFieldvelX, kx=1, ky=1), 
                    RectBivariateSpline(LCoorX[:,0], LCoorY[0,:], LFieldvelY, kx=1, ky=1)]
    for OBJitem in ListTimeProfileObj:
        CompDispX,CompDispY = GetLocDataTilt(OBJitem.Coord,OBJitem.TwinCoord, SplineDispl, True)
        CompvelX,CompVelY = GetLocDataTilt(OBJitem.Coord,OBJitem.TwinCoord, SplineVel,  True)
        OBJitem.appendFieldValues(TimeStep, CompDispX, CompDispY, CompvelX, CompVelY)
 
def FillObjectInTime(ListTimeProfileObj, freq, maxtimestep, fname, path, TiltAngle, MeshFilename = "default_mesh_coor.pbin"):
    TSList = np.arange(0, maxtimestep+1, freq).tolist()
    FilenameList = [os.path.join(path,fname.format(timestep=i)) for i in TSList]

    filename = os.path.join(path, MeshFilename)
    se2_coor = se2wave_load_coordinates(filename)

    [ExtractFieldsPerTS_tilting(ListTimeProfileObj, w_filename, se2_coor,TiltAngle) for w_filename in FilenameList]



start_time = time.time()

#%%%%%%%%%%%%%%%%%%%%%%% Parameter Input
freq = 10
maxtimestep = 2040 
thickness = 100 #25.0*1.001
TiltAngle = 20.00
InFolder = "P6_P3_100x100_Dg"

fname = "step-{timestep:04}_wavefield.pbin"
#path = "/import/freenas-m-03-geodynamics/jhayek/TEAR/Results/{}/se2wave/".format(InFolder)
path = "/home/nico/Documents/TEAR/Codes_TEAR/se2dr/Runs/{}".format(InFolder)
path = "/media/nico/Elements/Simulations/20210730-SSCTilt20-100delta/{}".format(InFolder)

OutputFolder = "/home/nico/Documents/TEAR/Codes_TEAR/PythonCodes/[SSC]Sigmoid/ProcessedData/" + GetTodayDate() + "-Tilting/"
#%%%%%%%%%%%%%%%%%%%%%%% Parameter Input

OutFileName = "{InFolder}-PdepBlend-TPList_t{timestep}_d{d}.pickle".format(InFolder=InFolder, timestep = maxtimestep, d = thickness)
print("\nSTART: "+OutFileName)

# Locations relative to the fault
Locations = [[8000,thickness],[6000,thickness],[4000,thickness],[2000,thickness],[0,thickness]]

#Locations in the real domain
TwinLocations = [list(ApplyTilting(TiltAngle,Loc[0],-Loc[1])) for Loc in Locations]
Locations = [list(ApplyTilting(TiltAngle,Loc[0],Loc[1])) for Loc in Locations]

ListTimeProfileObj = [SingleTimeProfile(Loc) for Loc in Locations]
[ListTimeProfileObj[idx].AddTwin(TLoc) for idx,TLoc in enumerate(TwinLocations)]
FillObjectInTime(ListTimeProfileObj, freq, maxtimestep, fname, path, -TiltAngle)

SavePickleFile(OutputFolder, OutFileName, ListTimeProfileObj)
print("--- %s seconds ---" % (time.time() - start_time))

