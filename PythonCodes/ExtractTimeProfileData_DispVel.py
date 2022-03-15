import os, sys, time
import numpy as np
from scipy.interpolate import RectBivariateSpline

sys.path.insert(0,"/home/nico/Documents/TEAR/Codes_TEAR/PythonCodes/LibFolder")
from Lib_GeneralFunctions import *
from Lib_ProfilePlotting import *
from Lib_ProfileProcessing import *

from se2waveload import *

def ExtractFieldsPerTS(ListTimeProfileObj, w_filename, se2_coor):
    se2_field = se2wave_load_wavefield(w_filename,True,True)
    TimeStep = se2_field["time"].item()

    LCoorX, LCoorY = SeparateList(se2_coor['coor'], se2_coor['nx'].item(), se2_coor['ny'].item())
    LFieldX, LFieldY = SeparateList(se2_field['displ'], se2_field['nx'].item(), se2_field['ny'].item())
    LFieldvelX, LFieldvelY = SeparateList(se2_field['vel'], se2_field['nx'].item(), se2_field['ny'].item())

    SplineDispl = [RectBivariateSpline(LCoorX[:,0], LCoorY[0,:], LFieldX, kx=1, ky=1), 
                    RectBivariateSpline(LCoorX[:,0], LCoorY[0,:], LFieldY, kx=1, ky=1)]
    SplineVel = [RectBivariateSpline(LCoorX[:,0], LCoorY[0,:], LFieldvelX, kx=1, ky=1), 
                    RectBivariateSpline(LCoorX[:,0], LCoorY[0,:], LFieldvelY, kx=1, ky=1)]
    for OBJitem in ListTimeProfileObj:
        CompDispX,CompDispY = GetOnlyLocData(OBJitem.Coord, SplineDispl, True)
        CompvelX,CompVelY = GetOnlyLocData(OBJitem.Coord, SplineVel, True)
        OBJitem.appendFieldValues(TimeStep, CompDispX, CompDispY, CompvelX, CompVelY)
 
def FillObjectInTime(ListTimeProfileObj, freq, maxtimestep, fname, path, MeshFilename = "default_mesh_coor.pbin"):
    TSList = np.arange(0, maxtimestep+1, freq).tolist()
    FilenameList = [os.path.join(path,fname.format(timestep=i)) for i in TSList]

    filename = os.path.join(path, MeshFilename)
    se2_coor = se2wave_load_coordinates(filename)

    [ExtractFieldsPerTS(ListTimeProfileObj, w_filename, se2_coor) for w_filename in FilenameList]



start_time = time.time()

freq = 1
maxtimestep =  901
thickness = 150.15

fname = "step-{timestep:04}_wavefield.pbin"
#path = "/home/nico/Documents/TEAR/Codes_TEAR/plot-utils_se2wave/se2wave/"
path = "/home/nico/Documents/TEAR/Codes_TEAR/se2dr/se2wave/"
#path = "/media/nico/Elements/Simulations/20200728/SSCdeg2/"
OutputFolder = "/home/nico/Documents/TEAR/Codes_TEAR/ProfilePicking/Output/" + GetTodayDate() + "/"
OutputFolder = "/home/nico/Documents/TEAR/Codes_TEAR/ProfilePicking/[TPV3]Results/" + GetTodayDate() + "/"


#Locations = [[8000,thickness],[6000,thickness],[4000,thickness],[2000,thickness],[0,thickness]]
Locations = [[12000,-3000], [12000,3000], [3000,thickness], [6000,thickness]]

ListTimeProfileObj = [SingleTimeProfile(Loc) for Loc in Locations]
FillObjectInTime(ListTimeProfileObj, freq, maxtimestep, fname, path)

SavePickleFile(OutputFolder, "TPList_t{timestep}_d{d}.pickle".format(timestep = maxtimestep, d = thickness), 
               ListTimeProfileObj)
print("--- %s seconds ---" % (time.time() - start_time))

