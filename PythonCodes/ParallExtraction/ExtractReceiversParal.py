import os, time, sys
import numpy as np
from scipy.interpolate import RectBivariateSpline
from functools import partial
import copy
import multiprocessing as mp

sys.path.insert(0,"/import/freenas-m-03-geodynamics/jhayek/petsc-3.12.5/lib/petsc/bin/")
sys.path.insert(0,"/import/freenas-m-03-geodynamics/jhayek/TEAR/se2wave/utils/python")
sys.path.insert(0,"/import/freenas-m-03-geodynamics/jhayek/TEAR/processing/TEAR/PythonCodes/LibFolder")


from Lib_GeneralFunctions import *
from Lib_ProfilePlotting import *
from Lib_ProfileProcessing import *

from se2waveload import *

def ExtractFieldsPerTS(w_filename, ListTimeProfileObj,  se2_coor):
    se2_field = se2wave_load_wavefield(w_filename,True,True)
    TimeStep = se2_field["time"].item()

    LCoorX, LCoorY = SeparateList(se2_coor['coor'], se2_coor['nx'].item(), se2_coor['ny'].item())
    LFieldX, LFieldY = SeparateList(se2_field['displ'], se2_field['nx'].item(), se2_field['ny'].item())
    LFieldvelX, LFieldvelY = SeparateList(se2_field['vel'], se2_field['nx'].item(), se2_field['ny'].item())

    SplineDispl = [RectBivariateSpline(LCoorX[:,0], LCoorY[0,:], LFieldX, kx=1, ky=1), 
                    RectBivariateSpline(LCoorX[:,0], LCoorY[0,:], LFieldY, kx=1, ky=1)]
    SplineVel = [RectBivariateSpline(LCoorX[:,0], LCoorY[0,:], LFieldvelX, kx=1, ky=1), 
                    RectBivariateSpline(LCoorX[:,0], LCoorY[0,:], LFieldvelY, kx=1, ky=1)]

    Child_ListTimeProfileObj = copy.deepcopy(ListTimeProfileObj)
    for OBJitem in Child_ListTimeProfileObj:
        CompDispX,CompDispY = GetOnlyLocData(OBJitem.Coord, SplineDispl)
        CompvelX,CompVelY = GetOnlyLocData(OBJitem.Coord, SplineVel)
        OBJitem.appendFieldValues(TimeStep, CompDispX, CompDispY, CompvelX, CompVelY)

    print("t:{} - List Receivers Obtained".format(TimeStep))

    return Child_ListTimeProfileObj
 
def FillObjectInTime(ListTimeProfileObj, freq, maxtimestep, fname, path, MeshFilename = "default_mesh_coor.pbin", NumProc = 1):
    TSList = np.arange(0, maxtimestep+1, freq).tolist()
    FilenameList = [os.path.join(path,fname.format(timestep=i)) for i in TSList]

    filename = os.path.join(path, MeshFilename)
    se2_coor = se2wave_load_coordinates(filename)

    

    ExtractFieldsPerTS_Partial = partial(ExtractFieldsPerTS, ListTimeProfileObj = ListTimeProfileObj, se2_coor= se2_coor)

    proc_Results = []

    with mp.Pool(processes=NumProc) as pool:
        processes = tqdm(pool.map_async(ExtractFieldsPerTS_Partial, FilenameList))
        proc_Results = processes.get()

    print(proc_Results)

    for proc in proc_Results:
        [ListTimeProfileObj[idx_rec].AssimilateSingleTimeProfile(SingleObj) for idx_rec, SingleObj in enumerate(proc)]




##################################################################################################################################

start_time = time.time()

##########################################
NProc = 64
ThickVal = "025"
thickness = float(ThickVal)*1.001

InFolder = "TEAR49_TPV_T0_P3_025x025_A12phi65_Delt1.001_7s"


fname = "step-{timestep:04}_wavefield.pbin"
NameWrapper = "{}/".format(InFolder)
path = "/import/freenas-m-03-geodynamics/jhayek/TEAR/Results/T2/Runs/{}".format(NameWrapper)

TimeStepList = GetListPatternFiles(path,fname,"{timestep:04}")


freq = int(TimeStepList[1])-int(TimeStepList[0])
maxtimestep = int(TimeStepList[-1]) 

OutputFolder = "/import/freenas-m-03-geodynamics/jhayek/SharedWolfel/PaperData/ReceiverData/" + GetTodayDate() + "Paral/"

OutFileName = "Receivers_{InFolder}_{timestep}_d{d}.pickle".format(InFolder=InFolder, timestep = maxtimestep, d = thickness)

#############################

print("\n>>START: "+OutFileName+"\n")

# Locations
Locations = [[0,thickness/2.0],[2000,thickness/2.0],[4000,thickness/2.0],[6000,thickness/2.0],[8000,thickness/2.0],
             [0,25],[2000,25],[4000,25],[6000,25],[8000,25],
             [0,50],[2000,50],[4000,50],[6000,50],[8000,50],
             [0,100],[2000,100],[4000,100],[6000,100],[8000,100],
             [0,200],[2000,200],[4000,200],[6000,200],[8000,200],
             [0,300],[2000,300],[4000,300],[6000,300],[8000,300],
             [0,400],[2000,400],[4000,400],[6000,400],[8000,400],
             [0,500],[2000,500],[4000,500],[6000,500],[8000,500],
            ]


ListTimeProfileObj = [SingleTimeProfile(Loc) for Loc in Locations]
FillObjectInTime(ListTimeProfileObj, freq, maxtimestep, fname, path, NumProc=NProc)

SavePickleFile(OutputFolder, OutFileName, ListTimeProfileObj)
print("--- %s seconds ---" % (time.time() - start_time))
