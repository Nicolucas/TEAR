import sys
#sys.path.insert(0,"/home/nico/Documents/TEAR/Codes_TEAR/PythonCodes/LibFolder")
sys.path.insert(0,"/import/freenas-m-03-geodynamics/jhayek/TEAR/processing/TEAR/PythonCodes/LibFolder")

from Lib_GeneralFunctions import *
from Lib_SigmoidProcessing import *


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LocOfInterest = [2000,4000,6000,8000]

dim = 25
Thickness = dim*2.001

NumPoints = 1200001
xx = np.linspace(-1.e4, 1.e4, NumPoints)

UniFolder = "T1"
fname = "step-{timestep:04}_wavefield.pbin"

path = "/import/freenas-m-03-geodynamics/jhayek/TEAR/Results/{}/se2wave/".format(UniFolder)

MaxTimeStep = 55400
freq = 100

OutFolder = "/import/freenas-m-03-geodynamics/jhayek/TEAR/processing/TEAR/PythonCodes/[SSC]Sigmoid/ProcessedData/"
OutFile = GetTodayDate()+"-{UniFolder}-{dim}x{dim}-P3-".format(UniFolder = UniFolder, dim = dim)+str(Thickness)

print("\n>>"+OutFile)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TimesList = np.arange(0,MaxTimeStep+1, freq).tolist()

MeshFilename = os.path.join(path, "default_mesh_coor.pbin")
se2_coor = se2wave_load_coordinates(MeshFilename)

ListFaultDataObj = Init_ListFaultDataObj(Thickness, LocOfInterest, xx, func(xx), func_der(xx), centerIdx = int((NumPoints-1)/2.0))

for i in TimesList:  
	FieldFilename = os.path.join(path,fname.format(timestep=i))

	ListFaultDataObj = PopulateListFaultDataObj_w_Fields(ListFaultDataObj, FieldFilename, se2_coor)

SavePickleFile(OutFolder,OutFile,ListFaultDataObj)
