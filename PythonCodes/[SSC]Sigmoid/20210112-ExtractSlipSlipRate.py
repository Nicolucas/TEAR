import sys
#sys.path.insert(0,"/home/nico/Documents/TEAR/Codes_TEAR/PythonCodes/LibFolder")
sys.path.insert(0,"/import/freenas-m-03-geodynamics/jhayek/TEAR/processing/TEAR/PythonCodes/LibFolder")

from Lib_GeneralFunctions import *
from Lib_SigmoidProcessing import *


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LocOfInterest = [2000,4000,6000,8000]
Thickness = 50.05
xx = np.linspace(-1.e4, 1.e4, 6001)

fname = "step-{timestep:04}_wavefield.pbin"

path = "/import/freenas-m-03-geodynamics/jhayek/TEAR/Results/T5/se2wave/"

MaxTimeStep = 2261
freq = 10

OutFolder = "/import/freenas-m-03-geodynamics/jhayek/TEAR/processing/TEAR/PythonCodes/[SSC]Sigmoid/ProcessedData/"
OutFile = GetTodayDate()+"-T5-25x25-P1-"+str(Thickness)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TimesList = np.arange(0,MaxTimeStep+1, freq).tolist()

MeshFilename = os.path.join(path, "default_mesh_coor.pbin")
se2_coor = se2wave_load_coordinates(MeshFilename)

ListFaultDataObj = Init_ListFaultDataObj(Thickness, LocOfInterest, xx, func(xx), func_der(xx))

for i in TimesList:  
	FieldFilename = os.path.join(path,fname.format(timestep=i))

	ListFaultDataObj = PopulateListFaultDataObj_w_Fields(ListFaultDataObj, FieldFilename, se2_coor)

SavePickleFile(OutFolder,OutFile,ListFaultDataObj)
