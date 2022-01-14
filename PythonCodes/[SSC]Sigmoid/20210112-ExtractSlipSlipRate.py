import sys
sys.path.insert(0,"/import/freenas-m-03-geodynamics/jhayek/petsc-3.12.5/lib/petsc/bin/")
sys.path.insert(0,"/import/freenas-m-03-geodynamics/jhayek/TEAR/se2wave/utils/python")
sys.path.insert(0,"/import/freenas-m-03-geodynamics/jhayek/TEAR/processing/TEAR/PythonCodes/LibFolder")

from Lib_GeneralFunctions import *
from Lib_SigmoidProcessing import *


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LocOfInterest = [0, 2000,4000,6000,8000]

dim = 25
Thickness = dim*2.501

NumPoints = 1200001
xx = np.linspace(-1.e4, 1.e4, NumPoints)

UniFolder = "TEAR2_Kos_Sig_P3_025x025_d2.5001_tanh12ph65"
fname = "step-{timestep:04}_wavefield.pbin"

path = "/import/freenas-m-03-geodynamics/jhayek/TEAR/Results/T2/Runs/{}/".format(UniFolder)

MaxTimeStep = 8180
freq = 10

OutFolder = "/import/freenas-m-03-geodynamics/jhayek/SharedWolfel/PaperData/Output/Sigmoid-{}/".format(GetTodayDate())
OutFile = "{UniFolder}-{dim}x{dim}-P3-".format(UniFolder = UniFolder, dim = dim)+str(Thickness)

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
