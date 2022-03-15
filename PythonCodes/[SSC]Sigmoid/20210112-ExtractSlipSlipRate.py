import sys

from Lib_GeneralFunctions import *
from Lib_SigmoidProcessing import *


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LocOfInterest = [0, 2000,4000,6000,8000]

dim = 25
Thickness = dim*2.0

NumPoints = 1200001
xx = np.linspace(-1.e4, 1.e4, NumPoints)

UniFolder = "TEAR35_Kos_Sig_P3_025x025_A12phi65_Delta2_4s"
fname = "step-{timestep:04}_wavefield.pbin"

path = "/import/freenas-m-03-geodynamics/jhayek/TEAR/Results/T2/Runs/{}/".format(UniFolder)


OutFolder = "/import/freenas-m-03-geodynamics/jhayek/SharedWolfel/PaperData/CorrectedSimulations/Sigmoid-{}/".format(GetTodayDate())
OutFile = "{UniFolder}-{dim}x{dim}-P3-".format(UniFolder = UniFolder, dim = dim)+str(Thickness)

print("\n>>"+OutFile)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TimeStepList = GetListPatternFiles(path,fname,"{timestep:04}")


freq = int(TimeStepList[1])-int(TimeStepList[0])
MaxTimeStep = int(TimeStepList[-1]) 


MeshFilename = os.path.join(path, "default_mesh_coor.pbin")
se2_coor = se2wave_load_coordinates(MeshFilename)

ListFaultDataObj = Init_ListFaultDataObj(Thickness, LocOfInterest, xx, func(xx), func_der(xx), centerIdx = int((NumPoints-1)/2.0))

print("\n>>START: "+OutFile+"\n")

for j in TimeStepList:  
	i=int(j)
	FieldFilename = os.path.join(path,fname.format(timestep=i))

	ListFaultDataObj = PopulateListFaultDataObj_w_Fields(ListFaultDataObj, FieldFilename, se2_coor)

SavePickleFile(OutFolder,OutFile,ListFaultDataObj)
