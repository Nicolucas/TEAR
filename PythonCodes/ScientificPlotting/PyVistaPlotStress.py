import matplotlib.pyplot as plt
import pyvista as pv
import numpy as np
import json
import glob

import os, sys, time

sys.path.insert(0,"/import/freenas-m-03-geodynamics/jhayek/TEAR/processing/TEAR/PythonCodes/LibFolder")
from Lib_GeneralFunctions import *
from Lib_PyVista import *


folder='/import/freenas-m-03-geodynamics/jhayek/TEAR/Results/T2/Runs/'

RunTest='TEAR51_TPV_T0_P3_025x025_A12phi65_Delta1.001_3s_NC/'
#RunTest='TEAR40_Kos_Sig_P3_025x025_A12phi65_Delta2.5_4s/'
#RunTest='TEAR40_Kos_T20_P3_025x025_A12phi65_Delta2.5_4s/'



TimeStep = 0800
#TimeStep = 5210 #TPV3

Stressfilename= 'Sigma-Aligned-step-{TimeStep:04d}.vtu'.format(TimeStep=TimeStep)
JSONfilename= 'step-{TimeStep:04d}_wavefield.json'.format(TimeStep=TimeStep)


Stressfilename= 'Sigma-Aligned-step-{TimeStep:04d}.vtu'.format(TimeStep=TimeStep)
JSONfilename= 'step-{TimeStep:04d}_wavefield.json'.format(TimeStep=TimeStep)


OutFolder = folder + RunTest + 'Out/'
CreateFolder(OutFolder)

print(JSONfilename)
with open(folder+RunTest+JSONfilename, 'r') as json_file:
  LoadedJson = json.load(json_file)

reader = pv.get_reader(folder+RunTest+Stressfilename)
StressMesh = reader.read()
data = StressMesh.get_array('sigma_xy',preference='point')
