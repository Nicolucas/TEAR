import os as os
from matplotlib import pyplot as plt
from scipy.interpolate import RectBivariateSpline
from Lib_ProfileFromImage import *
from se2waveload import *


LocIni,LocEnd = [-500, -500], [500, 500]
NumPoints = 1000

path = "/home/nico/Documents/TEAR/Codes_TEAR/plot-utils_se2wave/se2wave/"
filename = os.path.join(path,"default_mesh_coor.pbin")
se2_coor = se2wave_load_coordinates(filename);

w_filename = os.path.join(path,"step-1400_wavefield.pbin")
se2_field = se2wave_load_wavefield(w_filename,True,True);

LCoorX, LCoorY = SeparateList(se2_coor['coor'], se2_coor)
LFieldX, LFieldY = SeparateList(se2_field['displ'], se2_field)

SplineFunction = [RectBivariateSpline(LCoorX[:,0], LCoorY[0,:], LFieldX), 
                  RectBivariateSpline(LCoorX[:,0], LCoorY[0,:], LFieldY)]

ArrayDist, CompX, CompY = GetProfileData(LocIni,LocEnd,NumPoints, SplineFunction)



BuildAndSaveDomainFig(LCoorX,LCoorY,LFieldX, LocIni, LocEnd,
                      "Displacement field X-component [m]", "DispField_XComp.pdf")
BuildAndSaveDomainFig(LCoorX,LCoorY,LFieldY, LocIni, LocEnd,
                      "Displacement field Y-component [m]", "DispField_YComp.pdf")

PlotProfileInter(ArrayDist, CompX, "Displacement field X-component [m]", 
                "Profile_XComp.pdf")

PlotProfileInter(ArrayDist, CompY, "Displacement field Y-component [m]", 
                "Profile_YComp.pdf")