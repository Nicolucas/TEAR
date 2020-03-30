from matplotlib import pyplot as plt

from Lib_ProfileFromImage import *

from scipy.interpolate import RectBivariateSpline


filename = "/home/nico/Documents/TEAR/Codes_TEAR/plot-utils_se2wave/se2wave/default_mesh_coor.pbin"
se2_coor = se2wave_load_coordinates(filename);


w_filename = "/home/nico/Documents/TEAR/Codes_TEAR/plot-utils_se2wave/se2wave/step-1400_wavefield.pbin"
se2_field = se2wave_load_wavefield(w_filename,True,True);


LCoorX, LCoorY = SeparateList(se2_coor['coor'], se2_coor)
LFieldX, LFieldY = SeparateList(se2_field['displ'], se2_field)



SplineFunction = [RectBivariateSpline(LCoorX[:,0], LCoorY[0,:], LFieldX), 
                  RectBivariateSpline(LCoorX[:,0], LCoorY[0,:], LFieldY)]


PlotImage_GUI(LCoorX, LCoorY, LFieldY, SplineFunction)
