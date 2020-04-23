import os
import numpy as np
from matplotlib import pyplot as plt
from scipy.interpolate import RectBivariateSpline
from Lib_ProfileFromImage import *
from se2waveload import *

def ExtractCompPerTS(Loc, w_filename, se2_coor):
    se2_field = se2wave_load_wavefield(w_filename,True,True);

    LCoorX, LCoorY = SeparateList(se2_coor['coor'], se2_coor['nx'].item(), se2_coor['ny'].item())
    LFieldX, LFieldY = SeparateList(se2_field['displ'], se2_field['nx'].item(), se2_field['ny'].item())

    SplineFunction = [RectBivariateSpline(LCoorX[:,0], LCoorY[0,:], LFieldX), 
                    RectBivariateSpline(LCoorX[:,0], LCoorY[0,:], LFieldY)]

    CompX,CompY = GetLocData(Loc, SplineFunction)
    return [CompX,CompY]


def ExtracNPlotInTime(Loc, freq, maxtimestep, fname, path):
    TSList = np.arange(0, maxtimestep+1, freq).tolist()
    FilenameList = [os.path.join(path,fname.format(timestep=i)) for i in TSList]

    filename = os.path.join(path,"default_mesh_coor.pbin")
    se2_coor = se2wave_load_coordinates(filename)

    Fieldvalue = [ExtractCompPerTS(Loc, w_filename, se2_coor) for w_filename in FilenameList]

    CompX, CompY = zip(*Fieldvalue)

    PlotTimeProfile(TSList, CompX, "Time Profile X-Component @ ({x},{y})".format(x=Loc[0],y=Loc[1]), "{x}-{y}-XComp-TimeProfile.pdf".format(x=Loc[0],y=Loc[1]))
    PlotTimeProfile(TSList, CompY, "Time Profile Y-Component @ ({x},{y})".format(x=Loc[0],y=Loc[1]), "{x}-{y}-YComp-TimeProfile.pdf".format(x=Loc[0],y=Loc[1]))

    return CompX,CompY

# Timesteps from the simulation and output filename
freq = 1
maxtimestep = 241
fname = "step-{timestep:04}_wavefield.pbin"

path = "/home/nico/Documents/TEAR/Codes_TEAR/plot-utils_se2wave/se2wave/"

Loc = [[0,220],[0,200],[0,180],[0,150],[0,50],[0,0]]

CompList = [ExtracNPlotInTime(LocI, freq, maxtimestep, fname, path) for LocI in Loc]
CompXList, CompYList = zip(*CompList)

TSList = np.arange(0, maxtimestep+1, freq).tolist()
PlotTimeProfile(TSList, CompXList, "Time Profile X-Component - Combined", "Combined-XComp-TimeProfile.pdf",True,Loc)
PlotTimeProfile(TSList, CompYList, "Time Profile Y-Component - Combined", "Combined-YComp-TimeProfile.pdf",True,Loc)