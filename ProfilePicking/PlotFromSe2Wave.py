import importlib.util
spec = importlib.util.spec_from_file_location("PetscBinaryIO", "/home/nico/Documents/TEAR/Codes_TEAR/plot-utils_se2wave/se2wave/utils/python/PetscBinaryIO.py")
pio = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pio)

import numpy as np
from matplotlib import pyplot as plt
import math

from scipy.interpolate import RectBivariateSpline


def se2wave_load_coordinates(filename):
  debug = True
  io = pio.PetscBinaryIO()
  data = dict()
  with open(filename) as fp:
    v = io.readInteger(fp)
    data['mx'] = v
    
    v = io.readInteger(fp)
    data['my'] = v
    
    v = io.readInteger(fp)
    data['nx'] = v
    
    v = io.readInteger(fp)
    data['ny'] = v
    
    objecttype = io.readObjectType(fp)
    v = io.readVec(fp)
    data['coor'] = v
  
  if debug:
    print('#elements-x',data['mx'])
    print('#elements-y',data['my'])
    print('#basis-x',data['nx'])
    print('#basis-y',data['ny'])

  return data

def se2wave_load_wavefield(filename,has_displacement,has_velocity):
  debug = True
  io = pio.PetscBinaryIO()
  data = dict()
  with open(filename) as fp:
    v = io.readInteger(fp)
    data['mx'] = v
    
    v = io.readInteger(fp)
    data['my'] = v
    
    v = io.readInteger(fp)
    data['nx'] = v
    
    v = io.readInteger(fp)
    data['ny'] = v

    v = io.readInteger(fp)
    data['step'] = v

    v = io.readReal(fp)
    data['time'] = v

    if has_displacement:
      objecttype = io.readObjectType(fp)
      v = io.readVec(fp)
      data['displ'] = v

    if has_velocity:
      objecttype = io.readObjectType(fp)
      v = io.readVec(fp)
      data['vel'] = v

  if debug:
    print('#elements-x',data['mx'])
    print('#elements-y',data['my'])
    print('#basis-x',data['nx'])
    print('#basis-y',data['ny'])
    print('step',data['step'])
    print('time',data['time'])
  
  return data

def SeparateList(List2Sep,data):
    TotNum = len(List2Sep)
    xComponent = [List2Sep[i] for i in range(TotNum) if i%2==1]
    yComponent = [List2Sep[i] for i in range(TotNum) if i%2==0]
    
    xComponent = np.reshape(xComponent,(data['nx'][0],data['ny'][0]))
    yComponent = np.reshape(yComponent,(data['nx'][0],data['ny'][0]))
    return xComponent,yComponent

def GetProfileData(LocIni,LocEnd,NumPoints, SplineFunction):
    x0,y0,x1,y1 = LocIni[0], LocIni[1], LocEnd[0], LocEnd[1]
    x, y = np.linspace(x0, x1, NumPoints),   np.linspace(y0, y1, NumPoints)

    CompX = [SplineFunction[0](x[i], y[i])[0][0] for i in range(NumPoints)]
    CompY = [SplineFunction[1](x[i], y[i])[0][0] for i in range(NumPoints)]

    Dist = math.sqrt((x1 - x0)**2.0 + (y1 - y0)**2.0)
    ArrayDist = np.linspace(0, Dist, NumPoints)

    return ArrayDist, CompX, CompY

def FigSetupAndPlotImage(LCoorX,LCoorY,Field, SplineFunction):
    fig = plt.figure(figsize = (10,5) ,constrained_layout=True)
    gs = fig.add_gridspec(2, 4)
    ax1 = fig.add_subplot(gs[0:2, 0:2])
    ax1.set_aspect('equal', 'box')
    axList= [fig.add_subplot(gs[i,2:4]) for i in range(2)]

    img = ax1.pcolormesh(LCoorX,LCoorY,Field)

    return img, axList

def PlotProfiles(ArrayDist, CompX, CompY, axList):
    axList[0].plot(ArrayDist,CompX)
    axList[1].plot(ArrayDist,CompY)

def PlotLocLine(img,LocIni,LocEnd):
    img.axes.annotate("",
                    xy = (LocIni[0], LocIni[1]),
                    xycoords = 'data',
                    xytext = (LocEnd[0], LocEnd[1]), 
                    textcoords = 'data',
                    arrowprops = dict(
                        arrowstyle = "<-",
                        connectionstyle = "arc3", 
                        color = 'white',
                        alpha = 0.7,
                        linewidth = 3
                        ),
                    )


LocIni,LocEnd = [-500, -500], [500, 500]
NumPoints = 1000

filename = "/home/nico/Documents/TEAR/Codes_TEAR/plot-utils_se2wave/se2wave/default_mesh_coor.pbin"
se2_coor = se2wave_load_coordinates(filename);

w_filename = "/home/nico/Documents/TEAR/Codes_TEAR/plot-utils_se2wave/se2wave/step-1400_wavefield.pbin"
se2_field = se2wave_load_wavefield(w_filename,True,True);

LCoorX, LCoorY = SeparateList(se2_coor['coor'], se2_coor)
LFieldX, LFieldY = SeparateList(se2_field['displ'], se2_field)

SplineFunction = [RectBivariateSpline(LCoorX[:,0], LCoorY[0,:], LFieldX), 
                  RectBivariateSpline(LCoorX[:,0], LCoorY[0,:], LFieldY)]

ArrayDist, CompX, CompY = GetProfileData(LocIni,LocEnd,NumPoints, SplineFunction)

img, axList = FigSetupAndPlotImage(LCoorX, LCoorY, LFieldX, SplineFunction)

PlotLocLine(img,LocIni,LocEnd)

PlotProfiles(ArrayDist, CompX, CompY, axList)
plt.show()



