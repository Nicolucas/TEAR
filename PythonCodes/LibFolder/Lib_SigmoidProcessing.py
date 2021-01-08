
import os, sys, time
import numpy as np
from scipy.interpolate import RectBivariateSpline
from sklearn.metrics.pairwise import euclidean_distances

import matplotlib.pylab as plt

from se2waveload import *


## Zero level set definition
# Sigmoid or any function of interest to represent the center of the fault / Zero level set function
def func(x, k=-0.0002, amp = 2.0):
    fx = amp * (x - x * k) / (k - abs(x) * 2.0 * k + 1.0)
    return fx

# The respective derivative ofthe previous zero level set function
def func_der(x, k=-0.0002, amp = 2.0):
    fx_prime = amp * (1 - k * k) / ((k - abs(x) * 2.0 * k + 1.0)*(k - abs(x) * 2.0 * k + 1.0))
    return fx_prime


# Functions from Lib_ProfileProcessing adapted for sigmoid use and more general
def SeparateList(List2Sep,nx,ny):
    TotNum = len(List2Sep)
    xComponent = List2Sep[0:TotNum:2]
    yComponent = List2Sep[1:TotNum:2]

    xComponent = np.reshape(xComponent, (nx, ny), "F")
    yComponent = np.reshape(yComponent, (nx, ny), "F")
    return xComponent,yComponent


def ExtractFields(w_filename, se2_coor):
    se2_field = se2wave_load_wavefield(w_filename,True,True)
    TimeStep = se2_field["time"].item()

    LCoorX, LCoorY         = SeparateList(se2_coor['coor'],   se2_coor['nx'].item(),  se2_coor['ny'].item())
    LFieldX, LFieldY       = SeparateList(se2_field['displ'], se2_field['nx'].item(), se2_field['ny'].item())
    LFieldvelX, LFieldvelY = SeparateList(se2_field['vel'],   se2_field['nx'].item(), se2_field['ny'].item())
    
    return TimeStep, LCoorX, LCoorY, LFieldX, LFieldY, LFieldvelX, LFieldvelY


def GetBivariateSplineFuncFromFields(LCoorX, LCoorY, LFieldX, LFieldY,LFieldvelX, LFieldvelY):
    SplineDispl = [RectBivariateSpline(LCoorX[:,0], LCoorY[0,:], LFieldX, kx=1, ky=1), 
                   RectBivariateSpline(LCoorX[:,0], LCoorY[0,:], LFieldY, kx=1, ky=1)]
    SplineVel = [RectBivariateSpline(LCoorX[:,0], LCoorY[0,:], LFieldvelX, kx=1, ky=1), 
                 RectBivariateSpline(LCoorX[:,0], LCoorY[0,:], LFieldvelY, kx=1, ky=1)]
    
    return SplineDispl, SplineVel


def GetLocData(Loc, SplineFunction, GetSlip=False):

    CompX = SplineFunction[0](Loc[0],Loc[1])[0][0]
    CompY = SplineFunction[1](Loc[0],Loc[1])[0][0]
    return CompX, CompY


# Wrappers
def GetSplineFunctions(w_filename, se2_coor):
    TimeStepVal, LCoorX, LCoorY, LFieldX, LFieldY, LFieldvelX, LFieldvelY =  ExtractFields(w_filename, se2_coor)
    SplineDisplPair, SplineVelPair = GetBivariateSplineFuncFromFields(LCoorX, LCoorY, 
                                                              LFieldX, LFieldY,
                                                              LFieldvelX, LFieldvelY)
    return TimeStepVal, SplineDisplPair, SplineVelPair



def GetListCoords_FaultGeometry(xx):
    return list(map(list,zip(xx,func(xx))))

def GetClosestIndex(Coordx, Coordy, FaultNodes_CoordList):
    IndexMin = np.argmin(euclidean_distances(FaultNodes_CoordList,[[Coordx,Coordy]]))
    return IndexMin
    
def CreateMeshGrid(Nx,Ny):
    x = np.linspace(0, 1, Nx)
    y = np.linspace(0, 1, Ny)
    
    xv, yv = np.meshgrid(x, y)
    return xv, yv


def CumulativeDistOnFault(idx:int, xArray, centerIdx:int):
    Distance = 0.0
    if (idx > centerIdx):
        for i in range(centerIdx,idx):
            Distance += euclidean_distances([[xArray[i+1],func(xArray[i+1])]],
                                            [[xArray[i],func(xArray[i])]]
                                           )[0][0]
        
    elif (idx < centerIdx):
        for i in range(idx,centerIdx):
            Distance += euclidean_distances([[xArray[i+1],func(xArray[i+1])]],
                                            [[xArray[i],func(xArray[i])]]
                                           )[0][0]
    return Distance

def get_SDF(xCoord,yCoord,idx,xArray):
    return euclidean_distances([[xCoord,yCoord]],[[xArray[idx],func(xArray[idx])]])[0][0]

def FindDistIndx(DistOfInterest, xArray, centerIdx:int):
    Distance = 0.0
    for i in range(centerIdx,len(xArray)):
        Distance += euclidean_distances([[xArray[i+1],func(xArray[i+1])]],
                                        [[xArray[i],func(xArray[i])]]
                                        )[0][0]
        if(Distance>DistOfInterest):
            break
    return i
 

class FaultData:
    def __init__(self, LocOfInterest, HalfThickness, Xval, Fxval, FxPrimeVal):
        self.DistAlongFault = LocOfInterest
        self.HalfThickness = HalfThickness
        self.Xval = Xval
        self.Fxval = Fxval
        self.FxPrimeVal = FxPrimeVal
        
        self.Normal = np.array(self.NormalVector(self.FxPrimeVal))
        self.Tangent = np.array(self.TangentVector(self.FxPrimeVal))
        
        self.RecSide1X = self.Xval + (np.array(self.NormalVector(self.FxPrimeVal))*HalfThickness)[0]
        self.RecSide1Y = self.Fxval + (np.array(self.NormalVector(self.FxPrimeVal))*HalfThickness)[1]
        self.RecSide2X = self.Xval - (np.array(self.NormalVector(self.FxPrimeVal))*HalfThickness)[0]
        self.RecSide2Y = self.Fxval - (np.array(self.NormalVector(self.FxPrimeVal))*HalfThickness)[1]
        
        
        self.Time = []
        self.Slip = []
        self.SlipRate = []
        
    def __repr__(self):
        return "FaultData Object, distance {} - half thickness {}".format(self.DistAlongFault, self.HalfThickness)
    
    def __str__(self):
        return "Fault Data at distance {}, half thickness {}".format(self.DistAlongFault, self.HalfThickness)
    
    def GetReceiverCoords(self):
        return [self.RecSide1X, self.RecSide1Y]
    
    def GetTwinReceiverCoords(self):
        return [self.RecSide2X, self.RecSide2Y]
    
    def appendFaultValues(self, time, Slip, SlipRate):
        self.Time.append(time)
        self.Slip.append(Slip)
        self.SlipRate.append(SlipRate)
        
    def ExtractTangentFieldComponentDiff(self, SplineFunctionPair):
        SetCoords1 = [self.RecSide1X, self.RecSide1Y]
        SetCoords2 = [self.RecSide2X, self.RecSide2Y]
        
        Comp1X, Comp1Y = GetLocData(SetCoords1, SplineFunctionPair)
        Comp2X, Comp2Y = GetLocData(SetCoords2, SplineFunctionPair)
        
        TanDisp1 = self.Tangent[0]*Comp1X + self.Tangent[1]*Comp1Y 
        TanDisp2 = self.Tangent[0]*Comp2X + self.Tangent[1]*Comp2Y 
        
        return TanDisp1 - TanDisp2
        
    
    # Tangent vector for a given derivative
    def TangentVector(self, fPrimeX, **kwargs):
        mag = np.sqrt(1.0 + fPrimeX * fPrimeX)

        TangentX = 1.0/mag
        TangentY = fPrimeX/mag
        return TangentX, TangentY

    # Normal vector for a given derivative
    def NormalVector(self, fPrimeX, **kwargs):
        mag = np.sqrt(1.0 + fPrimeX * fPrimeX)

        NormalX = -fPrimeX/mag
        NormalY = 1.0/mag

        return NormalX, NormalY

def TransposeListOfCoordinates(List2D):
    return np.array(List2D).T.tolist()

def Init_ListFaultDataObj(Thickness, LocOfInterest, XArray, FXArray, FXPrimeArray, centerIdx = 3000):
    idx_Interest = [FindDistIndx(DistOfInterest = a, xArray = XArray, centerIdx = centerIdx) for a in LocOfInterest]

    ListFaultDataObj = [FaultData(LocOfInterest[i], Thickness, XArray[idx_Interest[i]], FXArray[idx_Interest[i]], 
                                  FXPrimeArray[idx_Interest[i]]) for i in range(len(idx_Interest))]
    
    return ListFaultDataObj

def PopulateListFaultDataObj_w_Fields(ListFaultDataObj, w_filename, se2_coor):
    TimeStep, SplineDisplPair, SplineVelPair = GetSplineFunctions(w_filename, se2_coor)
    
    for FaultDataObj in ListFaultDataObj:
        Slip = FaultDataObj.ExtractTangentFieldComponentDiff(SplineDisplPair)
        SlipRate = FaultDataObj.ExtractTangentFieldComponentDiff(SplineVelPair)
        
        FaultDataObj.appendFaultValues(TimeStep, Slip, SlipRate)
        
    return ListFaultDataObj