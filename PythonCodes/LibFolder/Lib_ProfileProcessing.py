import os
import math
import numpy as np

def ApplyTilting(angle,xComp,yComp):
    NewXComp = np.array(xComp) * np.cos(angle * np.pi / 180.0) - np.array(yComp) * np.sin(angle * np.pi / 180.0) 
    NewYComp = np.array(xComp) * np.sin(angle * np.pi / 180.0) + np.array(yComp) * np.cos(angle * np.pi / 180.0) 
    return NewXComp, NewYComp
    

def Tilt_SeparateList(List2Sep,nx,ny,angle):
    TotNum = len(List2Sep)
    xComponent = List2Sep[0:TotNum:2]
    yComponent = List2Sep[1:TotNum:2]
    
    xComponent,yComponent = ApplyTilting(angle,xComponent,yComponent)
    
    xComponent = np.reshape(xComponent, (nx, ny), "F")
    yComponent = np.reshape(yComponent, (nx, ny), "F")
    return xComponent,yComponent

def GetLocDataTilt(Loc, TwinLoc, SplineFunction, GetSlip=False):
    x0,y0 = Loc[0],Loc[1]

    CompX = SplineFunction[0](Loc[0],Loc[1])[0][0]
    CompY = SplineFunction[1](Loc[0],Loc[1])[0][0]
    if (GetSlip and (Loc[1]!=0)): 
        TwinX = SplineFunction[0](TwinLoc[0],TwinLoc[1])[0][0]
        TwinY = SplineFunction[1](TwinLoc[0],TwinLoc[1])[0][0]
        CompX=CompX-TwinX
        CompY=CompY-TwinY

    return CompX, CompY

def SeparateList(List2Sep,nx,ny):
    TotNum = len(List2Sep)
    xComponent = List2Sep[0:TotNum:2]
    yComponent = List2Sep[1:TotNum:2]

    xComponent = np.reshape(xComponent, (nx, ny), "F")
    yComponent = np.reshape(yComponent, (nx, ny), "F")
    return xComponent,yComponent

def GetProfileData(LocIni,LocEnd,NumPoints, SplineFunction):
    x0,y0,x1,y1 = LocIni[0], LocIni[1], LocEnd[0], LocEnd[1]
    x, y = np.linspace(x0, x1, NumPoints),   np.linspace(y0, y1, NumPoints)

    CompX = [SplineFunction[0](x[i], y[i])[0][0] for i in range(NumPoints)]
    CompY = [SplineFunction[1](x[i], y[i])[0][0] for i in range(NumPoints)]

    Dist = math.sqrt((x1 - x0)**2.0 + (y1 - y0)**2.0)
    ArrayDist = np.linspace(0, Dist, NumPoints)

    return ArrayDist, CompX, CompY

def GetLocData(Loc, SplineFunction, GetSlip=False):

    CompX = SplineFunction[0](Loc[0],Loc[1])[0][0]
    CompY = SplineFunction[1](Loc[0],Loc[1])[0][0]
    if (GetSlip and (Loc[1]!=0)): #dirty fix for getting slip
        TwinX = SplineFunction[0](Loc[0],-Loc[1])[0][0]
        TwinY = SplineFunction[1](Loc[0],-Loc[1])[0][0]
        CompX=CompX-TwinX
        CompY=CompY-TwinY

    return CompX, CompY

def GetOnlyLocData(Loc, SplineFunction):
    CompX = SplineFunction[0](Loc[0],Loc[1])[0][0]
    CompY = SplineFunction[1](Loc[0],Loc[1])[0][0]
    
    return CompX, CompY





class SingleTimeProfile:
    def __init__(self,Coord):
        self.Coord = Coord
        self.Time = []
        self.DispX = []
        self.DispY = []
        self.VelX = []
        self.VelY = []
        self.TwinCoord = []
    
    def __repr__(self):
        return "SingleTimeProfile({})".format(self.Coord)
    
    def __str__(self):
        return "Time Profile Object at coord: {}".format(self.Coord)
    
    def AddTwin(self,TwinCoord):
        self.TwinCoord = TwinCoord

    def appendFieldValues(self, time, dispx, dispy, velx, vely):
        self.Time.append(time)
        self.DispX.append(dispx)
        self.DispY.append(dispy)
        self.VelX.append(velx)
        self.VelY.append(vely)

    def AssimilateSingleTimeProfile(self, AnotherSTP_obj):
        self.Time.extend( AnotherSTP_obj.Time)
        self.DispX.extend(AnotherSTP_obj.DispX)
        self.DispY.extend(AnotherSTP_obj.DispY)
        self.VelX.extend( AnotherSTP_obj.VelX)
        self.VelY.extend( AnotherSTP_obj.VelY)

    
    def PrintValues(self):
        print("Coordinate:")
        print(self.Coord)
        print("Displacement:")
        print(self.DispX)
        print(self.DispY)
        print("Velocity:")
        print(self.VelX)
        print(self.VelY)


