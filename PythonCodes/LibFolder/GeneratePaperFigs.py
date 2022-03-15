import numpy as np
import matplotlib.pyplot as plt

import os, sys, time
from scipy.interpolate import RectBivariateSpline
from sklearn.metrics.pairwise import euclidean_distances

from matplotlib.ticker import FuncFormatter, MaxNLocator
import matplotlib.lines as mlines


from se2waveload import *
from Lib_GeneralFunctions import *
from Lib_GeneralSignalProcNAnalysis import *
from Lib_SigmoidProcessing import *
import pandas as pd

from mpl_toolkits.axes_grid1 import ImageGrid
from mpl_toolkits.axes_grid1.inset_locator import (inset_axes, InsetPosition,mark_inset)
import itertools


import string

def LabelizeAxisList(AxisList,Pos=[-0.1, 1.1], OffsetLabel=0,**kwargs):
    for n, ax in enumerate(AxisList):  
        ax.text(Pos[0], Pos[1], string.ascii_uppercase[n+OffsetLabel], transform=ax.transAxes, **kwargs)



# Sigmoid or any function of interest to represent the center of the fault / Zero level set function
def func(x, k=-0.0002, amp = 2.0):
    fx = amp * (x - x * k) / (k - abs(x) * 2.0 * k + 1.0)
    return fx

# The respective derivative ofthe previous zero level set function
def func_der(x, k=-0.0002, amp = 2.0):
    fx_prime = amp * (1 - k * k) / ((k - abs(x) * 2.0 * k + 1.0)*(k - abs(x) * 2.0 * k + 1.0))
    return fx_prime


# Sigmoid or any function of interest to represent the center of the fault / Zero level set function
def Tiltfunc(x, theta = 45*np.pi/180):
    fx = x*np.tan(theta)
    return fx

def Tiltfunc_der(x, theta = 45*np.pi/180):
    fx_prime = x.copy()
    fx_prime.fill(np.tan(theta))
    return fx_prime

class ZeroLevelSet:
    def __init__(self, Xval, Fxval, FxPrimeVal, GeometryDescription):
        self.Xval = Xval
        self.Fxval = Fxval
        self.FxPrimeVal = FxPrimeVal
        
        self.GeometryDescription = GeometryDescription
        
        self.Normal = np.array(self.NormalVector(self.FxPrimeVal))
        self.Tangent = np.array(self.TangentVector(self.FxPrimeVal))
        
    def __repr__(self):
        return "Zero level set: {GeometryDescription} geometry".format(GeometryDescription=self.GeometryDescription)

    def __str__(self):
        return "Zero level set: {GeometryDescription} geometry".format(GeometryDescription=self.GeometryDescription)
    
    
    def PlotZeroLevelSet(self):
        plt.plot(self.Xval,self.Fxval,"k-")
    
    # Tangent vector for a given derivative
    def TangentVector(self, fPrimeX):
        mag = np.sqrt(1.0 + fPrimeX * fPrimeX)

        TangentX = 1.0/mag
        TangentY = fPrimeX/mag
        return TangentX, TangentY

    # Normal vector for a given derivative
    def NormalVector(self,fPrimeX):
        mag = np.sqrt(1.0 + fPrimeX * fPrimeX)

        NormalX = -fPrimeX/mag
        NormalY = 1.0/mag

        return NormalX, NormalY
    
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

def FormatAx(ax,MaxNLoc=5,Axins=False):
    ax.set_aspect("equal")
    ax.xaxis.set_major_locator(MaxNLocator(MaxNLoc))
    ax.yaxis.set_major_locator(MaxNLocator(MaxNLoc))
    ax.yaxis.major.formatter.set_powerlimits((0,0)) 
    ax.xaxis.major.formatter.set_powerlimits((0,0)) 
    
    if Axins:
        ax.get_xaxis().get_offset_text().set_visible(False)
        ax_max = max(ax.get_xticks())
        exponent_axis = np.floor(np.log10(ax_max)).astype(int)
        ax.annotate(r'$\times$10$^{%i}$'%(exponent_axis),
                    xy=(0.5, 0.1), xycoords='axes fraction')

def FormatAxNormal(ax):
    ax.set_aspect("equal")
    ax.xaxis.set_major_locator(MaxNLocator(5))
    ax.yaxis.set_major_locator(MaxNLocator(5))

def PlotDomain(CoorX, CoorY, Field, FieldName,TimeTxt,**kwargs):
    try:
      fig = plt.figure(figsize = (10, 10),dpi=300, constrained_layout=True)
      gs = fig.add_gridspec(1, 1)
      ax = fig.add_subplot(gs[:, :])
    except:
      fig = plt.figure(figsize = (10, 10),dpi=300)
      ax = fig.add_subplot(1,1,1)
    ax.set_title("{FName}".format(FName = FieldName[0]))
    ax.set_xlabel("X-Coordinate [m]"), ax.set_ylabel("Y-Coordinate [m]")
    ax.set_aspect('equal', 'box')
    img = ax.pcolormesh(CoorX, CoorY, Field,**kwargs)

    ax.annotate(text="time [s]: {0:.2f}".format(TimeTxt),xy=[0.8,0.1], xycoords= "axes fraction")
    cbar = fig.colorbar(img, shrink=.5)
    cbar.ax.set_ylabel(FieldName[1])
    
    return fig, img,ax


def PlotFullSetup(CoorX, CoorY, Field1, Field2, StressFromPickle, FieldNames,TimeTxt,InsetZoom=[6250,6750,3400,3900],**kwargs):    
    fig = plt.figure(figsize = (12, 8),dpi=300) #constrained_layout=True
    gs = fig.add_gridspec(2, 3, wspace=0.15,hspace=0.2)

    
    ax01 = fig.add_subplot(gs[0, 0])
    ax02 = fig.add_subplot(gs[0, 1])
    ax03 = fig.add_subplot(gs[0, 2])
    ax1 = fig.add_subplot(gs[-1, 0])
    ax2 = fig.add_subplot(gs[-1, 1])
    ax3 = fig.add_subplot(gs[-1, 2])
       
    
    #Plot
    #ax1.set_title("{FName}".format(FName = FieldNames[0]))
    ax1.set_xlabel("$x$ [m]")
    ax2.set_xlabel("$x$ [m]")
    ax3.set_xlabel("$x$ [m]")
    ax1.set_ylabel("$y$ [m]")
    
    
    FormatAx(ax1)
    FormatAx(ax2)
    FormatAx(ax3)
        
    img1 = ax1.pcolormesh(CoorX, CoorY, Field1,**kwargs)
    img2 = ax2.pcolormesh(CoorX, CoorY, Field2,**kwargs)
    img3 = ax3.pcolormesh(StressFromPickle[0], StressFromPickle[1], StressFromPickle[2], shading="flat",
                          vmax = 2e7, vmin= -2e7, **kwargs)
    
    #ax2.tick_params(labelleft=False)
    #ax3.tick_params(labelleft=False)
    
    #ax2.yaxis.get_major_formatter().set_scientific(False)
    #ax3.yaxis.get_major_formatter().set_scientific(False)
    ax1.annotate(text="time [s]: {0:.2f}".format(TimeTxt),xy=[0.05,0.9], xycoords= "axes fraction")
    
    # Colorbar for the ax1
    cbaxes = inset_axes(ax1,width="40%",height="4%",loc=3, borderpad=2)
    plt.colorbar(img1,cax=cbaxes,orientation="horizontal", label=r"$u_{x}$ [m]")
    cbaxes.xaxis.set_label_position('top')
    
    # Colorbar for the ax2
    cbaxes = inset_axes(ax2,width="40%",height="4%",loc=3, borderpad=2)
    plt.colorbar(img2,cax=cbaxes,orientation="horizontal", label=r"$v_{x}$ [m/s]")
    cbaxes.xaxis.set_label_position('top')
    
    # Colorbar for the ax3
    cbaxes = inset_axes(ax3,width="40%",height="4%",loc=3, borderpad=2)
    plt.colorbar(img3,cax=cbaxes,orientation="horizontal", label=r"$\sigma_{12}$ [Pa]")
    cbaxes.xaxis.set_label_position('top')

    # Give the number of ticks in the colorbar
    cbaxes.xaxis.set_major_locator(MaxNLocator(4))

    # Give an offset for the scientific notation exponent
    cbaxes.get_xaxis().get_offset_text().set_visible(False)
    ax_max = max(cbaxes.get_xticks())
    exponent_axis = np.floor(np.log10(ax_max)).astype(int)
    cbaxes.annotate(r'$\times$10$^{%i}$'%(exponent_axis),
                xy=(1.01, -.01), xycoords='axes fraction')
    
    InsetLoc = [0.67, 0.08, 0.3, 0.3]    
    # Inset plot for the ax2
    axins = ax2.inset_axes(InsetLoc)
    axins.pcolormesh(CoorX, CoorY, Field2, edgecolors='silver',lw='0.1', **kwargs)
    
    axins.set_xlim(InsetZoom[0], InsetZoom[1])
    axins.set_ylim(InsetZoom[2], InsetZoom[3])
    #axins.set_xticklabels('')
    #axins.set_yticklabels('')
    FormatAx(axins,MaxNLoc=2,Axins=True)
    
    mark_inset(ax2, axins,loc1=2, loc2=1, edgecolor="black",ec=".5",linewidth=.5)


    # Inset plot for the ax3
    axins = ax3.inset_axes(InsetLoc)
    axins.pcolormesh(StressFromPickle[0], StressFromPickle[1], StressFromPickle[2], shading="flat", edgecolors='silver',lw='0.1',
                     vmax = 2e7, vmin= -2e7, **kwargs)
    
    axins.set_xlim(InsetZoom[0], InsetZoom[1])
    axins.set_ylim(InsetZoom[2], InsetZoom[3])
    #axins.set_xticklabels('')
    #axins.set_yticklabels('')
    FormatAx(axins,MaxNLoc=2,Axins=True)
    mark_inset(ax3, axins,loc1=2, loc2=1, edgecolor="black",ec=".5",linewidth=.5)

    gs.tight_layout(fig)
    gs.update(top=0.95)
    
    ax = [ax01,ax02,ax03,ax1,ax2,ax3]
    #cbar.ax.set_ylabel(FieldName[1])
    return fig, ax

def PlotHalfSetup(CoorX, CoorY, Field1, Field2, StressFromPickle, FieldNames,TimeTxt,InsetZoom=[6250,6750,3400,3900],**kwargs):    
    fig = plt.figure(figsize = (12, 4),dpi=300) #constrained_layout=True
    gs = fig.add_gridspec(1, 3, wspace=0.15,hspace=0.2)


    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[0, 2])
       
    
    #Plot
    #ax1.set_title("{FName}".format(FName = FieldNames[0]))
    ax2.set_xlabel("$x$ [m]")
    ax1.set_ylabel("$y$ [m]")
    
    FormatAx(ax1)
    FormatAx(ax2)
    FormatAx(ax3)
        
    img1 = ax1.pcolormesh(CoorX, CoorY, Field1,**kwargs)
    img2 = ax2.pcolormesh(CoorX, CoorY, Field2,**kwargs)
    img3 = ax3.pcolormesh(StressFromPickle[0], StressFromPickle[1], StressFromPickle[2], shading="flat",
                          vmax = 2e7, vmin= -2e7, **kwargs)
    
    ax2.tick_params(labelleft=False)
    ax3.tick_params(labelleft=False)
    
    ax2.yaxis.get_major_formatter().set_scientific(False)
    ax3.yaxis.get_major_formatter().set_scientific(False)
    ax1.annotate(text="t: {0:.2f} s".format(TimeTxt),xy=[0.05,0.9], xycoords= "axes fraction")
    
    # Colorbar for the ax1
    cbaxes = inset_axes(ax1,width="40%",height="4%",loc=3, borderpad=2)
    plt.colorbar(img1,cax=cbaxes,orientation="horizontal", label=r"$u_{x}$ [m]")
    cbaxes.xaxis.set_label_position('top')
    
    # Colorbar for the ax2
    cbaxes = inset_axes(ax2,width="40%",height="4%",loc=3, borderpad=2)
    plt.colorbar(img2,cax=cbaxes,orientation="horizontal", label=r"$v_{x}$ [m/s]")
    cbaxes.xaxis.set_label_position('top')
    
    # Colorbar for the ax3
    cbaxes = inset_axes(ax3,width="40%",height="4%",loc=3, borderpad=2)
    plt.colorbar(img3,cax=cbaxes,orientation="horizontal", label=r"$\sigma_{12}$ [Pa]")
    cbaxes.xaxis.set_label_position('top')

    # Give the number of ticks in the colorbar
    cbaxes.xaxis.set_major_locator(MaxNLocator(4))

    # Give an offset for the scientific notation exponent
    cbaxes.get_xaxis().get_offset_text().set_visible(False)
    ax_max = max(cbaxes.get_xticks())
    exponent_axis = np.floor(np.log10(ax_max)).astype(int)
    cbaxes.annotate(r'$\times$10$^{%i}$'%(exponent_axis),
                xy=(1.01, -.01), xycoords='axes fraction')
    
    
    InsetLoc = [0.67, 0.08, 0.3, 0.3]
    # Inset plot for the ax2
    axins = ax2.inset_axes(InsetLoc)
    axins.pcolormesh(CoorX, CoorY, Field2, edgecolors='silver',lw='0.1', **kwargs)
    
    axins.set_xlim(InsetZoom[0], InsetZoom[1])
    axins.set_ylim(InsetZoom[2], InsetZoom[3])
    #axins.set_xticklabels('')
    #axins.set_yticklabels('')
    
    mark_inset(ax2, axins,loc1=2, loc2=1, edgecolor="black",ec=".5",linewidth=.5)


    # Inset plot for the ax3
    axins = ax3.inset_axes(InsetLoc)
    axins.pcolormesh(StressFromPickle[0], StressFromPickle[1], StressFromPickle[2], shading="flat", edgecolors='silver',lw='0.1',
                     vmax = 2e7, vmin= -2e7, **kwargs)
    
    axins.set_xlim(InsetZoom[0], InsetZoom[1])
    axins.set_ylim(InsetZoom[2], InsetZoom[3])
    axins.set_xticklabels('')
    axins.set_yticklabels('')
    mark_inset(ax3, axins,loc1=2, loc2=1, edgecolor="black",ec=".5",linewidth=.5)

    #gs.tight_layout(fig)
    #gs.update(top=0.95)
    
    ax = [ax1,ax2,ax3]
    #cbar.ax.set_ylabel(FieldName[1])
    return fig, ax

def Plot4KomaSetup(CoorX, CoorY, Field1, Field2, FieldNames,TimeTxt,InsetZoom=[6250,6750,3400,3900],**kwargs):    
    fig = plt.figure(figsize = (8, 8),dpi=300) #constrained_layout=True
    gs = fig.add_gridspec(2, 2, wspace=0.15,hspace=0.2)

    
    ax01 = fig.add_subplot(gs[0, 0])
    ax02 = fig.add_subplot(gs[0, 1])
    ax1 = fig.add_subplot(gs[-1, 0])
    ax2 = fig.add_subplot(gs[-1, 1])
       
    ax = [ax01,ax02,ax1,ax2]
    #Plot
    #ax1.set_title("{FName}".format(FName = FieldNames[0]))
    ax2.set_xlabel("$x$ [m]")
    ax1.set_xlabel("$x$ [m]")
    ax1.set_ylabel("$y$ [m]")
    
    FormatAx(ax1)
    FormatAx(ax2)
        
    img1 = ax1.pcolormesh(CoorX, CoorY, Field1,**kwargs)
    img2 = ax2.pcolormesh(CoorX, CoorY, Field2,**kwargs)
    
    #ax2.tick_params(labelleft=False,labelright=True)
    
    #ax2.yaxis.get_major_formatter().set_scientific(False)
    ax1.annotate(text="time [s]: {0:.2f}".format(TimeTxt),xy=[0.05,0.9], xycoords= "axes fraction")
    
    cbaxes = inset_axes(ax1, width="40%",height="4%",loc=3, borderpad=2)
    plt.colorbar(img1,cax=cbaxes,orientation="horizontal", label=r"$u_{x}$ [m]")
    cbaxes.xaxis.set_label_position('top')
    
    cbaxes = inset_axes(ax2, width="40%",height="4%",loc=3, borderpad=2)
    plt.colorbar(img2,cax=cbaxes,orientation="horizontal", label=r"$v_{x}$ [m/s]")
    cbaxes.xaxis.set_label_position('top')
    
      
    axins = ax2.inset_axes([0.67, 0.1, 0.3, 0.3])
    axins.pcolormesh(CoorX, CoorY, Field2, edgecolors='silver',lw='0.1', **kwargs)
    
    axins.set_xlim(InsetZoom[0], InsetZoom[1])
    axins.set_ylim(InsetZoom[2], InsetZoom[3])
    FormatAx(axins,MaxNLoc=2,Axins=True)

    
    #axins.grid(True, which='both', axis='both', linestyle='-', color='k')
    mark_inset(ax2, axins,loc1=2, loc2=1, edgecolor="black",ec="0.5",linewidth=.5)

    gs.tight_layout(fig)
    #gs.update(top=0.95)
    
    #cbar.ax.set_ylabel(FieldName[1])
    return fig, ax

def PlotF4Setup(CoorX, CoorY, Field1, StressFromPickle, FieldNames,TimeTxt,InsetZoom=[6250,6750,3400,3900],**kwargs):    
    fig = plt.figure(figsize = (8, 8),dpi=300) #constrained_layout=True
    gs = fig.add_gridspec(2, 2, wspace=0.15,hspace=0.2)

    
    ax01 = fig.add_subplot(gs[0, 0])
    ax02 = fig.add_subplot(gs[0, 1])
    ax1 = fig.add_subplot(gs[-1, 0])
    ax2 = fig.add_subplot(gs[-1, 1])
       
    ax = [ax01,ax02,ax1,ax2]
    #Plot
    #ax1.set_title("{FName}".format(FName = FieldNames[0]))
    ax2.set_xlabel("$x$ [m]")
    ax1.set_xlabel("$x$ [m]")
    ax1.set_ylabel("$y$ [m]")
    
    FormatAx(ax1)
    FormatAx(ax2)
        
    img1 = ax1.pcolormesh(CoorX, CoorY, Field1,**kwargs)
    img2 = ax2.pcolormesh(StressFromPickle[0], StressFromPickle[1], StressFromPickle[2], shading="flat",
                          vmax = 2e7, vmin= -2e7, **kwargs)
    
    #ax2.tick_params(labelleft=True)
    
    #ax2.yaxis.get_major_formatter().set_scientific(False)
    ax1.annotate(text="time [s]: {0:.2f}".format(TimeTxt),xy=[0.05,0.9], xycoords= "axes fraction")
    
    # Colorbar for the ax1
    cbaxes = inset_axes(ax1, width="35%",height="4%",loc=3, borderpad=2)
    plt.colorbar(img1,cax=cbaxes,orientation="horizontal", label=r"$v_{x}$ [m/s]")
    cbaxes.xaxis.set_label_position('top')
    
    # Colorbar for the ax2
    cbaxes = inset_axes(ax2, width="35%",height="4%",loc=3, borderpad=2)
    plt.colorbar(img2,cax=cbaxes,orientation="horizontal", label=r"$\sigma_{12}$ [Pa]")
    cbaxes.xaxis.set_label_position('top')

    # Give the number of ticks in the colorbar
    cbaxes.xaxis.set_major_locator(MaxNLocator(4))

    # Give an offset for the scientific notation exponent
    cbaxes.get_xaxis().get_offset_text().set_visible(False)
    ax_max = max(cbaxes.get_xticks())
    exponent_axis = np.floor(np.log10(ax_max)).astype(int)
    cbaxes.annotate(r'$\times$10$^{%i}$'%(exponent_axis),
                xy=(1.01, -.01), xycoords='axes fraction')
    
    # Inset plot for the ax1
    axins = ax1.inset_axes([0.67, 0.1, 0.3, 0.3])
    axins.pcolormesh(CoorX, CoorY, Field1, edgecolors='silver',lw='0.1', **kwargs)
    
    axins.set_xlim(InsetZoom[0], InsetZoom[1])
    axins.set_ylim(InsetZoom[2], InsetZoom[3])
    #axins.set_xticklabels('')
    #axins.set_yticklabels('')
    FormatAx(axins,MaxNLoc=2,Axins=True)
    
    mark_inset(ax1, axins,loc1=2, loc2=1, edgecolor="black",ec="0.5",linewidth=.5)
    
    # Inset plot for the ax2
    axins2 = ax2.inset_axes([0.67, 0.1, 0.3, 0.3])
    axins2.pcolormesh(StressFromPickle[0], StressFromPickle[1], StressFromPickle[2], shading="flat", edgecolors='silver',lw='0.1',
                      vmax = 2e7, vmin= -2e7, **kwargs)
    
    axins2.set_xlim(InsetZoom[0], InsetZoom[1])
    axins2.set_ylim(InsetZoom[2], InsetZoom[3])
    #axins2.set_xticklabels('')
    #axins2.set_yticklabels('')
    FormatAx(axins2,MaxNLoc=2,Axins=True)
    
    mark_inset(ax2, axins2,loc1=2, loc2=1, edgecolor="black",ec="0.5",linewidth=.5)

    gs.tight_layout(fig)
    #gs.update(top=0.95)
    
    return fig, ax

def PlotVistas(CoorX, CoorY, Field1, Field2, Field3, Field4, TimeTxt, InsetZoom=[6250,6750,3400,3900],**kwargs):    
    fig = plt.figure(figsize = (8, 8),dpi=300) #constrained_layout=True
    gs = fig.add_gridspec(2, 2, wspace=0.15,hspace=0.2)

    
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[-1, 0])
    ax4 = fig.add_subplot(gs[-1, 1])
       
    ax = [ax1,ax2,ax3,ax4]

    #Plot
    ax4.set_xlabel("$x$ [m]")
    ax3.set_xlabel("$x$ [m]")
    ax1.set_ylabel("$y$ [m]")
    ax3.set_ylabel("$y$ [m]")
    
    FormatAx(ax1)
    FormatAx(ax2)
        
    img1 = ax1.pcolormesh(CoorX, CoorY, Field1,**kwargs)
    img2 = ax2.pcolormesh(CoorX, CoorY, Field2,**kwargs)
    img3 = ax3.pcolormesh(CoorX, CoorY, Field3,**kwargs)
    img4 = ax4.pcolormesh(CoorX, CoorY, Field4,**kwargs)
    
    ax2.tick_params(labelleft=False)
    
    ax2.yaxis.get_major_formatter().set_scientific(False)
    ax1.annotate(text="time [s]: {0:.2f}".format(TimeTxt),xy=[0.05,0.9], xycoords= "axes fraction")
    
    cbaxes = inset_axes(ax1, width="40%",height="4%",loc=3, borderpad=2)
    plt.colorbar(img1,cax=cbaxes,orientation="horizontal", label=r"$u_{x}$ [m]")
    cbaxes.xaxis.set_label_position('top')
    
    cbaxes = inset_axes(ax2, width="40%",height="4%",loc=3, borderpad=2)
    plt.colorbar(img2,cax=cbaxes,orientation="horizontal", label=r"$v_{x}$ [m/s]")
    cbaxes.xaxis.set_label_position('top')

    cbaxes = inset_axes(ax3, width="40%",height="4%",loc=3, borderpad=2)
    plt.colorbar(img3,cax=cbaxes,orientation="horizontal", label=r"$u_{y}$ [m]")
    cbaxes.xaxis.set_label_position('top')
    
    cbaxes = inset_axes(ax4, width="40%",height="4%",loc=3, borderpad=2)
    plt.colorbar(img4,cax=cbaxes,orientation="horizontal", label=r"$v_{y}$ [m/s]")
    cbaxes.xaxis.set_label_position('top')
    
      
    axins = ax2.inset_axes([0.65, 0.05, 0.3, 0.3])
    axins.pcolormesh(CoorX, CoorY, Field2, edgecolors='silver',lw='0.1', **kwargs)
    
    axins.set_xlim(InsetZoom[0], InsetZoom[1])
    axins.set_ylim(InsetZoom[2], InsetZoom[3])
    axins.set_xticklabels('')
    axins.set_yticklabels('')
    
    #axins.grid(True, which='both', axis='both', linestyle='-', color='k')
    mark_inset(ax2, axins,loc1=2, loc2=1, edgecolor="black",ec="0.5",linewidth=.5)

    gs.tight_layout(fig)
    gs.update(top=0.95)
    
    #cbar.ax.set_ylabel(FieldName[1])
    return fig, ax

# Save into a class the 
class SSCreference:
    def __init__(self, filename, coordinates, RefSource="SEM2DPACK"):
        
        line = pd.read_csv(filename.format("slip"), header=None)
        self.Time = line[0]
        self.Slip = line[1]
        
        line = pd.read_csv(filename.format("sr"), header=None)
        self.SlipRate = line[1]
        
        self.Coord = coordinates #Only used for labels and printing
        self.RefSource = RefSource
    #end __init__
    
    # Default object printing information
    def __repr__(self):
        return "The TPV3reference object was generated from: {} and the receiver is located at {}".format(self.RefSource, self.Coord)
    #end __repr__
    
    def __str__(self):
        return "The TPV3reference object was generated from: {} and the receiver is located at {}".format(self.RefSource, self.Coord)
    #end __str__
    
    def PlotReference(self, ax, SlipSlipRate, filtering=True, **kwargs):
        
        if SlipSlipRate=="Slip":
            if(filtering):
                ax.plot(self.Time, Butterworth(self.Slip, **kwargs), label = "", c = "k", ls = "--", zorder=1)
            else:
                ax.plot(self.Time, self.Slip, label = "", c = "k", ls = "--", zorder=1)
        elif SlipSlipRate=="SlipRate":
            if(filtering):
                ax.plot(self.Time, Butterworth(self.SlipRate, **kwargs), label = "", c = "k", ls = "--", zorder=1)
            else:
                ax.plot(self.Time, self.SlipRate, label = "", c = "k", ls = "--", zorder=1)
            
        return ax
# Save into a class the 


class TPV3reference:
    def __init__(self, filename, coordinates, RefSource="SEM2DPACK"):
        
        line = pd.read_csv(filename.format("slip"), header=None)
        self.Time = line[0]
        self.Slip = line[1]
        
        line = pd.read_csv(filename.format("sr"), header=None)
        self.SlipRate = line[1]
        
        self.Coord = coordinates #Only used for labels and 
        self.RefSource = RefSource
    #end __init__
    
    # Default object printing information
    def __repr__(self):
        return "The TPV3 reference object was generated from: {} and the receiver is located at {}".format(self.RefSource, self.Coord)
    #end __repr__
    
    def __str__(self):
        return "The TPV3 reference object was generated from: {} and the receiver is located at {}".format(self.RefSource, self.Coord)
    #end __str__
    
    def PlotReference(self, ax, SlipSlipRate, filtering=True, **kwargs):
        
        if SlipSlipRate=="Slip":
            if(filtering):
                ax.plot(self.Time, Butterworth(self.Slip, **kwargs), label = "", c = "k", ls = "--", zorder=1)
            else:
                ax.plot(self.Time, self.Slip, label = "", c = "k", ls = "--", zorder=1)
        elif SlipSlipRate=="SlipRate":
            if(filtering):
                ax.plot(self.Time, Butterworth(self.SlipRate, **kwargs), label = "", c = "k", ls = "--", zorder=1)
            else:
                ax.plot(self.Time, self.SlipRate, label = "", c = "k", ls = "--", zorder=1)
            
        return ax
    
def GenericFigAxis():
    fig = plt.figure(figsize=[15,5])
    gs = GridSpec(1, 2)

    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    
    
    return fig, [ax1, ax2]
    

def formatGivenAxes(AxesList,inverted=True):
    for i, ax in enumerate(AxesList):
        ax.set_xlim(-.2,4)
        ax.set_ylim(-.5,8.25)
        
        
        #ax.xaxis.set_label_position('top') 
        ax.xaxis.set_ticks_position('both')
        ax.xaxis.set_major_locator(MaxNLocator(5))
        
        ax.set_xlabel("time [s]")
                
    Lines = AxesList[-1].get_lines()
    
    ReceiversLabelList = ['0 km','2 km','4 km', '6 km', '8 km']
    if (inverted):
        ReceiversLabelList.reverse()
      
    
    legend2 = AxesList[0].legend(Lines, ReceiversLabelList , loc=2)
    AxesList[0].add_artist(legend2)
    
    
    LinesContDisc = []
    LinesContDisc.append(mlines.Line2D([], [], color = "k", ls = "--",
                         linewidth = 1, label ="SEM2DPACK" ))
    LinesContDisc.append(mlines.Line2D([], [], color = "k", ls = "-",
                         linewidth = 1, label ="se2dr" ))
    
    legendContDisc = AxesList[1].legend(LinesContDisc, ["SEM2DPACK","se2dr"], loc = 2)
    AxesList[1].add_artist(legendContDisc)
    
    AxesList[1].set_ylabel("Slip rate [m/s]")
    AxesList[0].set_ylabel("Slip [m]")
    

def flip(items, ncol):
    return itertools.chain(*[items[i::ncol] for i in range(ncol)])

def Format_LabelsOnFig_formatAxis(fig, AxesList,inverted=False, ncols=5,ncolSim=1, **kwargs):
    for i, ax in enumerate(AxesList):
        ax.set_xlim(-.2,4)
        ax.set_ylim(-.5,8.25)
        
        
        #ax.xaxis.set_label_position('top') 
        ax.xaxis.set_ticks_position('both')
        ax.xaxis.set_major_locator(MaxNLocator(5))
        
        ax.set_xlabel("time [s]")
        
    AxesList[1].set_ylabel("Slip rate [m/s]")
    AxesList[0].set_ylabel("Slip [m]")
                
    Lines = AxesList[-1].get_lines()[:5]
    
    ReceiversLabelList = ['0 km','2 km','4 km', '6 km', '8 km']
    if (inverted):
        Lines.reverse()
    
    
    AxReference = kwargs.pop("AxLabelLocRef", AxesList[1])
    
    Lines = list(flip(Lines, ncols))
    ReceiverLabelLines = list(flip(ReceiversLabelList, ncols))
    ReceiverArtist = AxReference.legend(Lines,     # The line objects
                                ReceiverLabelLines,   # The labels for each line
                                loc='upper right',   # Position of legend
                                borderaxespad=0.1,    # Small spacing around legend box
                                ncol = ncols,  bbox_to_anchor=(1.0, 1.1),**kwargs
                                )

    fig.add_artist(ReceiverArtist)
        
    
    LinesContDisc = []
    LinesContDisc.append(mlines.Line2D([], [], color = "k", ls = "--",
                         linewidth = 1, label ="SEM2DPACK" ))
    LinesContDisc.append(mlines.Line2D([], [], color = "k", ls = "-",
                         linewidth = 1, label ="se2dr" ))
    
    SimulationArtist = AxesList[0].legend(LinesContDisc,     # The line objects
                                  ["SEM2DPACK","se2dr"],   # The labels for each line
                                  loc='upper left',   # Position of legend
                                  borderaxespad=0.1,    # Small spacing around legend box
                                  ncol = ncolSim, bbox_to_anchor=(0.0, 1.1)
                                  )

    fig.add_artist(SimulationArtist)
    
    


def Multi_format_axes(fig,cmap, LabelsPerColor):
    """
    Format a figure that contains different files with 
    information from several receivers for simulations under sets of blending parameters.
    """
    ColorDict = dict(enumerate(LabelsPerColor)) 
    
    
    for i, ax in enumerate(fig.axes):
        ax.set_xlim(-0.5,4)
        ax.set_ylim(-0.5,8)
        ax.set_xlabel("time(s)")
    Lines = []
    for idx,colcol in enumerate(cmap.colors):
        Lines.append(mlines.Line2D([], [], color = colcol,
                     linewidth = 3, label = ColorDict.get(-idx)))
    
    legend2 = fig.axes[-1].legend(Lines, LabelsPerColor, loc = 2)
    fig.axes[-1].add_artist(legend2)
    fig.axes[-1].set_ylabel("Slip rate (m/s)")
    fig.axes[0].set_ylabel("Slip (m)")
