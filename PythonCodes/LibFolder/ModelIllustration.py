import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import math
    

def DrawPlotAxes(BGLoc, BGScale, axesDirList, ax, Label=""):
    
    if Label=="":
        Labels=[""]*(len(axesDir))
    else :
        Labels = Label
    
    for idx,Dir in enumerate(axesDirList):
        ax.annotate("", 
                    xytext=(BGLoc[0]-0.5*(Dir[0])*BGScale,BGLoc[1]-0.5*(Dir[1])*BGScale), 
                    xy=(BGLoc[0]+0.5*(Dir[0])*BGScale,BGLoc[1]+0.5*(Dir[1])*BGScale),
                    arrowprops=dict(arrowstyle="->"),horizontalalignment='left')
        ax.annotate(Labels[idx], 
                    xytext=(BGLoc[0]-0.7*(Dir[0])*BGScale,BGLoc[1]-0.7*(Dir[1])*BGScale),
                    xy=(BGLoc[0],BGLoc[1]),
                    horizontalalignment='center',
                    verticalalignment = 'center'
                   )


def LocateSlipReceivers_MeshAligned(ax,loc,Color):
    OffsetY = 2*loc[1]

    MarkerSpecs = {"marker":".","facecolors":Color,"edgecolors":"k","s":150,"zorder":9}
    ax.scatter(loc[0],loc[1], **MarkerSpecs)
    ax.scatter(loc[0],loc[1]-OffsetY, **MarkerSpecs)
        
def GenKostrovCase(ax):
    Delta = 0.04
    
    ####### Add Grid and axis specs
    ax.set_xlim(-.05,.85)
    ax.set_ylim(-.45,.45)
    ax.xaxis.set_major_locator(plt.MultipleLocator(1.0))
    ax.xaxis.set_minor_locator(plt.MultipleLocator(.02))
    ax.yaxis.set_major_locator(plt.MultipleLocator(1.0))
    ax.yaxis.set_minor_locator(plt.MultipleLocator(.02))
    ax.grid(which='major', axis='x', linewidth=0.75, linestyle='-', color='0.75')
    ax.grid(which='minor', axis='x', linewidth=0.25, linestyle='-', color='0.75')
    ax.grid(which='major', axis='y', linewidth=0.75, linestyle='-', color='0.75')
    ax.grid(which='minor', axis='y', linewidth=0.25, linestyle='-', color='0.75')
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_aspect('equal')
    #####################################


    ####### Add patches

    ax.add_patch(
                 patches.Rectangle((-1, -Delta), 2., 2*Delta, fill=True, color="lightgray" ) 
                ) 
    ax.add_patch(
                 patches.Rectangle((-1, -Delta), 2., 2*Delta, fill=False ) 
                ) 
    #####################################


    ####### Add axes of shear and normal background stress

    BGLoc = [0.1,0.25]
    BGScale = 0.15
    axesDir = [[0,-1],[1,0]]
    Label = ['$-\sigma^b_{22}$','$\sigma^b_{12}$']
    DrawPlotAxes(BGLoc, BGScale, axesDir,ax, Label)
    #####################################

    
    ###### Outside Arrows
    OffsetArr = 0.02
    
    ax.annotate('', xy=(-OffsetArr, 0), xycoords='axes fraction', xytext=(-OffsetArr,1), 
            arrowprops=dict(arrowstyle="<->", color='k'))
    
    ax.annotate('', xy=(0, -OffsetArr), xycoords='axes fraction', xytext=(1, -OffsetArr), 
            arrowprops=dict(arrowstyle="<->", color='k'))
    

    ########################################

    ####### 
    BGLoc = [.5, 0]
    BGScale = 2*Delta
    Dir = [0,-1]
    Label = ['2$\delta$']
    ax.annotate("", 
                xytext=(BGLoc[0]-0.5*(Dir[0])*BGScale,BGLoc[1]-0.5*(Dir[1])*BGScale), 
                xy=(BGLoc[0]+0.5*(Dir[0])*BGScale,BGLoc[1]+0.5*(Dir[1])*BGScale),
                arrowprops=dict(arrowstyle="<->"),horizontalalignment='left')
    ax.annotate(Label[0], 
                xytext=(BGLoc[0]+.02,BGLoc[1]), 
                xy=(BGLoc[0],BGLoc[1]),
                horizontalalignment='left',verticalalignment = 'center')
    #####################################

def GenTPV3Case(ax):
    ####### Add Grid and axis specs
    ax.set_xlim(-.2,.6)
    ax.set_ylim(-.4,.4)
    ax.xaxis.set_major_locator(plt.MultipleLocator(1.0))
    ax.xaxis.set_minor_locator(plt.MultipleLocator(.02))
    ax.yaxis.set_major_locator(plt.MultipleLocator(1.0))
    ax.yaxis.set_minor_locator(plt.MultipleLocator(.02))
    ax.grid(which='major', axis='x', linewidth=0.75, linestyle='-', color='0.75')
    ax.grid(which='minor', axis='x', linewidth=0.25, linestyle='-', color='0.75')
    ax.grid(which='major', axis='y', linewidth=0.75, linestyle='-', color='0.75')
    ax.grid(which='minor', axis='y', linewidth=0.25, linestyle='-', color='0.75')
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_aspect('equal')
    #####################################


    ####### Add maximum stress direction

    ax.add_patch(
                 patches.Rectangle((-1, -0.04), 2., 0.08, fill=True, color="lightgray" ) 
                ) 
    ax.add_patch(
                 patches.Rectangle((-1, -0.04), 2., 0.08, fill=False ) 
                ) 

    ax.add_patch(
                 patches.Rectangle((-.15, -0.04), .3, 0.08, fill=True, color="gray" ) 
                ) 
    ax.add_patch(
                 patches.Rectangle((-.15, -0.04), .3, 0.08, fill=False ) 
                ) 


    #####################################

    ####### Add axes of shear and normal background stress
    ####### 
    BGLoc = [0.0,0.25]
    BGScale = 0.15
    axesDir = [[0,-1],[1,0]]
    Label = ['$-\sigma^b_{22}$','$\sigma^b_{12}$']
    DrawPlotAxes(BGLoc, BGScale, axesDir,ax,Label)
    #####################################

    ####### Thickness Indicator
    BGLoc = [.35, 0]
    BGScale = 0.08
    Dir = [0,-1]
    Label = ['2$\delta$']
    ax.annotate("", 
                xytext=(BGLoc[0]-0.5*(Dir[0])*BGScale,BGLoc[1]-0.5*(Dir[1])*BGScale), 
                xy=(BGLoc[0]+0.5*(Dir[0])*BGScale,BGLoc[1]+0.5*(Dir[1])*BGScale),
                arrowprops=dict(arrowstyle="<->"),horizontalalignment='left')
    ax.annotate(Label[0], 
                xytext=(BGLoc[0]+.02,BGLoc[1]), 
                xy=(BGLoc[0],BGLoc[1]),
                horizontalalignment='left',verticalalignment = 'center')
    #####################################

    ####### 
    BGLoc = [0, -0.07]
    BGScale = 0.3
    Dir = [1,0]
    Label = ['2$L_{nuc}$']
    ax.annotate("", 
                xytext=(BGLoc[0]-0.5*(Dir[0])*BGScale,BGLoc[1]-0.5*(Dir[1])*BGScale), 
                xy=(BGLoc[0]+0.5*(Dir[0])*BGScale,BGLoc[1]+0.5*(Dir[1])*BGScale),
                arrowprops=dict(arrowstyle="<->"),horizontalalignment='left')
    ax.annotate(Label[0], 
                xytext=(BGLoc[0],BGLoc[1]-.04), 
                xy=(BGLoc[0],BGLoc[1]),
                horizontalalignment='center',verticalalignment = 'center')
    #####################################
    
    ###### Outside Arrows
    OffsetArr = 0.02
    
    ax.annotate('', xy=(-OffsetArr, 0), xycoords='axes fraction', xytext=(-OffsetArr,1), 
            arrowprops=dict(arrowstyle="<->", color='k'))
    
    ax.annotate('', xy=(0, -OffsetArr), xycoords='axes fraction', xytext=(1, -OffsetArr), 
            arrowprops=dict(arrowstyle="<->", color='k'))
    
