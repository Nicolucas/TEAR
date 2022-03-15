# -*- coding:utf-8 -*-
# @Script: RecordViews.py
# @Author: J.N. Hayek
# @Email: jhayek@geophysik.uni-muenchen.de
# @Description: Saves several views of the stress and velocity field (inc. mesh visualization)).



# trace generated using paraview version 5.10.0
#import paraview
#paraview.compatibility.major = 5
#paraview.compatibility.minor = 10

#### import the simple module from the paraview
from paraview.simple import *
import os, sys
import numpy as np
from glob import glob
import json

##### Paths 
Path='/import/freenas-m-03-geodynamics/jhayek/TEAR/Results/T2/Runs/'
fNameFolder='TEAR33_TPV3_T0_P3_050x050_A12phi65_Delta1.001_3s_Stressfull_KV0.3_NES_F/'
fNamePattern='Sigma-Aligned-step-{timestep:04}.vtu'
#####################################################



#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')

# destroy renderView1
Delete(renderView1)
del renderView1


#================================================================
# Extra functions
#================================================================
def GetListPatternFiles(path,fname,ToFormat):
    FindFilename_ = fname.replace(ToFormat,"*")
    PathNFile_ = os.path.join(path,FindFilename_)
    print(PathNFile_)

    FileList_ = glob(PathNFile_)
    list_ = [int(i.replace(PathNFile_.split('*')[0],'').replace(PathNFile_.split('*')[1],'')) for i in FileList_]

    return sorted(list_)

def CreateFolder(PathFolder):
    try:
        os.makedirs(PathFolder)
    except:
        pass
#================================================================
# END: Extra functions
#================================================================

#================================================================
# Get Layout
#================================================================
layout1 = GetLayoutByName("Layout #1")

# split cell
layout1.SplitHorizontal(0, 0.5)
layout1.SplitVertical(1, 0.5)
layout1.SplitVertical(2, 0.5)
layout1.SplitVertical(6, 0.5)

# resize frame
layout1.SetSplitFraction(2, 0.333)

#================================================================
# Display view on each window
#================================================================

def Display(layout,hint):
    # set active view
    SetActiveView(None)

    # get the material library
    materialLibrary1 = GetMaterialLibrary()

    # Create a new 'Render View'
    renderView1 = CreateView('RenderView')
    renderView1.AxesGrid = 'GridAxes3DActor'
    renderView1.StereoType = 'Crystal Eyes'
    renderView1.CameraFocalDisk = 1.0
    renderView1.BackEnd = 'OSPRay raycaster'
    renderView1.OSPRayMaterialLibrary = materialLibrary1

    # assign view to a particular cell in the layout
    AssignViewToLayout(view=renderView1, layout=layout, hint=hint)
    return renderView1

renderViewList = []

renderViewList.append(Display(layout1,1))
renderViewList.append(Display(layout1,2))
renderViewList.append(Display(layout1,3))
renderViewList.append(Display(layout1,4))
renderViewList.append(Display(layout1,5))

#================================================================
# END: Display view on each window
#================================================================

#================================================================
# Load a stress field and the velocity field
#================================================================

def VisualizeStressField(renderView2,sigmaAlignedstep0,StressComponent):
    # set active view
    SetActiveView(renderView2)

    # set active source
    SetActiveSource(sigmaAlignedstep0)

    # show data in view
    sigmaAlignedstep0Display = Show(sigmaAlignedstep0, renderView2, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    sigmaAlignedstep0Display.Representation = 'Surface'
    sigmaAlignedstep0Display.ColorArrayName = [None, StressComponent]
    sigmaAlignedstep0Display.SelectTCoordArray = 'None'
    sigmaAlignedstep0Display.SelectNormalArray = 'None'
    sigmaAlignedstep0Display.SelectTangentArray = 'None'
    sigmaAlignedstep0Display.OSPRayScaleArray = StressComponent
    sigmaAlignedstep0Display.OSPRayScaleFunction = 'PiecewiseFunction'
    sigmaAlignedstep0Display.SelectOrientationVectors = 'None'
    sigmaAlignedstep0Display.ScaleFactor = 500.0
    sigmaAlignedstep0Display.SelectScaleArray = 'None'
    sigmaAlignedstep0Display.GlyphType = 'Arrow'
    sigmaAlignedstep0Display.GlyphTableIndexArray = 'None'
    sigmaAlignedstep0Display.GaussianRadius = 100.0
    sigmaAlignedstep0Display.SetScaleArray = ['POINTS', StressComponent]
    sigmaAlignedstep0Display.ScaleTransferFunction = 'PiecewiseFunction'
    sigmaAlignedstep0Display.OpacityArray = ['POINTS', StressComponent]
    sigmaAlignedstep0Display.OpacityTransferFunction = 'PiecewiseFunction'
    sigmaAlignedstep0Display.DataAxesGrid = 'GridAxesRepresentation'
    sigmaAlignedstep0Display.PolarAxes = 'PolarAxesRepresentation'
    sigmaAlignedstep0Display.ScalarOpacityUnitDistance = 250.4710267850061
    sigmaAlignedstep0Display.OpacityArrayName = ['POINTS', StressComponent]

    # init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
    sigmaAlignedstep0Display.OSPRayScaleFunction.Points = [-1.014485005353776, 0.0, 0.5, 0.0, 1.0144850053537764, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    sigmaAlignedstep0Display.ScaleTransferFunction.Points = [-463400.0, 0.0, 0.5, 0.0, 469440.0, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    sigmaAlignedstep0Display.OpacityTransferFunction.Points = [-463400.0, 0.0, 0.5, 0.0, 469440.0, 1.0, 0.5, 0.0]

    #changing interaction mode based on data extents
    renderView2.InteractionMode = '2D'
    #renderView2.CameraPosition = [0.0, 0.0, 10000.0]
    #renderView2.CameraParallelScale = 2420.416891764817
    renderView2.CameraPosition = [1910.145574444335, 0.0, 10000.0]
    renderView2.CameraFocalPoint = [1910.145574444335, 0.0, 0.0]
    renderView2.CameraParallelScale = 526.7532364066545

    # get display properties
    stressFieldDisplay = GetDisplayProperties(sigmaAlignedstep0, view=renderView2)

    # set scalar coloring
    ColorBy(stressFieldDisplay, ('POINTS', StressComponent))

    # rescale color and/or opacity maps used to include current data range
    stressFieldDisplay.RescaleTransferFunctionToDataRange(True, False)

    # show color bar/color legend
    stressFieldDisplay.SetScalarBarVisibility(renderView2, True)

    # get color transfer function/color map for 'sigma_xx'
    sigma_xxLUT = GetColorTransferFunction(StressComponent)

    # get opacity transfer function/opacity map for 'sigma_xx'
    sigma_xxPWF = GetOpacityTransferFunction(StressComponent)


def VisualizeVelocityField(renderView1,Motionstep0,MotionField):
    # set active view
    SetActiveView(renderView1)

    # set active source
    SetActiveSource(Motionstep0)

    # show data in view
    step0Display = Show(Motionstep0, renderView1, 'StructuredGridRepresentation')

    # trace defaults for the display properties.
    step0Display.Representation = 'Surface'
    step0Display.ColorArrayName = [None, MotionField]
    step0Display.SelectTCoordArray = 'None'
    step0Display.SelectNormalArray = 'None'
    step0Display.SelectTangentArray = 'None'
    step0Display.OSPRayScaleArray = MotionField
    step0Display.OSPRayScaleFunction = 'PiecewiseFunction'
    step0Display.SelectOrientationVectors = 'None'
    step0Display.ScaleFactor = 2000.0
    step0Display.SelectScaleArray = 'None'
    step0Display.GlyphType = 'Arrow'
    step0Display.GlyphTableIndexArray = 'None'
    step0Display.GaussianRadius = 100.0
    step0Display.SetScaleArray = ['POINTS', MotionField]
    step0Display.ScaleTransferFunction = 'PiecewiseFunction'
    step0Display.OpacityArray = ['POINTS', MotionField]
    step0Display.OpacityTransferFunction = 'PiecewiseFunction'
    step0Display.DataAxesGrid = 'GridAxesRepresentation'
    step0Display.PolarAxes = 'PolarAxesRepresentation'
    step0Display.ScalarOpacityUnitDistance = 250.4710267850061

    # init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
    step0Display.OSPRayScaleFunction.Points = [-1.014485005353776, 0.0, 0.5, 0.0, 1.0144850053537764, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    step0Display.ScaleTransferFunction.Points = [-0.7560609106722755, 0.0, 0.5, 0.0, 0.7560609106722715, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    step0Display.OpacityTransferFunction.Points = [-0.7560609106722755, 0.0, 0.5, 0.0, 0.7560609106722715, 1.0, 0.5, 0.0]

    #changing interaction mode based on data extents
    renderView1.InteractionMode = '2D'
    renderView1.CameraPosition = [1910.145574444335, 0.0, 10000.0]
    renderView1.CameraFocalPoint = [1910.145574444335, 0.0, 0.0]
    renderView1.CameraParallelScale = 526.7532364066545

    # set scalar coloring
    ColorBy(step0Display, ('POINTS', MotionField))

    # rescale color and/or opacity maps used to include current data range
    step0Display.RescaleTransferFunctionToDataRange(True, False)

    # show color bar/color legend
    step0Display.SetScalarBarVisibility(renderView1, True)

    # get color transfer function/color map for 'velo_x'
    velo_xLUT = GetColorTransferFunction(MotionField)

    # get opacity transfer function/opacity map for 'velo_x'
    velo_xPWF = GetOpacityTransferFunction(MotionField)




### Get available timesteps
TimeStepList = GetListPatternFiles(Path+fNameFolder,fNamePattern,"{timestep:04}")

freq = int(TimeStepList[1])-int(TimeStepList[0])
maxtimestep = int(TimeStepList[-1]) 

TSList = np.arange(freq, maxtimestep+1, freq).tolist()
### Produce paths of filenames with the available filenames
StressFname = [os.path.join(Path+fNameFolder,fNamePattern.format(timestep=i)) for i in TSList]

fNamePattern='step-{timestep:04}.vts'
DispFname=[os.path.join(Path+fNameFolder,fNamePattern.format(timestep=i)) for i in TSList]

# create a new 'XML Unstructured Grid Reader'
sigmaAlignedstep = XMLUnstructuredGridReader(registrationName='StressField', FileName=StressFname)
sigmaAlignedstep.PointArrayStatus = ['sigma_xx', 'sigma_yy', 'sigma_xy']

# create a new 'XML Structured Grid Reader'
step0 = XMLStructuredGridReader(registrationName='MotionField', FileName=DispFname)
step0.PointArrayStatus = ['disp._x', 'disp._y', 'velo._x', 'velo._y']

# get animation scene
animationScene1 = GetAnimationScene()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

LastFrameIdx = len(TSList)-2

animationScene1.AnimationTime = LastFrameIdx

# get the time-keeper
timeKeeper1 = GetTimeKeeper()

VisualizeStressField(renderViewList[1],sigmaAlignedstep,StressComponent='sigma_xy')
VisualizeStressField(renderViewList[3],sigmaAlignedstep,StressComponent='sigma_xx')
VisualizeStressField(renderViewList[4],sigmaAlignedstep,StressComponent='sigma_yy')

VisualizeVelocityField(renderViewList[0],step0,'velo._x')
VisualizeVelocityField(renderViewList[2],step0,'velo._x')
#================================================================
# END: Load a stress field and the velocity field
#================================================================

#================================================================
# Ruler as scale
#================================================================

ruler1 = Ruler(registrationName='Ruler1')
# Properties modified on ruler1
ruler1.Point1 = [-1500.0, -50.0, 0.0]
ruler1.Point2 = [1500.0, -50.0, 0.0]


ruler2 = Ruler(registrationName='Ruler2')
# Properties modified on ruler1
ruler2.Point1 = [1500.0, -50.0, 0.0]
ruler2.Point2 = [1500.0, 50.0, 0.0]

# show data in view
ruler1Display = Show(ruler1, renderViewList[2], 'RulerSourceRepresentation')
ruler2Display = Show(ruler2, renderViewList[2], 'RulerSourceRepresentation')

ruler1Display = Show(ruler1, renderViewList[1], 'RulerSourceRepresentation')

#================================================================
# END: Ruler as scale
#================================================================

jsonPattern="step-{timestep:04}_wavefield.json"

StampDataList = [os.path.join(Path+fNameFolder,jsonPattern.format(timestep=i)) for i in TSList]
#Prelude to extract frequency and timestamp
with open(StampDataList[0], 'rb') as f:
    StampData = json.load(f)
    fq = StampData['se2wave']['time']

#================================================================
# Timestamp
#================================================================
# create a new 'Annotate Time Filter'
annotateTimeFilter1 = AnnotateTimeFilter(registrationName='AnnotateTimeFilter1', Input=step0)
# Properties modified on annotateTimeFilter1
annotateTimeFilter1.Scale = fq
annotateTimeFilter1.Shift = fq
# show data in view
annotateTimeFilter1Display = Show(annotateTimeFilter1, renderViewList[2], 'TextSourceRepresentation')
#================================================================
# END: Timestamp
#================================================================


def RescaleFieldColors(FieldSource,renderView):
    
    # get display properties
    SetActiveView(renderView)
    Field = FindSource(FieldSource)
    SetActiveSource(Field)

    FieldDisplay = GetDisplayProperties(Field, view=renderView)
    # rescale color and/or opacity maps used to exactly fit the current data range
    FieldDisplay.RescaleTransferFunctionToDataRange(False, True)


RescaleFieldColors('MotionField',renderViewList[0])
RescaleFieldColors('MotionField',renderViewList[2])
RescaleFieldColors('StressField',renderViewList[1])
RescaleFieldColors('StressField',renderViewList[3])
RescaleFieldColors('StressField',renderViewList[4])

#================================================================
# Representation change
#================================================================
def ChangeRepresentation(FieldSource,renderView):
    
    # get display properties
    SetActiveView(renderView)
    Field = FindSource(FieldSource)
    SetActiveSource(Field)

    FieldDisplay = GetDisplayProperties(Field, view=renderView)
    # rescale color and/or opacity maps used to exactly fit the current data range
    FieldDisplay.SetRepresentationType('Surface With Edges')

ChangeRepresentation('MotionField',renderViewList[2])
#================================================================
# END: Representation
#================================================================
def SaveMyStuff():
    DefaultPath="/import/freenas-m-03-geodynamics/jhayek/SharedWolfel/Video/"

    CreateFolder(DefaultPath + fNameFolder)
    # save animation
    SaveAnimation(DefaultPath + fNameFolder + 'Fig_.png', layout1, SaveAllViews=1, ImageResolution=[1523, 632], FrameWindow=[0, LastFrameIdx])

#SaveMyStuff()

#SaveMyStuff('/import/freenas-m-03-geodynamics/jhayek/SharedWolfel/Video/OutputVideo/')
#================================================================
# Save figures into video
#================================================================

#================================================================
# End: Save figures into video
#================================================================

#================================================================
# Clean up layout function
#================================================================
def clean(renderViewList=renderViewList):
    layout1 = GetLayoutByName("Layout #1")
    for render in renderViewList:
        # destroy render
        Delete(render)
        del render
    
    for idx in [5,5,2,1]:
        layout1.Collapse(idx)

    stressField = FindSource('StressField')
    Delete(stressField)
    del stressField

    motionField = FindSource('MotionField')
    Delete(motionField)
    del motionField

    R1 = FindSource('Ruler1')
    Delete(R1)
    del R1

    R2 = FindSource('Ruler2')
    Delete(R2)
    del R2

    Time = FindSource('AnnotateTimeFilter1')
    Delete(Time)
    del Time

