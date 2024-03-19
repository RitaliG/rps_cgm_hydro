# script-version: 2.0
# Catalyst state generated using paraview version 5.10.1

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# ----------------------------------------------------------------
# setup views used in the visualization
# ----------------------------------------------------------------

# get the material library
materialLibrary1 = GetMaterialLibrary()

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [1215, 725]
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.CenterOfRotation = [0.0, -0.000701904296875, 500.0]
renderView1.StereoType = 'Crystal Eyes'
renderView1.CameraPosition = [1964.259862960441, -0.000701904296875, 500.0]
renderView1.CameraFocalPoint = [0.0, -0.000701904296875, 500.0]
renderView1.CameraViewUp = [0.0, 0.0, 1.0]
renderView1.CameraViewAngle = 21.60919540229885
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 508.3878620646296
renderView1.BackEnd = 'OSPRay raycaster'
renderView1.OSPRayMaterialLibrary = materialLibrary1

SetActiveView(None)

# ----------------------------------------------------------------
# setup view layouts
# ----------------------------------------------------------------

# create new layout object 'Layout #1'
layout1 = CreateLayout(name='Layout #1')
layout1.AssignView(0, renderView1)
layout1.SetSize(1215, 725)

# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView1)
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# create a new 'XDMF Reader'
grid = XDMFReader(registrationName='grid', FileNames=['/scratch/ritali/project-ram-pressure/model-for-rps/final-setup-params/check/check-again/temp2e6/overpressurised/output/data.0000.dbl.xmf', '/scratch/ritali/project-ram-pressure/model-for-rps/final-setup-params/check/check-again/temp2e6/overpressurised/output/data.0001.dbl.xmf'])
grid.CellArrayStatus = ['PbykB', 'Temp', 'X', 'Y', 'Z', 'mach', 'ndens', 'prs', 'rho', 'tr1', 'tr2', 'tr3', 'vx1', 'vx2', 'vx3']
grid.GridStatus = ['node_mesh']

# create a new 'Slice'
slice1 = Slice(registrationName='Slice1', Input=grid)
slice1.SliceType = 'Plane'
slice1.HyperTreeGridSlicer = 'Plane'
slice1.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice1.SliceType.Origin = [0.0, -0.00069427490234375, 500.0]

# init the 'Plane' selected for 'HyperTreeGridSlicer'
slice1.HyperTreeGridSlicer.Origin = [0.0, -0.00069427490234375, 500.0]

# create a new 'Annotate Time Filter'
annotateTimeFilter1 = AnnotateTimeFilter(registrationName='AnnotateTimeFilter1', Input=grid)
annotateTimeFilter1.Format = 'Time: {TEXT_time:.2f} Myr'
annotateTimeFilter1.Scale = 9.81

# create a new 'Ruler'
ruler1 = Ruler(registrationName='Ruler1')
ruler1.Point1 = [0.0, 130.0, 290.0]
ruler1.Point2 = [0.0, 230.0, 290.0]

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from slice1
slice1Display = Show(slice1, renderView1, 'GeometryRepresentation')

# get color transfer function/color map for 'tr3'
tr3LUT = GetColorTransferFunction('tr3')
tr3LUT.RGBPoints = [0.0, 0.0, 0.0, 0.34902, 0.03125, 0.039216, 0.062745, 0.380392, 0.0625, 0.062745, 0.117647, 0.411765, 0.09375, 0.090196, 0.184314, 0.45098, 0.125, 0.12549, 0.262745, 0.501961, 0.15625, 0.160784, 0.337255, 0.541176, 0.1875, 0.2, 0.396078, 0.568627, 0.21875, 0.239216, 0.454902, 0.6, 0.25, 0.286275, 0.521569, 0.65098, 0.28125, 0.337255, 0.592157, 0.701961, 0.3125, 0.388235, 0.654902, 0.74902, 0.34375, 0.466667, 0.737255, 0.819608, 0.375, 0.572549, 0.819608, 0.878431, 0.40625, 0.654902, 0.866667, 0.909804, 0.4375, 0.752941, 0.917647, 0.941176, 0.46875, 0.823529, 0.956863, 0.968627, 0.5, 0.941176, 0.984314, 0.988235, 0.5, 0.988235, 0.960784, 0.901961, 0.52, 0.988235, 0.945098, 0.85098, 0.54, 0.980392, 0.898039, 0.784314, 0.5625, 0.968627, 0.835294, 0.698039, 0.59375, 0.94902, 0.733333, 0.588235, 0.625, 0.929412, 0.65098, 0.509804, 0.65625, 0.909804, 0.564706, 0.435294, 0.6875, 0.878431, 0.458824, 0.352941, 0.71875, 0.839216, 0.388235, 0.286275, 0.75, 0.760784, 0.294118, 0.211765, 0.78125, 0.701961, 0.211765, 0.168627, 0.8125, 0.65098, 0.156863, 0.129412, 0.84375, 0.6, 0.094118, 0.094118, 0.875, 0.54902, 0.066667, 0.098039, 0.90625, 0.501961, 0.05098, 0.12549, 0.9375, 0.45098, 0.054902, 0.172549, 0.96875, 0.4, 0.054902, 0.192157, 1.0, 0.34902, 0.070588, 0.211765]
tr3LUT.ColorSpace = 'Lab'
tr3LUT.NanColor = [0.25, 0.0, 0.0]
tr3LUT.ScalarRangeInitialized = 1.0

# trace defaults for the display properties.
slice1Display.Representation = 'Surface'
slice1Display.ColorArrayName = ['CELLS', 'tr3']
slice1Display.LookupTable = tr3LUT
slice1Display.SelectTCoordArray = 'None'
slice1Display.SelectNormalArray = 'None'
slice1Display.SelectTangentArray = 'None'
slice1Display.OSPRayScaleFunction = 'PiecewiseFunction'
slice1Display.SelectOrientationVectors = 'None'
slice1Display.ScaleFactor = 48.0
slice1Display.SelectScaleArray = 'rho'
slice1Display.GlyphType = 'Arrow'
slice1Display.GlyphTableIndexArray = 'rho'
slice1Display.GaussianRadius = 2.4
slice1Display.SetScaleArray = [None, '']
slice1Display.ScaleTransferFunction = 'PiecewiseFunction'
slice1Display.OpacityArray = [None, '']
slice1Display.OpacityTransferFunction = 'PiecewiseFunction'
slice1Display.DataAxesGrid = 'GridAxesRepresentation'
slice1Display.PolarAxes = 'PolarAxesRepresentation'

# show data from annotateTimeFilter1
annotateTimeFilter1Display = Show(annotateTimeFilter1, renderView1, 'TextSourceRepresentation')

# trace defaults for the display properties.
annotateTimeFilter1Display.WindowLocation = 'Upper Center'
annotateTimeFilter1Display.Position = [0.45061728395061723, 0.9375862068965517]
annotateTimeFilter1Display.Bold = 1
annotateTimeFilter1Display.FontSize = 23

# show data from ruler1
ruler1Display = Show(ruler1, renderView1, 'RulerSourceRepresentation')

# trace defaults for the display properties.
ruler1Display.LabelFormat = '%6.3g kpc'
ruler1Display.RulerMode = 1
ruler1Display.Graduation = 25.0
ruler1Display.AxisColor = [0.6666666666666666, 0.6666666666666666, 1.0]
ruler1Display.Color = [0.6666666666666666, 0.6666666666666666, 1.0]

# setup the color legend parameters for each legend in this view

# get color legend/bar for tr3LUT in view renderView1
tr3LUTColorBar = GetScalarBar(tr3LUT, renderView1)
tr3LUTColorBar.WindowLocation = 'Any Location'
tr3LUTColorBar.Position = [0.7680041152263374, 0.24551724137931033]
tr3LUTColorBar.Title = 'Wind Tracer'
tr3LUTColorBar.ComponentTitle = ''
tr3LUTColorBar.TitleBold = 1
tr3LUTColorBar.TitleFontSize = 18
tr3LUTColorBar.LabelBold = 1
tr3LUTColorBar.ScalarBarThickness = 18
tr3LUTColorBar.ScalarBarLength = 0.5286206896551726

# set color bar visibility
tr3LUTColorBar.Visibility = 1

# show color legend
slice1Display.SetScalarBarVisibility(renderView1, True)

# ----------------------------------------------------------------
# setup color maps and opacity mapes used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# get opacity transfer function/opacity map for 'tr3'
tr3PWF = GetOpacityTransferFunction('tr3')
tr3PWF.ScalarRangeInitialized = 1

# ----------------------------------------------------------------
# setup extractors
# ----------------------------------------------------------------

# create extractor
pNG1 = CreateExtractor('PNG', renderView1, registrationName='PNG1')
# trace defaults for the extractor.
pNG1.Trigger = 'TimeValue'

# init the 'TimeValue' selected for 'Trigger'
pNG1.Trigger.Length = 0.5

# init the 'PNG' selected for 'Writer'
pNG1.Writer.FileName = 'wind_tracer_{timestep:06d}.png'
pNG1.Writer.ImageResolution = [1215, 725]
pNG1.Writer.FontScaling = 'Do not scale fonts'
pNG1.Writer.OverrideColorPalette = 'BlackBackground'
pNG1.Writer.Format = 'PNG'

# init the 'PNG' selected for 'Format'
pNG1.Writer.Format.CompressionLevel = '0'

# ----------------------------------------------------------------
# restore active source
SetActiveSource(grid)
# ----------------------------------------------------------------

# ------------------------------------------------------------------------------
# Catalyst options
from paraview import catalyst
options = catalyst.Options()
options.ExtractsOutputDirectory = 'output/catalyst/wind_tracer'
options.GenerateCinemaSpecification = 1
options.GlobalTrigger = 'TimeValue'
options.EnableCatalystLive = 1
options.CatalystLiveTrigger = 'TimeStep'

# init the 'TimeValue' selected for 'GlobalTrigger'
options.GlobalTrigger.Length = 0.001

# ------------------------------------------------------------------------------
if __name__ == '__main__':
    from paraview.simple import SaveExtractsUsingCatalystOptions
    # Code for non in-situ environments; if executing in post-processing
    # i.e. non-Catalyst mode, let's generate extracts using Catalyst options
    SaveExtractsUsingCatalystOptions(options)
