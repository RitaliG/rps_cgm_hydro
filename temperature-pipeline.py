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
renderView1.CameraPosition = [1101.4095260343797, -0.000701904296875, 500.0]
renderView1.CameraFocalPoint = [0.0, -0.000701904296875, 500.0]
renderView1.CameraViewUp = [0.0, 0.0, 1.0]
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 505.011386031381
renderView1.Background = [0.0, 0.0, 0.0]
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
grid = XDMFReader(registrationName='grid')
grid.CellArrayStatus = ['PbykB', 'Temp', 'X', 'Y', 'Z', 'ndens', 'prs', 'rho', 'tr1', 'tr2', 'tr3', 'vx1', 'vx2', 'vx3']
grid.GridStatus = ['node_mesh']

# create a new 'Annotate Time Filter'
annotateTimeFilter1 = AnnotateTimeFilter(registrationName='AnnotateTimeFilter1', Input=grid)
annotateTimeFilter1.Format = 'Time: {TEXT_time:.2f} Myr'
annotateTimeFilter1.Scale = 9.81

# create a new 'Slice'
slice1 = Slice(registrationName='Slice1', Input=grid)
slice1.SliceType = 'Plane'
slice1.HyperTreeGridSlicer = 'Plane'
slice1.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice1.SliceType.Origin = [0.0, -0.00069427490234375, 500.0]

# init the 'Plane' selected for 'HyperTreeGridSlicer'
slice1.HyperTreeGridSlicer.Origin = [0.0, -0.00069427490234375, 500.0]

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from slice1
slice1Display = Show(slice1, renderView1, 'GeometryRepresentation')

# get color transfer function/color map for 'Temp'
tempLUT = GetColorTransferFunction('Temp')
tempLUT.RGBPoints = [40721.45439232455, 0.0, 0.0, 0.34902, 48363.77216378239, 0.039216, 0.062745, 0.380392, 57440.346687399455, 0.062745, 0.117647, 0.411765, 68220.3492398268, 0.090196, 0.184314, 0.45098, 81023.467280445, 0.12549, 0.262745, 0.501961, 96229.38497525085, 0.160784, 0.337255, 0.541176, 114289.042959162, 0.2, 0.396078, 0.568627, 135738.0112517663, 0.239216, 0.454902, 0.6, 161212.37190837497, 0.286275, 0.521569, 0.65098, 191467.58241594592, 0.337255, 0.592157, 0.701961, 227400.8792392356, 0.388235, 0.654902, 0.74902, 270077.8859078067, 0.466667, 0.737255, 0.819608, 320764.21472272265, 0.572549, 0.819608, 0.878431, 380962.99925054726, 0.654902, 0.866667, 0.909804, 452459.47065332375, 0.752941, 0.917647, 0.941176, 537373.8997924278, 0.823529, 0.956863, 0.968627, 638224.4751362028, 0.941176, 0.984314, 0.988235, 638224.4751362028, 0.988235, 0.960784, 0.901961, 712491.1295875724, 0.988235, 0.945098, 0.85098, 795399.7841161406, 0.980392, 0.898039, 0.784314, 900258.4918262891, 0.968627, 0.835294, 0.698039, 1069212.7095392642, 0.94902, 0.733333, 0.588235, 1269875.073237175, 0.929412, 0.65098, 0.509804, 1508196.346004903, 0.909804, 0.564706, 0.435294, 1791244.0885259416, 0.878431, 0.458824, 0.352941, 2127412.251845296, 0.839216, 0.388235, 0.286275, 2526670.1050362843, 0.760784, 0.294118, 0.211765, 3000857.8798710075, 0.701961, 0.211765, 0.168627, 3564037.9000148815, 0.65098, 0.156863, 0.129412, 4232911.60769284, 0.6, 0.094118, 0.094118, 5027314.855003637, 0.54902, 0.066667, 0.098039, 5970806.1480442295, 0.501961, 0.05098, 0.12549, 7091365.2885775035, 0.45098, 0.054902, 0.172549, 8422223.13188209, 0.4, 0.054902, 0.192157, 10002847.067752536, 0.34902, 0.070588, 0.211765]
tempLUT.UseLogScale = 1
tempLUT.ColorSpace = 'Lab'
tempLUT.NanColor = [0.25, 0.0, 0.0]
tempLUT.ScalarRangeInitialized = 1.0

# trace defaults for the display properties.
slice1Display.Representation = 'Surface'
slice1Display.ColorArrayName = ['CELLS', 'Temp']
slice1Display.LookupTable = tempLUT
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
annotateTimeFilter1Display.Bold = 1
annotateTimeFilter1Display.FontSize = 23

# setup the color legend parameters for each legend in this view

# get color legend/bar for tempLUT in view renderView1
tempLUTColorBar = GetScalarBar(tempLUT, renderView1)
tempLUTColorBar.WindowLocation = 'Any Location'
tempLUTColorBar.Position = [0.8469135802469137, 0.05517241379310339]
tempLUTColorBar.Title = 'Temperature (K)'
tempLUTColorBar.ComponentTitle = ''
tempLUTColorBar.TitleBold = 1
tempLUTColorBar.TitleFontSize = 20
tempLUTColorBar.LabelBold = 1
tempLUTColorBar.LabelFontSize = 20
tempLUTColorBar.AutomaticLabelFormat = 0
tempLUTColorBar.LabelFormat = '%-#6.1e'
tempLUTColorBar.ScalarBarThickness = 18
tempLUTColorBar.ScalarBarLength = 0.8941379310344832

# set color bar visibility
tempLUTColorBar.Visibility = 1

# show color legend
slice1Display.SetScalarBarVisibility(renderView1, True)

# ----------------------------------------------------------------
# setup color maps and opacity mapes used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# get opacity transfer function/opacity map for 'Temp'
tempPWF = GetOpacityTransferFunction('Temp')
tempPWF.Points = [40721.454392324544, 0.0, 0.5, 0.0, 10002847.067752533, 1.0, 0.5, 0.0]
tempPWF.ScalarRangeInitialized = 1

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
pNG1.Writer.FileName = 'temperature_{timestep:06d}{camera}.png'
pNG1.Writer.ImageResolution = [1215, 725]
pNG1.Writer.FontScaling = 'Do not scale fonts'
pNG1.Writer.OverrideColorPalette = 'BlackBackground'
pNG1.Writer.Format = 'PNG'

# init the 'PNG' selected for 'Format'
pNG1.Writer.Format.CompressionLevel = '0'

# ----------------------------------------------------------------
# restore active source
SetActiveSource(pNG1)
# ----------------------------------------------------------------

# ------------------------------------------------------------------------------
# Catalyst options
from paraview import catalyst
options = catalyst.Options()
options.ExtractsOutputDirectory = 'output/catalyst/temperature'
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
