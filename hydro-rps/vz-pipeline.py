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

# create a new 'Annotate Time Filter'
annotateTimeFilter1 = AnnotateTimeFilter(registrationName='AnnotateTimeFilter1', Input=grid)
annotateTimeFilter1.Format = 'Time: {TEXT_time:.2f} Myr'
annotateTimeFilter1.Scale = 9.81

# create a new 'Ruler'
ruler1 = Ruler(registrationName='Ruler1')
ruler1.Point1 = [0.0, 130.0, 290.0]
ruler1.Point2 = [0.0, 230.0, 290.0]

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

# get color transfer function/color map for 'vx3'
vx3LUT = GetColorTransferFunction('vx3')
vx3LUT.RGBPoints = [0.0, 0.831373, 0.909804, 0.980392, 0.1, 0.74902, 0.862745, 0.960784, 0.2, 0.694118, 0.827451, 0.941176, 0.4, 0.568627, 0.760784, 0.921569, 0.6, 0.45098, 0.705882, 0.901961, 0.8, 0.345098, 0.643137, 0.858824, 1.0, 0.247059, 0.572549, 0.819608, 1.2, 0.180392, 0.521569, 0.780392, 1.28, 0.14902, 0.490196, 0.74902, 1.44, 0.129412, 0.447059, 0.709804, 1.6, 0.101961, 0.427451, 0.690196, 1.68, 0.094118, 0.403922, 0.658824, 1.76, 0.090196, 0.392157, 0.639216, 1.84, 0.082353, 0.368627, 0.619608, 1.92, 0.070588, 0.352941, 0.6, 2.0, 0.066667, 0.329412, 0.568627, 2.08, 0.07451, 0.313725, 0.541176, 2.16, 0.086275, 0.305882, 0.509804, 2.24, 0.094118, 0.286275, 0.478431, 2.32, 0.101961, 0.278431, 0.45098, 2.4, 0.109804, 0.266667, 0.411765, 2.48, 0.113725, 0.258824, 0.380392, 2.56, 0.113725, 0.25098, 0.34902, 2.64, 0.109804, 0.266667, 0.321569, 2.72, 0.105882, 0.301961, 0.262745, 2.8, 0.094118, 0.309804, 0.243137, 2.88, 0.082353, 0.321569, 0.227451, 2.96, 0.07451, 0.341176, 0.219608, 3.04, 0.070588, 0.360784, 0.211765, 3.12, 0.066667, 0.380392, 0.215686, 3.2, 0.062745, 0.4, 0.176471, 3.4, 0.07451, 0.419608, 0.145098, 3.6, 0.086275, 0.439216, 0.117647, 3.8, 0.121569, 0.470588, 0.117647, 4.0, 0.184314, 0.501961, 0.14902, 4.2, 0.254902, 0.541176, 0.188235, 4.4, 0.32549, 0.580392, 0.231373, 4.6, 0.403922, 0.619608, 0.278431, 4.8, 0.501961, 0.670588, 0.333333, 5.04, 0.592157, 0.729412, 0.4, 5.2, 0.741176, 0.788235, 0.490196, 5.36, 0.858824, 0.858824, 0.603922, 5.6, 0.921569, 0.835294, 0.580392, 6.0, 0.901961, 0.729412, 0.494118, 6.4, 0.858824, 0.584314, 0.388235, 6.8, 0.8, 0.439216, 0.321569, 7.2, 0.678431, 0.298039, 0.203922, 7.6, 0.54902, 0.168627, 0.109804, 7.8, 0.478431, 0.082353, 0.047059, 8.0, 0.45098, 0.007843, 0.0]
vx3LUT.ColorSpace = 'RGB'
vx3LUT.NanColor = [0.25, 0.0, 0.0]
vx3LUT.ScalarRangeInitialized = 1.0

# trace defaults for the display properties.
slice1Display.Representation = 'Surface'
slice1Display.ColorArrayName = ['CELLS', 'vx3']
slice1Display.LookupTable = vx3LUT
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

# get color legend/bar for vx3LUT in view renderView1
vx3LUTColorBar = GetScalarBar(vx3LUT, renderView1)
vx3LUTColorBar.WindowLocation = 'Any Location'
vx3LUTColorBar.Position = [0.7885802469135802, 0.20137931034482756]
vx3LUTColorBar.Title = 'Vz [100 km/s]'
vx3LUTColorBar.ComponentTitle = ''
vx3LUTColorBar.TitleBold = 1
vx3LUTColorBar.TitleFontSize = 18
vx3LUTColorBar.LabelBold = 1
vx3LUTColorBar.ScalarBarThickness = 18
vx3LUTColorBar.ScalarBarLength = 0.586551724137931

# set color bar visibility
vx3LUTColorBar.Visibility = 1

# show color legend
slice1Display.SetScalarBarVisibility(renderView1, True)

# ----------------------------------------------------------------
# setup color maps and opacity mapes used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# get opacity transfer function/opacity map for 'vx3'
vx3PWF = GetOpacityTransferFunction('vx3')
vx3PWF.Points = [0.0, 0.0, 0.5, 0.0, 8.0, 1.0, 0.5, 0.0]
vx3PWF.ScalarRangeInitialized = 1

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
pNG1.Writer.FileName = 'vz_{timestep:06d}.png'
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
options.ExtractsOutputDirectory = 'output/catalyst/vz'
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
