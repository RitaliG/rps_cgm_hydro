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
renderView1.ViewSize = [1611, 725]
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
layout1.SetSize(1611, 725)

# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView1)
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# create a new 'Ruler'
ruler1 = Ruler(registrationName='Ruler1')
ruler1.Point1 = [0.0, 130.0, 290.0]
ruler1.Point2 = [0.0, 230.0, 290.0]

# create a new 'XDMF Reader'
grid = XDMFReader(registrationName='grid', FileNames=['/scratch/ritali/project-ram-pressure/model-for-rps/final-setup-params/check/check-again/temp2e6/overpressurised/output/data.0000.dbl.xmf', '/scratch/ritali/project-ram-pressure/model-for-rps/final-setup-params/check/check-again/temp2e6/overpressurised/output/data.0001.dbl.xmf'])
grid.CellArrayStatus = ['PbykB', 'Temp', 'X', 'Y', 'Z', 'mach', 'ndens', 'prs', 'rho', 'tr1', 'tr2', 'tr3', 'vx1', 'vx2', 'vx3']
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

# get color transfer function/color map for 'vx2'
vx2LUT = GetColorTransferFunction('vx2')
vx2LUT.RGBPoints = [-1.8, 0.831373, 0.909804, 0.980392, -1.7551758539493045, 0.74902, 0.862745, 0.960784, -1.7103517078986092, 0.694118, 0.827451, 0.941176, -1.6207034157972182, 0.568627, 0.760784, 0.921569, -1.5310551236958276, 0.45098, 0.705882, 0.901961, -1.4414068315944366, 0.345098, 0.643137, 0.858824, -1.3517585394930458, 0.247059, 0.572549, 0.819608, -1.262110247391655, 0.180392, 0.521569, 0.780392, -1.2262509305510987, 0.14902, 0.490196, 0.74902, -1.154532296869986, 0.129412, 0.447059, 0.709804, -1.0828136631888732, 0.101961, 0.427451, 0.690196, -1.046954346348317, 0.094118, 0.403922, 0.658824, -1.0110950295077605, 0.090196, 0.392157, 0.639216, -0.9752357126672042, 0.082353, 0.368627, 0.619608, -0.9393763958266479, 0.070588, 0.352941, 0.6, -0.9035170789860916, 0.066667, 0.329412, 0.568627, -0.8676577621455352, 0.07451, 0.313725, 0.541176, -0.8317984453049788, 0.086275, 0.305882, 0.509804, -0.7959391284644224, 0.094118, 0.286275, 0.478431, -0.7600798116238663, 0.101961, 0.278431, 0.45098, -0.72422049478331, 0.109804, 0.266667, 0.411765, -0.6883611779427536, 0.113725, 0.258824, 0.380392, -0.6525018611021971, 0.113725, 0.25098, 0.34902, -0.6166425442616408, 0.109804, 0.266667, 0.321569, -0.5807832274210845, 0.105882, 0.301961, 0.262745, -0.5449239105805284, 0.094118, 0.309804, 0.243137, -0.5090645937399718, 0.082353, 0.321569, 0.227451, -0.4732052768994155, 0.07451, 0.341176, 0.219608, -0.4373459600588594, 0.070588, 0.360784, 0.211765, -0.40148664321830285, 0.066667, 0.380392, 0.215686, -0.3656273263777463, 0.062745, 0.4, 0.176471, -0.27597903427635573, 0.07451, 0.419608, 0.145098, -0.1863307421749647, 0.086275, 0.439216, 0.117647, -0.09668245007357412, 0.121569, 0.470588, 0.117647, -0.007034157972183097, 0.184314, 0.501961, 0.14902, 0.08261413412920793, 0.254902, 0.541176, 0.188235, 0.17226242623059873, 0.32549, 0.580392, 0.231373, 0.2619107183319891, 0.403922, 0.619608, 0.278431, 0.3515590104333801, 0.501961, 0.670588, 0.333333, 0.4591369609550491, 0.592157, 0.729412, 0.4, 0.5308555946361622, 0.741176, 0.788235, 0.490196, 0.6025742283172748, 0.858824, 0.858824, 0.603922, 0.7101521788389433, 0.921569, 0.835294, 0.580392, 0.8894487630417254, 0.901961, 0.729412, 0.494118, 1.0687453472445074, 0.858824, 0.584314, 0.388235, 1.2480419314472886, 0.8, 0.439216, 0.321569, 1.4273385156500706, 0.678431, 0.298039, 0.203922, 1.6066350998528518, 0.54902, 0.168627, 0.109804, 1.6962833919542428, 0.478431, 0.082353, 0.047059, 1.7859316840556339, 0.45098, 0.007843, 0.0]
vx2LUT.ColorSpace = 'RGB'
vx2LUT.NanColor = [0.25, 0.0, 0.0]
vx2LUT.ScalarRangeInitialized = 1.0

# trace defaults for the display properties.
slice1Display.Representation = 'Surface'
slice1Display.ColorArrayName = ['CELLS', 'vx2']
slice1Display.LookupTable = vx2LUT
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

# get color legend/bar for vx2LUT in view renderView1
vx2LUTColorBar = GetScalarBar(vx2LUT, renderView1)
vx2LUTColorBar.WindowLocation = 'Any Location'
vx2LUTColorBar.Position = [0.7405338299193048, 0.2358620689655173]
vx2LUTColorBar.Title = 'Vy [100 km/s]'
vx2LUTColorBar.ComponentTitle = ''
vx2LUTColorBar.TitleBold = 1
vx2LUTColorBar.TitleFontSize = 18
vx2LUTColorBar.LabelBold = 1
vx2LUTColorBar.ScalarBarThickness = 18
vx2LUTColorBar.ScalarBarLength = 0.5506896551724138

# set color bar visibility
vx2LUTColorBar.Visibility = 1

# show color legend
slice1Display.SetScalarBarVisibility(renderView1, True)

# ----------------------------------------------------------------
# setup color maps and opacity mapes used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# get opacity transfer function/opacity map for 'vx2'
vx2PWF = GetOpacityTransferFunction('vx2')
vx2PWF.Points = [-1.8, 0.0, 0.5, 0.0, 1.7859316840556339, 1.0, 0.5, 0.0]
vx2PWF.ScalarRangeInitialized = 1

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
pNG1.Writer.FileName = 'vy_{timestep:06d}.png'
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
options.ExtractsOutputDirectory = 'output/catalyst/vy'
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
