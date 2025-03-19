# ====================================================================================================
# slicer.mrmlScene.GetRootDirectory()
# slicer.app.slicerHome -> .slicerrc.py (getSlicerRCFileName())
# shortCuts = [
#   ("Ctrl+b", lambda: slicer.app.layoutManager().setLayout(slicer.vtkMRMLLayoutNode.SlicerLayoutOneUpRedSliceView)),
#   ("Ctrl+n", lambda: slicer.app.layoutManager().setLayout(slicer.vtkMRMLLayoutNode.SlicerLayoutOneUpYellowSliceView)),
#   ("Ctrl+m", lambda: slicer.app.layoutManager().setLayout(slicer.vtkMRMLLayoutNode.SlicerLayoutOneUpGreenSliceView)),
#   ("Ctrl+,", lambda: slicer.app.layoutManager().setLayout(slicer.vtkMRMLLayoutNode.SlicerLayoutFourUpView))
#   ("Ctrl+e", lambda: slicer.util.restart()) # must be free key
#   ]

# for (shortCutKey, callBack) in shortCuts:
#   shortCut = qt.QShortcut(slicer.util.mainWindow())
#   shortCut.setKey(qt.QKeySequence(shortCutKey))
#   shortCut.connect( "activated()", callBack)

# Startup Module: Application Settings -> Modules -> Change Startup Module
# Favorite Module: Drag-and-drop modules from the Modules list to the Favorite Modules list to add a module.

# ====================================================================================================
#from __main__ import vtk, qt, ctk, slicer
from __main__ import qt, ctk, slicer

# ====================================================================================================
# Python modules firstly necessary to be installed from CLI
# import pip
# pip.main(['install', 'moduleName'])
# import 'moduleName'

# ====================================================================================================
# modules must be in subfolders including CMakeList.txt
import PelvicFloorLib

# ====================================================================================================
# version ('development' or 'release')

# ====================================================================================================
class PelvicFloor:
  def __init__(self, parent): 
    parent.title = "Pelvic Floor Model"
    parent.categories = ["Registration"]
    parent.dependencies = []
    parent.contributors = ["Ludek Hyncik (UWB)", "Adam Wittek (UWA)"]
    parent.helpText = """
    A module to morph pelvic floor based on chosen landmarks.
    """
    parent.acknowledgementText = """
    This file was originally developed by Ludek Hyncik, University of West Bohemia
    and Adam Wittek, The University of Western Australia and was partially funded
    by Gledden Fellowship Grant at The University of Western Australia."""
    self.parent = parent

# ====================================================================================================
class PelvicFloorWidget:
  def __init__(self, parent = None): # constructor 
    if not parent:
      self.parent = slicer.qMRMLWidget()
      self.parent.setLayout(qt.QVBoxLayout())
      self.parent.setMRMLScene(slicer.mrmlScene)
    else:
      self.parent = parent
    self.layout = self.parent.layout()
    if not parent:
      self.setup()
      self.parent.show()

  # ====================================================================================================
  # Instantiate and connect widgets
  def setup(self):

    # ====================================================================================================
	  # Collapsible Button for 'Female Pelvic Floor'
    pelvisFloorCollapsibleButton = ctk.ctkCollapsibleButton()
    pelvisFloorCollapsibleButton.text = "Model Manipulation"
    self.layout.addWidget(pelvisFloorCollapsibleButton)

    # ----------------------------------------------------------------------------------------------------
    # Layout within the 'Female Pelvic Floor' collapsible button
    #PelvisFloorFormLayout = qt.QFormLayout(PelvisFloorCollapsibleButton)
    pelvisFloorLayout = qt.QGridLayout(pelvisFloorCollapsibleButton)

    rowHeight = 1 # height of single row
    totalRowWidth = 60 # total row width

    noOfColumns = 4 # number of columns in row
    colWidth = totalRowWidth / noOfColumns # width of single column

    row = 0 # row position
    column = 0 # column position

    OpenButton = qt.QPushButton("Open")
    OpenButton.toolTip = "Open VTK file."
    OpenButton.connect('clicked(bool)', self.onOpenButtonClicked)
    self.OpenButton = OpenButton
    pelvisFloorLayout.addWidget(OpenButton, row, column, rowHeight, colWidth)

    column += colWidth # column position

    TemplateButton = qt.QPushButton("Template")
    TemplateButton.setStyleSheet('background-color: green')
    TemplateButton.toolTip = "Read template model as VTK and display as Model Node 'Model'."
    TemplateButton.connect('clicked(bool)', self.onTemplateButtonClicked)
    # by doing this, one is able to refer to TemplateButton instead
    # of self.TemplateButton within the scope of PelvicFloorWidget
    self.TemplateButton = TemplateButton # set local var as instance attribute
    # (item, row, column[, rowSpan = 1[, columnSpan = 1]])
    pelvisFloorLayout.addWidget(TemplateButton, row, column, rowHeight, colWidth)
    self.layout.addStretch(1) # add vertical spacer, otherwise spread

    column += colWidth # column position

    ImagesButton = qt.QPushButton("MRI")
    ImagesButton.setStyleSheet('background-color: grey')
    ImagesButton.toolTip = "Read MRI from directory."
    ImagesButton.connect('clicked(bool)', self.onImagesButtonClicked)
    self.ImagesButton = ImagesButton
    pelvisFloorLayout.addWidget(ImagesButton, row, column, rowHeight, colWidth)

    column += colWidth # column position

    ExportButton = qt.QPushButton("Export")
    ExportButton.toolTip = "Export model as unstructured VTK."
    ExportButton.connect('clicked(bool)', self.onExportButtonClicked)
    self.ExportButton = ExportButton
    pelvisFloorLayout.addWidget(ExportButton, row, column, rowHeight, colWidth)

    noOfColumns = 3 # number of columns in row
    colWidth = totalRowWidth / noOfColumns # width of single column

    row += 1 # row position
    column = 0 # column position

    RotLButton = qt.QPushButton("Rotate L")
    RotLButton.toolTip = "Rotate model by 90 degrees around axis L at LPS = (0, 0, 0)."
    RotLButton.connect('clicked(bool)', self.onRotLButtonClicked)
    self.RotLButton = RotLButton
    pelvisFloorLayout.addWidget(RotLButton, row, column, rowHeight, colWidth)

    column += colWidth # column position

    RotPButton = qt.QPushButton("Rotate P")
    RotPButton.toolTip = "Rotate model by 90 degrees around axis P at LPS = (0, 0, 0)."
    RotPButton.connect('clicked(bool)', self.onRotPButtonClicked)
    self.RotPButton = RotPButton
    pelvisFloorLayout.addWidget(RotPButton, row, column, rowHeight, colWidth)

    column += colWidth # column position

    RotSButton = qt.QPushButton("Rotate S")
    RotSButton.toolTip = "Rotate model by 90 degrees around axis S at LPS = (0, 0, 0)."
    RotSButton.connect('clicked(bool)', self.onRotSButtonClicked)
    self.RotSButton = RotSButton
    pelvisFloorLayout.addWidget(RotSButton, row, column, rowHeight, colWidth)

    noOfColumns = 4 # number of columns in row
    colWidth = totalRowWidth / noOfColumns # width of single column

    row += 1 # row position
    column = 0 # column position

    LandmarksButton = qt.QPushButton("Landmarks")
    LandmarksButton.setStyleSheet('background-color: yellow')
    LandmarksButton.toolTip = "Read landmarks from file."
    LandmarksButton.connect('clicked(bool)', self.onLandmarksButtonClicked)
    self.LandmarksButton = LandmarksButton
    pelvisFloorLayout.addWidget(LandmarksButton, row, column, rowHeight, colWidth)

    column += colWidth # column position

    AlignButton = qt.QPushButton("Align")
    AlignButton.setStyleSheet('background-color: blue')
    AlignButton.toolTip = "Align model to PSI. Name 'PSI' must exist as a Markups Node."
    AlignButton.connect('clicked(bool)', self.onAlignButtonClicked)
    self.AlignButton = AlignButton
    pelvisFloorLayout.addWidget(AlignButton, row, column, rowHeight, colWidth)

    column += colWidth # column position

    SegmentationButton = qt.QPushButton("Segmentation")
    SegmentationButton.toolTip = "Show model segmentation in MRML scene as Model Segmentation."
    SegmentationButton.connect('clicked(bool)', self.onSegmentationButtonClicked)
    self.SegmentationButton = SegmentationButton
    pelvisFloorLayout.addWidget(SegmentationButton, row, column, rowHeight, colWidth)

    column += colWidth # column position

    LabelsButton = qt.QPushButton("Labels")
    LabelsButton.setStyleSheet('background-color: blue')
    LabelsButton.toolTip = "Show model labeled segmentation in MRML scene as Model Segmentation."
    LabelsButton.connect('clicked(bool)', self.onLabelsButtonClicked)
    self.LabelsButton = LabelsButton
    pelvisFloorLayout.addWidget(LabelsButton, row, column, rowHeight, colWidth)

    noOfColumns = 3 # number of columns in row
    colWidth = totalRowWidth / noOfColumns # width of single column

    row += 1 # row position
    column = 0 # column position

    ScaleButton = qt.QPushButton("Scale")
    ScaleButton.toolTip = "Scale model based on landmarks identifed as Markups Nodes. " + \
      "Landmarks defining APD and TD must be present"
    ScaleButton.connect('clicked(bool)', self.onScaleButtonClicked)
    self.ScaleButton = ScaleButton
    pelvisFloorLayout.addWidget(ScaleButton, row, column, rowHeight, colWidth)

    column += colWidth # column position

    MorphButton = qt.QPushButton("Morph")
    MorphButton.setStyleSheet('background-color: green')
    MorphButton.toolTip = "Morph model based on landmarks identifed as Markups Nodes. " + \
      "At least four landmarks not lying on a plane must be present"
    MorphButton.connect('clicked(bool)', self.onMorphButtonClicked)
    self.MorpButton = MorphButton
    pelvisFloorLayout.addWidget(MorphButton, row, column, rowHeight, colWidth)

    column += colWidth # column position

    BirthCanalButton = qt.QPushButton("Birth Canal")
    BirthCanalButton.toolTip = "View birth canal based on landmarks identifed as Markups Nodes."
    BirthCanalButton.connect('clicked(bool)', self.onBirthCanalButtonClicked)
    self.BirthCanalButton = BirthCanalButton
    pelvisFloorLayout.addWidget(BirthCanalButton, row, column, rowHeight, colWidth)

    noOfColumns = 2 # number of columns in row
    colWidth = totalRowWidth / noOfColumns # width of single column

    noOfColumns = 4 # number of columns in row
    colWidth = totalRowWidth / noOfColumns # width of single column

    row += 1 # row position
    column = 0 # column position

    DeleteButton = qt.QPushButton("Delete")
    DeleteButton.setStyleSheet('background-color: red')
    DeleteButton.toolTip = "Delete all models displayed in MRML scene as Model Nodes."
    DeleteButton.connect('clicked(bool)', self.onDeleteButtonClicked)
    self.DeleteButton = DeleteButton
    pelvisFloorLayout.addWidget(DeleteButton, row, column, rowHeight, colWidth)

    column += colWidth # column position

    DeleteMarkupsButton = qt.QPushButton("Landmarks")
    # DeleteMarkupsButton = qt.QPushButton("Markups")
    DeleteMarkupsButton.setStyleSheet('background-color: orange')
    DeleteMarkupsButton.toolTip = "Delete all models displayed in MRML scene as Model Nodes."
    DeleteMarkupsButton.connect('clicked(bool)', self.onDeleteMarkupsButtonClicked)
    self.DeleteMarkupsButton = DeleteMarkupsButton
    pelvisFloorLayout.addWidget(DeleteMarkupsButton, row, column, rowHeight, colWidth)

    column += colWidth # column position

    DeleteSegmentationsButton = qt.QPushButton("Segmentations")
    DeleteSegmentationsButton.setStyleSheet('background-color: orange')
    DeleteSegmentationsButton.toolTip = "Delete all models displayed in MRML scene as Model Nodes."
    DeleteSegmentationsButton.connect('clicked(bool)', self.onDeleteSegmentationsButtonClicked)
    self.DeleteSegmentationsButton = DeleteSegmentationsButton
    pelvisFloorLayout.addWidget(DeleteSegmentationsButton, row, column, rowHeight, colWidth)

    column += colWidth # column position

    DeleteModelsButton = qt.QPushButton("Models")
    DeleteModelsButton.setStyleSheet('background-color: orange')
    DeleteModelsButton.toolTip = "Delete all models displayed in MRML scene as Model Nodes."
    DeleteModelsButton.connect('clicked(bool)', self.onDeleteModelsButtonClicked)
    self.DeleteModelsButton = DeleteModelsButton
    pelvisFloorLayout.addWidget(DeleteModelsButton, row, column, rowHeight, colWidth)

    noOfColumns = 3 # number of columns in row
    colWidth = totalRowWidth / noOfColumns # width of single column

    row += 1 # row position
    column = 0 # column position

    UserScript1Button = qt.QPushButton("User Script 1")
    UserScript1Button.setStyleSheet('background-color: grey')
    UserScript1Button.connect('clicked(bool)', self.onUserScript1ButtonClicked)
    self.UserScript1Button = UserScript1Button
    pelvisFloorLayout.addWidget(UserScript1Button, row, column, rowHeight, colWidth)

    column += colWidth # column position

    UserScript2Button = qt.QPushButton("User Script 2")
    UserScript2Button.connect('clicked(bool)', self.onUserScript2ButtonClicked)
    self.UserScript2Button = UserScript2Button
    pelvisFloorLayout.addWidget(UserScript2Button, row, column, rowHeight, colWidth)

    column += colWidth # column position

    UserScript3Button = qt.QPushButton("User Script 3")
    UserScript3Button.connect('clicked(bool)', self.onUserScript3ButtonClicked)
    self.UserScript3Button = UserScript3Button
    pelvisFloorLayout.addWidget(UserScript3Button, row, column, rowHeight, colWidth)

    # ====================================================================================================
    # Collapsible button for 'Control'
    controlCollapsibleButton = ctk.ctkCollapsibleButton()
    controlCollapsibleButton.text = "Global Control"
    self.layout.addWidget(controlCollapsibleButton)
    self.layout.addStretch(1) # add vertical spacer, otherwise spread

    # ----------------------------------------------------------------------------------------------------
    # Layout within the 'Control' collapsible button
    ControlLayout = qt.QGridLayout(controlCollapsibleButton)

    rowHeight = 1 # height of single row
    totalRowWidth = 2 # total row width

    noOfColumns = 2 # number of columns in row
    colWidth = totalRowWidth / noOfColumns # width of single column

    row = 0 # row position
    column = 0 # column position

    RestartButton = qt.QPushButton("Restart")
    RestartButton.setStyleSheet('background-color: yellow')
    RestartButton.toolTip = "Restart Slicer."
    RestartButton.connect('clicked(bool)', self.onRestartButtonClicked)
    self.RestartButton = RestartButton
    ControlLayout.addWidget(RestartButton, row, column, rowHeight, colWidth)

    column += colWidth # column position

    ExitButton = qt.QPushButton("Exit")
    ExitButton.setStyleSheet('background-color: red')
    ExitButton.toolTip = "Exit Slicer."
    ExitButton.connect('clicked(bool)', self.onExitButtonClicked)
    self.ExitButton = ExitButton
    ControlLayout.addWidget(ExitButton, row, column, rowHeight, colWidth)

  # ====================================================================================================
  # Actions on click
  def onImagesButtonClicked(self): # read MRI
    PelvicFloorLib.imagesButton()
  
  def onOpenButtonClicked(self): # read template model
    PelvicFloorLib.openButton()

  def onTemplateButtonClicked(self): # read template model
    PelvicFloorLib.templateButton()

  def onRotLButtonClicked(self): # rotate model around L
    PelvicFloorLib.rotateButton('L')

  def onRotPButtonClicked(self): # rotate model around P
    PelvicFloorLib.rotateButton('P')

  def onRotSButtonClicked(self): # rotate model around S
    PelvicFloorLib.rotateButton('S')

  def onLandmarksButtonClicked(self): # read landmarks
    PelvicFloorLib.landmarksButton()

  def onAlignButtonClicked(self): # align model
    PelvicFloorLib.alignButton('PSI')

  def onSegmentationButtonClicked(self): # segment model
    PelvicFloorLib.segmentationButton()

  def onLabelsButtonClicked(self): # label model
    PelvicFloorLib.labelsButton()

  def onScaleButtonClicked(self): # scale model
    PelvicFloorLib.scaleButton()

  def onMorphButtonClicked(self): # morph model
    PelvicFloorLib.morphButton()

  def onBirthCanalButtonClicked(self): # birth canal view
    PelvicFloorLib.canalButton()

  def onExportButtonClicked(self): # birth canal view
    PelvicFloorLib.exportButton()

  def onDeleteButtonClicked(self): # delete models
    PelvicFloorLib.deleteButton()
    
  def onDeleteMarkupsButtonClicked(self): # delete models
    PelvicFloorLib.deleteMarkupsButton()
    
  def onDeleteSegmentationsButtonClicked(self): # delete models
    PelvicFloorLib.deleteSegmentationsButton()
    
  def onDeleteModelsButtonClicked(self): # delete models
    PelvicFloorLib.deleteModelsButton()
    
  def onUserScript1ButtonClicked(self): # delete models
    PelvicFloorLib.userScript1Button()
    
  def onUserScript2ButtonClicked(self): # delete models
    PelvicFloorLib.userScript2Button()
    
  def onUserScript3ButtonClicked(self): # delete models
    PelvicFloorLib.userScript3Button()
    
  def onRestartButtonClicked(self): # restart Slicer
    slicer.util.restart()

  def onExitButtonClicked(self): # exit Slicer
    slicer.util.exit()
