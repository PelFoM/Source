shortCuts = [
  ("Ctrl+b", lambda: slicer.app.layoutManager().setLayout(slicer.vtkMRMLLayoutNode.SlicerLayoutOneUpRedSliceView)),
  ("Ctrl+n", lambda: slicer.app.layoutManager().setLayout(slicer.vtkMRMLLayoutNode.SlicerLayoutOneUpYellowSliceView)),
  ("Ctrl+m", lambda: slicer.app.layoutManager().setLayout(slicer.vtkMRMLLayoutNode.SlicerLayoutOneUpGreenSliceView)),
  ("Ctrl+,", lambda: slicer.app.layoutManager().setLayout(slicer.vtkMRMLLayoutNode.SlicerLayoutFourUpView)),
  ("Ctrl+e", lambda: slicer.util.restart())
  ]

for (shortCutKey, callBack) in shortCuts:
  shortCut = qt.QShortcut(slicer.util.mainWindow())
  shortCut.setKey(qt.QKeySequence(shortCutKey))
  shortCut.connect( "activated()", callBack)

import numpy as np
from vtk.util import numpy_support
