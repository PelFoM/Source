# ====================================================================================================
from __main__ import qt
from PelvicFloorLib import pathName, fileName, extension, \
    read, update, \
        deleteModels, deleteMarkups, \
            registrationButton, exportButton
import slicer, numpy as np
from vtk.util import numpy_support

# ====================================================================================================
# user script button 1
def userScript1Button():
    """
    Scrip defined by user.
    
    Returns
    -------
    Defined by user.
    """

    # subjects = ['P041', 'P071', 'P251', 'P292', 'P341', 'P381', 'P401']
    # subjects = ['P041', 'P071', 'P251', 'P292', 'P381', 'P401']
    # subjects = ['P041']

    # subjects = ['P011', 'P031', 'P041', 'P061', 'P071', 'P081', 'P091', 
    #             'P101', 'P111', 'P121', 'P131', 'P141', 'P161', 'P171', 'P181', 
    #             'P231', 'P232', 'P241', 'P251', 'P252', 'P261', 'P271', 'P272', 'P281', 'P282', 'P283', 'P291', 'P292', 'P293', 
    #             'P301', 'P302', 'P311', 'P321', 'P322', 'P323', 'P331', 'P332', 'P333', 'P341', 'P361', 'P371', 'P381', 'P391', 
    #             'P401', 'P411', 'P421', 'P431', 'P441', 'P451', 
    #             'P511', 'P521', 'P531', 'P541', 'P551', 'P561', 'P571', 'P581', 'P591', 
    #             'P601', 'P611', 'P621', 'P631', 'P641', 'P651']

    subjects = ['P041', 'P381', 'P401', 'PS00'] # 'P071'
    landmarks = ['bone', 'palp']
    #methods = ['TPS']
    methods = ['WC2']
    path = '../../output/presentations/2026_CBM/Morphing/'

    for subject in subjects:

        for landmark in landmarks:

            for method in methods:

                deleteModels() # delete existing models and read template model
                modelNode = read(pathName + fileName + extension, subject, 'vtk')
                points = numpy_support.vtk_to_numpy(modelNode.GetMesh().GetPoints().GetData())
                points[:, [0, 1]] *= -1 # convert to LPS (internally RAS) coordinate system
                update(modelNode) # show model

                with open(pathName + path + subject + '_landmarks_' + landmark + '.fcsv') as lines:

                    for line in lines:
                        if line.startswith('#'):
                            continue

                        currentLine = line.strip().split(',')
                        
                        # LPS -> RAS
                        x = -float(currentLine[1]) 
                        y = -float(currentLine[2])
                        z =  float(currentLine[3])

                        node = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLMarkupsFiducialNode")
                        node.SetName(currentLine[11])
                        node.AddControlPoint(x, y, z)

                        display = node.GetDisplayNode()
                        display.SetSelectedColor(1, 0, 0) # displayed in red
                        display.SetColor(0, 1, 0) # green when selected

                # morphing and export
                registrationButton('rigid')
                registrationButton('morphing')
                exportButton(pathName + path + subject + '_morphed_' + landmark + '_' + method + '.vtk')
                
                # delete all
                deleteModels()
                deleteMarkups()

# ====================================================================================================
# user script button 2
def userScript2Button():
    # import vtkSlicerSegmentComparisonModuleLogicPython
    # s=slicer.util.getNode('vtkMRMLSegmentationNode1') # Whatever your input is
    # p1=s.GetClosedSurfaceRepresentation('Segment_1') # Again, depends on your input. If it's not a segmentation then you'll need to access the model nodes
    # p2=s.GetClosedSurfaceRepresentation('Segment_2')
    # pdf=vtkSlicerSegmentComparisonModuleLogicPython.vtkPolyDataDistanceHistogramFilter()
    # pdf.SetInputReferencePolyData(p1)
    # pdf.SetInputComparePolyData(p2)
    # pdf.Update()
    # pdf.GetNthPercentileHausdorffDistance(0)
    
    # noOfSegmentationNode = slicer.mrmlScene.GetNumberOfNodesByClass('vtkMRMLSegmentationNode')
    # -> must be 2
    segmentationNode1 = slicer.mrmlScene.GetNthNodeByClass(0, 'vtkMRMLSegmentationNode')
    segment_1 = segmentationNode1.GetClosedSurfaceInternalRepresentation('Model_3')

    segmentationNode2 = slicer.mrmlScene.GetNthNodeByClass(1, 'vtkMRMLSegmentationNode')
    segment_2 = segmentationNode2.GetClosedSurfaceInternalRepresentation('Segment_1')

    # np.size(numpy_support.vtk_to_numpy(segment_1.GetPoints().GetData()))
    X = numpy_support.vtk_to_numpy(segment_1.GetPoints().GetData())
    # os.getcwd()
    print(fileName)
    np.savetxt('matrix_1', X, fmt='%.2e')
    
    # np.size(numpy_support.vtk_to_numpy(segment_2.GetPoints().GetData()))
    Y = numpy_support.vtk_to_numpy(segment_2.GetPoints().GetData())
    np.savetxt('matrix_2.txt', Y, fmt='%.2e')
    
    # print(segmentationNode1.GetName())
    # print(segment_1)

    # import vtkSlicerSegmentComparisonModuleLogicPython
    # pdf = vtkSlicerSegmentComparisonModuleLogicPython.vtkPolyDataDistanceHistogramFilter()

    # pdf.SetInputReferencePolyData(segment_1)
    # pdf.SetInputComparePolyData(segment_2)
    # pdf.Update()
    # pdf.GetNthPercentileHausdorffDistance(0)

# ====================================================================================================
# user script button 3
def userScript3Button():
    """
    Scrip defined by user.
    
    Returns
    -------
    Defined by user.
    """

    import vtk
    #from scipy.interpolate import interp1d
    from scipy.interpolate import CubicSpline

    # ellipse from 4 points
    def ellipse_from_diameters(p1, p2, p3, p4, numberOfPoints=100):
        center = (p1 + p2 + p3 + p4) / 4
        axis_a = (p1 - p2) / 2
        axis_b = (p3 - p4) / 2
        theta = np.linspace(0, 2 * np.pi, numberOfPoints)
        ellipse = np.array([center + np.cos(t) * axis_a + np.sin(t) * axis_b for t in theta])
        return ellipse

    ellipses_points = [
        (np.array([1, 0, 0]), np.array([-1, 0, 0]), np.array([0, 0.5, 0]), np.array([0, -0.5, 0])), # pelvic inlet
        (np.array([0.8, 0.2, 1]), np.array([-1.2, 0.3, 1]), np.array([0.1, 0.7, 1]), np.array([-0.1, -0.7, 1])), # middle plane
        (np.array([0.6, 0.5, 2]), np.array([-0.9, 0.6, 2]), np.array([0.2, 0.9, 2]), np.array([-0.2, -0.9, 2])) # pelvic outlet
    ]

    ellipses = [ellipse_from_diameters(*pts) for pts in ellipses_points]
    ellipses = np.array(ellipses)
    Z = [ellipse[0, 2] for ellipse in ellipses]

    # cubic spline between ellipses
    numberOfLayers = 60
    Zi = np.linspace(min(Z), max(Z), numberOfLayers)
    numberOfPoints = ellipses.shape[1]

    ellipses_smooth = []
    for point in range(numberOfPoints):
        x = ellipses[:, point, 0]
        y = ellipses[:, point, 1]
        #z = ellipses[:, point, 2]
        #xi = interp1d(Z, x)
        xi = CubicSpline(Z, x)
        #yi = interp1d(Z, y)
        yi = CubicSpline(Z, y)
        Xi = xi(Zi)
        Yi = yi(Zi)
        ellipses_smooth.append(np.vstack((Xi, Yi, Zi)).T)

    surface = np.array(ellipses_smooth).transpose(1, 0, 2)  # (layers, points, xyz)

    # vtk polygon mesh
    points = vtk.vtkPoints()
    polys = vtk.vtkCellArray()

    numberOfLayers, numOfPoints, _ = surface.shape
    for layer in range(numberOfLayers):
        for point in range(numOfPoints):
            points.InsertNextPoint(surface[layer, point])

    # add triangles
    for layer in range(numberOfLayers - 1):
        for point in range(numOfPoints - 1):
            i0 = layer * numOfPoints + point
            i1 = i0 + 1
            i2 = i0 + numOfPoints
            i3 = i2 + 1

            # 2 triangle cells
            polys.InsertNextCell(3)
            polys.InsertCellPoint(i0)
            polys.InsertCellPoint(i2)
            polys.InsertCellPoint(i1)
            polys.InsertNextCell(3)
            polys.InsertCellPoint(i1)
            polys.InsertCellPoint(i2)
            polys.InsertCellPoint(i3)

    # closing by last point in circle
    for layer in range(numberOfLayers - 1):
        i0 = layer * numOfPoints + (numOfPoints - 1)
        i1 = layer * numOfPoints
        i2 = i0 + numOfPoints
        i3 = i2 - (numOfPoints - 1)
        polys.InsertNextCell(3)
        polys.InsertCellPoint(i0)
        polys.InsertCellPoint(i2)
        polys.InsertCellPoint(i1)
        polys.InsertNextCell(3)
        polys.InsertCellPoint(i1)
        polys.InsertCellPoint(i2)
        polys.InsertCellPoint(i3)

    polydata = vtk.vtkPolyData()
    polydata.SetPoints(points)
    polydata.SetPolys(polys)

    # add loft to MRML scene
    modelNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLModelNode", "Smooth Loft Surface")
    modelNode.SetAndObservePolyData(polydata)

    displayNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLModelDisplayNode")
    displayNode.SetColor(0.3, 0.7, 1.0)
    displayNode.SetOpacity(0.7)
    displayNode.SetBackfaceCulling(False)
    modelNode.SetAndObserveDisplayNodeID(displayNode.GetID())

    slicer.app.layoutManager().threeDWidget(0).threeDView().resetFocalPoint()

    # add ellipses to MRML scene
    for idx, ellipse in enumerate(ellipses):
        pointsEllipse = vtk.vtkPoints()
        for point in ellipse:
            pointsEllipse.InsertNextPoint(point)

        lines = vtk.vtkCellArray()
        polyline = vtk.vtkPolyLine()
        polyline.GetPointIds().SetNumberOfIds(len(ellipse))
        for i in range(len(ellipse)):
            polyline.GetPointIds().SetId(i, i)
        lines.InsertNextCell(polyline)

        polydataEllipse = vtk.vtkPolyData()
        polydataEllipse.SetPoints(pointsEllipse)
        polydataEllipse.SetLines(lines)

        # add ellise to MRML scene
        modelNodeEllipse = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLModelNode", f"Ellipse {idx+1}")
        modelNodeEllipse.SetAndObservePolyData(polydataEllipse)

        displayNodeEllipse = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLModelDisplayNode")
        displayNodeEllipse.SetColor(1, 0, 0)
        displayNodeEllipse.SetLineWidth(2)
        displayNodeEllipse.SetOpacity(0.5)
        displayNodeEllipse.SetBackfaceCulling(False)
        modelNodeEllipse.SetAndObserveDisplayNodeID(displayNodeEllipse.GetID())
