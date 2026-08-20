# ====================================================================================================
from __main__ import qt
from PelvicFloorLib import QMessageTitle, pathName, fileName, extension, \
    read, update, getLandmarks, \
        deleteModels, deleteMarkups, \
            registrationButton, exportButton
from scipy.interpolate import splprep, splev
import slicer, numpy as np
from vtk.util import numpy_support

# ====================================================================================================
# user script button 1
def userScript1Button():
    """
    Script defined by user.
    
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
    """
    Script defined by user.
    
    Returns
    -------
    Defined by user.
    """

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
    Script defined by user.
    
    Returns
    -------
    Defined by user.
    """
