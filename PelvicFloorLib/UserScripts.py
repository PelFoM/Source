# ====================================================================================================
from __main__ import qt
from PelvicFloorLib import read, align, update, move, \
    morphButton, exportButton, \
        deleteModels, deleteMarkups, \
            pathName, fileName, extension
import json, slicer, vtk, numpy as np, time
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
    subjects = ['P041', 'P071', 'P251', 'P292', 'P381', 'P401']
    # subjects = ['P041']

    # subjects = ['P011', 'P031', 'P041', 'P061', 'P071', 'P081', 'P091', 
    #             'P101', 'P111', 'P121', 'P131', 'P141', 'P161', 'P171', 'P181', 
    #             'P231', 'P232', 'P241', 'P251', 'P252', 'P261', 'P271', 'P272', 'P281', 'P282', 'P283', 'P291', 'P292', 'P293', 
    #             'P301', 'P302', 'P311', 'P321', 'P322', 'P323', 'P331', 'P332', 'P333', 'P341', 'P361', 'P371', 'P381', 'P391', 
    #             'P401', 'P411', 'P421', 'P431', 'P441', 'P451', 
    #             'P511', 'P521', 'P531', 'P541', 'P551', 'P561', 'P571', 'P581', 'P591', 
    #             'P601', 'P611', 'P621', 'P631', 'P641', 'P651']

    for subject in subjects:
        deleteModels() # delete existing models and read template model
        modelNode = read(pathName + fileName + extension, 'Template model', 'vtu')
        points = numpy_support.vtk_to_numpy(modelNode.GetMesh().GetPoints().GetData())
        points[:, [0, 1, 2]] = -points[:, [2, 0, 1]] # convert to LPS (internally RAS) coordinate system
        update(modelNode) # show model

        # read landmarks
        with open(pathName + '../TargetModels/' + subject + '_target.txt') as landmarks:
            landmarksDictXyz = json.loads(landmarks.read())
        for landmark in landmarksDictXyz: # separate markups
            node = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLMarkupsFiducialNode')
            node.SetName(landmark)
            slicer.util.updateMarkupsControlPointsFromArray(node, np.array([landmarksDictXyz[landmark]]))

        # align to PSI and remove PSI
        noOfModelNode = slicer.mrmlScene.GetNumberOfNodesByClass('vtkMRMLModelNode')
        align(slicer.mrmlScene.GetNthNodeByClass(noOfModelNode - 1, 'vtkMRMLModelNode'), 'PSI')
        node = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLMarkupsFiducialNode')
        slicer.mrmlScene.RemoveNode(slicer.mrmlScene.GetNthNodeByClass(0, 'vtkMRMLMarkupsFiducialNode'))

        # morphing and export
        morphButton()
        exportButton()
        # noOfModelNode = slicer.mrmlScene.GetNumberOfNodesByClass('vtkMRMLModelNode')
        # modelNode = slicer.mrmlScene.GetNthNodeByClass(noOfModelNode - 1, 'vtkMRMLModelNode')
        # modelWriter = vtk.vtkPolyDataWriter()
        # modelWriter.SetFileName(pathName + '../TargetModels/' + subject + '_morphed.vtk')
        # modelWriter.SetInputData(modelNode.GetMesh())
        # modelWriter.Write()
        
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

    deleteModels() # delete existing models and read template model
    modelNodePelvis = read(pathName + fileName + extension, 'Template model', 'vtu')
    points = numpy_support.vtk_to_numpy(modelNodePelvis.GetMesh().GetPoints().GetData())
    points[:, [0, 1, 2]] = -points[:, [2, 0, 1]] # convert to LPS (internally RAS) coordinate system
    update(modelNodePelvis) # show model

    modelNodeHead = read(pathName + 'head.vtk', 'Template model', 'vtu')
    points = numpy_support.vtk_to_numpy(modelNodeHead.GetMesh().GetPoints().GetData())
    points[:, [0, 1, 2]] = -points[:, [2, 0, 1]] # convert to LPS (internally RAS) coordinate system
    modelNodeHead = move(modelNodeHead, [100, -300, -150]) # show model

    # slicer.mrmlScene.GetSubjectHierarchyNode().GetNumberOfItems()
    # list(slicer.mrmlScene.GetNodesByClass('vtkMRMLDisplayableNode'))

    # https://github.com/Slicer/Slicer/blob/main/Modules/Scripted/ScreenCapture/ScreenCapture.py
    for t in range(10):
        # modelNodeMove = read(pathName + 'head.vtk', 'Time = ' + str(t), 'vtu')
        # move(modelNodeMove, [100, -300, -150 - t])
        # continue

        # slicer.mrmlScene.RemoveNode(modelNodeHead)
        # modelNodeHead = read(pathName + 'head.vtk', 'Template model', 'vtu')
        # points = numpy_support.vtk_to_numpy(modelNodeHead.GetMesh().GetPoints().GetData())
        # points[:, [0, 1, 2]] = -points[:, [2, 0, 1]] # convert to LPS (internally RAS) coordinate system
        # modelNodeHead = move(modelNodeHead, [100, -300, -150 - 3 * t]) # show model
    
        points = numpy_support.vtk_to_numpy(modelNodeHead.GetMesh().GetPoints().GetData())
        # points = points + vector # creates local variable points
        points += [0, 0, - 3 * t] # updated points by pointer
        modelNodeHead.GetDisplayNode().SetVisibility(1)

        # volNode.GetDisplayNode().SetVisibility(0)
        # volNode.SetDisplayVisibility(0)
        
        # slicer.util.delayDisplay("Head downloaded...",10)

        # volumeNode.AddAndObserveDisplayNodeID(displayNode.GetID())
        # logic.UpdateDisplayNodeFromVolumeNode(displayNode, volumeNode)

        # slicer.vtkMRMLDisplayNode("vtkMRMLSegmentationNode")
        # modelNodeHead.CreateDefaultDisplayNodes()
        # modelNodeHead.DisplayNodes()
        # modelNodeHead = update(modelNodeHead)
        # modelNodeHead.Update()
        # time.sleep(1)

        # # get model mesh and remove
        # volMesh = modelNodeMove.GetMesh()
        # slicer.mrmlScene.RemoveNode(modelNodeMove)

        # # show new model mesh
        # extractSurface = vtk.vtkGeometryFilter()
        # extractSurface.SetInputData(volMesh)
        # extractSurface.Update()
        # modelNodeMove = slicer.modules.models.logic().AddModel(extractSurface.GetOutput())
        # modelNodeMove.GetDisplayNode().SetOpacity(opacity)

    # noOfModelNode = slicer.mrmlScene.GetNumberOfNodesByClass('vtkMRMLModelNode')
    # align(slicer.mrmlScene.GetNthNodeByClass(noOfModelNode - 1, 'vtkMRMLModelNode'), 'PSI')
    # node = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLMarkupsFiducialNode')
    # slicer.mrmlScene.RemoveNode(slicer.mrmlScene.GetNthNodeByClass(0, 'vtkMRMLMarkupsFiducialNode'))

    # morphing and export
    # morphButton()   
    # noOfModelNode = slicer.mrmlScene.GetNumberOfNodesByClass('vtkMRMLModelNode')
    # modelNode = slicer.mrmlScene.GetNthNodeByClass(noOfModelNode - 1, 'vtkMRMLModelNode')
        
    # pass
