# ====================================================================================================
from __main__ import qt
from PelvicFloorLib import pathName, fileName, QMessageTitle, opacity
import json, slicer, vtk, numpy as np
from vtk.util import numpy_support

# ====================================================================================================
# read model
def read(pathFileName, modelName, reader):
    """
    Read model from VTK (polydata or unstructured grid)

    Returns
    -------
    The model read in the main MRML screen.
    """

    # define MRML window
    view = slicer.app.layoutManager().threeDWidget(0).threeDView()
    #view = slicer.app.layoutManager().sliceWidget("Red").sliceView()
    view.cornerAnnotation().SetText(vtk.vtkCornerAnnotation.UpperRight, modelName)
    view.cornerAnnotation().GetTextProperty().SetColor(1, 1, 1)
    view.forceRender() # update view

    # position = [0, 0, 0]
    # crosshairNode = slicer.util.getNode("Crosshair")
    # crosshairNode.SetCrosshairRAS(position) # set crosshair position
    # slicer.vtkMRMLSliceNode.JumpAllSlices(slicer.mrmlScene, *position, # center position in all slice views
    #                                       slicer.vtkMRMLSliceNode.CenteredJumpSlice)
    # # make crosshair visible
    # crosshairNode.SetCrosshairMode(slicer.vtkMRMLCrosshairNode.ShowBasic)

    # volMesh = vtk.vtkUnstructuredGrid()
    if reader == 'vtu':
        modelReader = vtk.vtkUnstructuredGridReader()
    else: # 'vtp'
        modelReader = vtk.vtkPolyDataReader()
        
    modelReader.SetFileName(pathFileName)
    modelReader.Update()

    # volMesh = modelReader.GetOutput()
    # volMesh.GetPoints()
    # volMesh.GetPoints().GetData()
    # numpy_support.vtk_to_numpy(volMesh.GetPoints().GetData())
    # volMesh.GetPolys()
    # volMesh.GetPolys().GetData()
    # numpy_support.vtk_to_numpy(volMesh.GetPolys().GetData())
    # volMesh.GetCellData()
    # volMesh.GetCellData().GetArray(0)
    # numpy_support.vtk_to_numpy(volMesh.GetCellData().GetArray(0))
    
    # make model active (necessary to access data)
    slicer.modules.models.logic().AddModel(modelReader.GetOutputPort())

    # number of model nodes
    noOfModelNode = slicer.mrmlScene.GetNumberOfNodesByClass('vtkMRMLModelNode')
    
    # last model node (last read)
    modelNode = slicer.mrmlScene.GetNthNodeByClass(noOfModelNode - 1, 'vtkMRMLModelNode')
    
    return(modelNode)

# ====================================================================================================
# export model
# def export(pathFileName, modelName):
#     """
#     Export model to PC.
#
#     Returns
#     -------
#     The model exported.
#     """

# ====================================================================================================
# write tetrahedral mesh mesh to PC file
def export(outputFileName, volMesh):
    """
    Write segmented tetrahedral mesh to 'fileName'.

    Returns
    -------
    The file save to 'fileName' in the PC format.
    """     
    
    fileOut = open(outputFileName + '.pc', 'w')

    nodeCoordinates = numpy_support.vtk_to_numpy(volMesh.GetPoints().GetData())
    for index in range(len(nodeCoordinates)): # node numbering start from 1
        fileOut.write("NODE  / %8d%16.3f%16.3f%16.3f\n" % (index + 1, 
                                                           nodeCoordinates[index][0], 
                                                           nodeCoordinates[index][1], 
                                                           nodeCoordinates[index][2]))
    triangleNumber = 0
    triangles = numpy_support.vtk_to_numpy(volMesh.GetPolys().GetData())
    triangleNodes = triangles.reshape(int(len(triangles) / 4), 4)
    for index in range(len(triangleNodes)):
        triangleNumber += 1 # triangle numbering start from 1
        fileOut.write("SHELL / %8d%8d%8d%8d%8d\n" % (triangleNumber, 
                                                     triangleNodes[index][0], 
                                                     triangleNodes[index][1] + 1, 
                                                     triangleNodes[index][2] + 1, 
                                                     triangleNodes[index][3] + 1))
        
    fileOut.close()

# ====================================================================================================
# create mesh from polydata
# https://discourse.slicer.org/t/create-a-model-mesh-in-python-from-array-of-vertices-and-triangles/5541
def createMesh(arrayVertices, arrayVertexNormals, arrayTriangles, labelsScalars, arrayScalars):
    """
    Create mesh from elementary VTK data.

    modelNode : a vtkMRMLModelNode in the Slicer scene which will hold the mesh
    arrayVertices : list of triples [[x1,y1,z2], [x2,y2,z2], ... ,[xn,yn,zn]] of vertex coordinates
    arrayVertexNormals : list of triples [[nx1,ny1,nz2], [nx2,ny2,nz2], ... ] of vertex normals
    arrayTriangles : list of triples of 0-based indices defining triangles
    labelsScalars : list of strings such as ["bipolar", "unipolar"] to label the individual scalars data sets
    arrayScalars : an array of n rows for n vertices and m colums for m inidividual scalar sets

    Returns
    -------
    mesh : vtkmodules.vtkCommonDataModel.vtkPolyData
        Polydata mesh.
    """

    def mkVtkIdList(it):
        vil = vtk.vtkIdList()
        for i in it:
            vil.InsertNextId(int(i))
        return vil

    # create the building blocks of polydata including data attributes
    mesh = vtk.vtkPolyData()
    points = vtk.vtkPoints()
    normals = vtk.vtkFloatArray()
    polys = vtk.vtkCellArray()
    
    # load the array data into the respective VTK data structures
    for i in range(len(arrayVertices)):
        points.InsertPoint(i, arrayVertices[i])
    
    for i in range(len(arrayTriangles)):
        polys.InsertNextCell(mkVtkIdList(arrayTriangles[i]))
    
    for i in range(len(arrayVertexNormals)):
        normals.InsertTuple3(i, arrayVertexNormals[i][0], arrayVertexNormals[i][1], arrayVertexNormals[i][2])
    
    # put together the mesh object
    mesh.SetPoints(points)
    mesh.SetPolys(polys)
    if(len(arrayVertexNormals) == len(arrayVertices)):
        mesh.GetPointData().SetNormals(normals)
    
    # add scalars
    scalars = []
    for j in range(len(labelsScalars)):
        scalars.append(vtk.vtkFloatArray())
        
        for i in range(len(arrayVertices)):
            scalars[j].InsertTuple1(i,arrayScalars[i][j])
        
        scalars[j].SetName(labelsScalars[j])
        mesh.GetPointData().AddArray(scalars[j])
    
    return(mesh)

# ====================================================================================================
# separate model by labels
def separate(volMesh):
    """
    Separate model mesh by labels.

    Returns
    -------
    segments : list
        The list of segment labels.
    segmentVertices : array
        The dictionary of the set of vertices by label.
    segmentTriangles : array
        The dictionary of the set of triangles by label.
    """

    # volMesh.GetCellData() -> vtkCellData
    # volMesh.GetCellData().GetArray(0) -> vtkIntArray
    triangleSegments = numpy_support.vtk_to_numpy(volMesh.GetCellData().GetArray(0))
    segments = list(set(triangleSegments)) # unique elements
    
    # volMesh.GetPoints() -> vtkPoints
    # volMesh.GetPoints().GetData() -> tkTypeFloat64Array
    vertices = numpy_support.vtk_to_numpy(volMesh.GetPoints().GetData())

    # volMesh.GetPolys() -> vtkCellArray
    # volMesh.GetPolys().GetData() -> vtkIdTypeArray
    # volMesh.GetPolys().GetNumberOfCells() = size(vtkIdTypeArray) / 4
    triangles = numpy_support.vtk_to_numpy(volMesh.GetPolys().GetData())
    triangles = triangles.reshape(int(len(triangles) / 4), 4)
    
    # triangleNormals = [] # triangle normals    
    # for triangle in range(len(triangles)):
    #     # type = triangles[triangle][0] # unstructured grid type (3 for triangle)
    #     X1 = vertices[triangles[triangle][1]] # triangle vertex 1 coordinates
    #     X2 = vertices[triangles[triangle][2]] # triangle vertex 2 coordinates
    #     X3 = vertices[triangles[triangle][3]] # triangle vertex 3 coordinates
    #     triangleNormals.append(np.cross(X2 - X1, X3 - X1))

    # vertexNormals = [] # vertex normals
    # for vertex in range(len(vertices)):
    #     vertexNormals.append(np.average(triangleNormals[triangle] for triangle in range(len(triangles))))

    # separate segments
    segmentVertices = {} # dictionary of vertices per segment
    segmentTriangles = {} # dictionary of triangles per segment
    
    for segment in segments: # initialize dictionaries
        segmentTriangles[segment] = triangles[triangleSegments == segment][:, 1 : ]
        vertexNumbers = np.unique(segmentTriangles[segment]) # triangle vertices numbers

        # renumber
        segmentVertices[segment] = vertices#[vertexNumbers] # triangle vertices unique numbers
        # ren = renumber(vertexNumbers, np.arange(len(vertexNumbers)))
        # segmentTriangles[segment] = ren[segmentTriangles[segment]]

    return(segments, segmentVertices, segmentTriangles)

# ====================================================================================================
def renumber(oldIndices, newIndices):
    """
    Prepare vector for renumbering from 'oldIndices' to 'newIndices'.

    Returns
    -------
    renumber : array
        The renumbering vector. Then 'newArgument = renumber[oldArrgument]'.
    """

    oldIndices = np.asarray(oldIndices)
    newIndices = np.asarray(newIndices)

    renumber = np.zeros(oldIndices.max() + 1, dtype=newIndices.dtype)
    renumber[oldIndices] = newIndices

    return renumber

# ====================================================================================================
# get coincident landmarks from model and markups node
def getLandmarks(modelNode, landmarksDict):
    """
    Identify existence of image landmarks from 'landmarksDict' in model landmarks.
    Model landmarks file '*_landmarks.dic' must be in model numbering.

    Returns
    -------
    allExist : logical value
        The existence of all landmarks.
    modelDictXyz : dictionary
        The dictionary of model landmarks (key = landmark, value = array).
    imageDictXyz : dictionary
        The dictionary of MRI landmarks (key = landmark, value = array).
    """

    # MRI landmarks coordinates
    imageDictXyz = {} # dictionary of MRI landmarks coordinates
    noOfMarkupsNode = slicer.mrmlScene.GetNumberOfNodesByClass('vtkMRMLMarkupsNode')

    # model landmark coordinates in local order numbering
    with open(pathName + fileName + '_landmarks.dic') as landmarks:
        modelDictXyzAll = json.loads(landmarks.read())

    modelDictXyz = {} # dictionary of model landmarks coordinates
    points = numpy_support.vtk_to_numpy(modelNode.GetMesh().GetPoints().GetData())    
    
    # markupsNode remain in the database after deleting including all properties,
    # therefore reversed loop search from the newest one
    for i in reversed(range(noOfMarkupsNode)): # from the newest to find the newest

        markupsNode = slicer.mrmlScene.GetNthNodeByClass(i, 'vtkMRMLMarkupsNode')
        if markupsNode.GetName() in landmarksDict.keys():
            landmark = markupsNode.GetName()
            landmarksDict[landmark] = True

            # add MRI and model landmarks
            imageDictXyz.update({landmark: np.array(markupsNode.GetNthControlPointPosition(0))})
            modelDictXyz.update({landmark: points[modelDictXyzAll[landmark]]})

    # all necessary landmarks exist
    if all(value == True for value in landmarksDict.values()):
        allExist = True           
    else:
        allExist = False
        for landmark in landmarksDict.keys():
            if not landmarksDict[landmark]:
                qt.QMessageBox.information(slicer.util.mainWindow(), 
                                           QMessageTitle, landmark + \
                                           ' does not exist.')

    return(allExist, modelDictXyz, imageDictXyz)

# ====================================================================================================
# align model
def align(modelNode, landmark):
    """
    Align model to 'landmark'.
    
    Returns
    -------
    The aligned model displayed in the main MRML screen.
    If 'landmark' does not exist, nothing happens.
    landmarkExists : True / False
        The logical variable if landmark exists.
    """

    # # get current model mesh points
    # points = numpy_support.vtk_to_numpy(modelNode.GetMesh().GetPoints().GetData())

    # check existence of landmark
    landmarksDict = {landmark: False}
    landmarkExists, modelDictXyz, imageDictXyz = getLandmarks(modelNode, landmarksDict)

    if landmarkExists: # align to landmark
        vector = imageDictXyz[landmark] - modelDictXyz[landmark] # direction vector
        move(modelNode, vector) # points are updated by memore inside move
    else:
        pass

    return(landmarkExists)

# ====================================================================================================
# move model
def move(modelNode, vector):
    """
    Move model node by 'vector'.
    
    Returns
    -------
    The model moved by 'vector' and displayed in the main MRML screen.
    """

    points = numpy_support.vtk_to_numpy(modelNode.GetMesh().GetPoints().GetData())
    #points = points + vector # creates local variable points
    points += vector # updated points by pointer

    # update model node
    movedModelNode = update(modelNode)
    #return(points)
    return(movedModelNode)

# ====================================================================================================
# rotate model around vector by angle
def rotate(modelNode, landmark, vector, theta):
    """
    Rotate model node in 'landmark' around 'vector' by 'angle' in radians.
    
    Returns
    -------
    The rotated model displayed the main MRML screen.
    """

    # rotatin unit vector
    urn = np.linalg.norm(vector)
    urx = vector[0] / urn
    ury = vector[1] / urn
    urz = vector[2] / urn

    # rotatin matrix (https://en.wikipedia.org/wiki/Rotation_matrix)
    rotationMatrix = [[np.cos(theta) + (urx ** 2) * (1.0 - np.cos(theta)),
                       urx * ury * (1.0 - np.cos(theta)) - urz * np.sin(theta),
                       urx * urz * (1.0 - np.cos(theta)) + ury * np.sin(theta)], 
                      [ury * urx * (1.0 - np.cos(theta)) + urz * np.sin(theta), 
                       np.cos(theta) + (ury ** 2) * (1.0 - np.cos(theta)), 
                       ury * urz * (1.0 - np.cos(theta)) - urx * np.sin(theta)], 
                      [urz * urx * (1.0 - np.cos(theta)) - ury * np.sin(theta), 
                       urz * ury * (1.0 - np.cos(theta)) + urx * np.sin(theta), 
                       np.cos(theta) + (urz ** 2) * (1.0 - np.cos(theta))]]

    # get current model mesh points
    points = numpy_support.vtk_to_numpy(modelNode.GetMesh().GetPoints().GetData())

    # check existence of landmark
    landmarksDict = {landmark: False}
    landmarkExists, modelDictXyz, imageDictXyz = getLandmarks(modelNode, landmarksDict)

    if landmarkExists: # rotate around vector going through landmark

        # (R * M')' ~ M * R'
        #points = np.dot(points, np.transpose(rotationMatrix))
        points = np.matmul(points - modelDictXyz[landmark], \
                           np.transpose(rotationMatrix)) + modelDictXyz[landmark]

    else: # rotate around vector going through (0, 0, 0)
        points = np.dot(points, np.transpose(rotationMatrix))
    
    # modelNode.GetMesh().GetPoints() is vtkmodules.vtkCommonCore.vtkPoints
    newPoints = vtk.vtkPoints()
    for i in range(len(points)):
        newPoints.InsertPoint(i, points[i])
    
    # update points
    # newPoints.InsertPoints(numpy_support.numpy_to_vtk(points))
    modelNode.GetMesh().SetPoints(newPoints)

    # update model node
    rotatedModelNode = update(modelNode)
    return(rotatedModelNode)

# ====================================================================================================
# scale model
def scale(modelNode, vector):
    """
    Scale model node by 'vector'.
    
    Returns
    -------
    The model scaled by 'vector' and displayed in the main MRML screen.
    """
    
    points = numpy_support.vtk_to_numpy(modelNode.GetMesh().GetPoints().GetData())
    points[:, [0]] *= vector[0]
    points[:, [1]] *= vector[1]
    points[:, [2]] *= vector[2]

    # update model node
    scaledModelNode = update(modelNode)
    #return(points)
    return(scaledModelNode)

# ====================================================================================================
# update model in MRML scene
def update(modelNode):
    """
    Delete 'current' and show 'new' model node based on updated mesh.

    Returns
    -------
    The new model node displayed in the main MRML screen.
    """

    # get model mesh and remove
    volMesh = modelNode.GetMesh()
    slicer.mrmlScene.RemoveNode(modelNode)

    # show new model mesh
    extractSurface = vtk.vtkGeometryFilter()
    extractSurface.SetInputData(volMesh)
    extractSurface.Update()
    modelNode = slicer.modules.models.logic().AddModel(extractSurface.GetOutput())
    modelNode.GetDisplayNode().SetOpacity(opacity)

    deleteSegmentations() # segmentation nodes
    # segmentationExist = deleteSegmentations()
    # if segmentationExist: 
    #     segmentation()

    # updated model
    noOfModelNode = slicer.mrmlScene.GetNumberOfNodesByClass('vtkMRMLModelNode')
    updatedModelNode = slicer.mrmlScene.GetNthNodeByClass(noOfModelNode - 1, 'vtkMRMLModelNode')

    return(updatedModelNode)

# ====================================================================================================
# delete all models
def deleteModels():
    """
    Delete all existing model nodes.

    Returns
    -------
    All models the main MRML screen deleted.
    modelNodeExists : True / False
        The logical variable if any model exists.
    """

    modelNodeExists = False # model existence logical variable
    noOfmodelNode = slicer.mrmlScene.GetNumberOfNodesByClass('vtkMRMLModelNode')
    for i in reversed(range(noOfmodelNode)): # must be from the last, otherwise problem in loop 
                                             # appears, because by deleting model, the length
                                             # of range(noOfModelNode) is one element shorter
        modelNode = slicer.mrmlScene.GetNthNodeByClass(i, 'vtkMRMLModelNode')
        if modelNode.GetName().find('Model') > -1:
            modelNodeExists = True
            slicer.mrmlScene.RemoveNode(modelNode)
    
    return(modelNodeExists)

# ====================================================================================================
# delete all segmentations
def deleteSegmentations():
    """
    Delete all existing segmentations.

    Returns
    -------
    All segmentations the main MRML screen deleted.
    segmentationNodeExists : True / False
        The logical variable if any segmentation exists.
    """

    segmentationNodeExists = False # segmentation existence logical variable
    noOfSegmentationNode = slicer.mrmlScene.GetNumberOfNodesByClass('vtkMRMLSegmentationNode')
    for i in reversed(range(noOfSegmentationNode)): # must be from the last, otherwise problem in loop 
                                                    # appears, because by deleting model, the length
                                                    # of range(noOfSegmentationNode) is one element shorter
        segmentationNode = slicer.mrmlScene.GetNthNodeByClass(i, 'vtkMRMLSegmentationNode')
        if segmentationNode.GetName().find('Segmentation') > -1:
            segmentationNodeExists = True
            slicer.mrmlScene.RemoveNode(segmentationNode)

    return(segmentationNodeExists)

# ====================================================================================================
# delete all markups
def deleteMarkups():
    """
    Delete all existing markups nodes.

    Returns
    -------
    All markups in the main MRML screen deleted.
    markupsNodeExists : True / False
        The logical variable if any markups node exists.
    """

    markupsNodeExists = False # object existence logical variable
    noOfMarkupsNode = slicer.mrmlScene.GetNumberOfNodesByClass('vtkMRMLMarkupsNode')
    if noOfMarkupsNode > 0:
        markupsNodeExists = True
        for i in reversed(range(noOfMarkupsNode)): # must be from the last, otherwise problem in loop 
                                                   # appears, because by deleting model, the length
                                                   # of range(noOfMarkupsNode) is one element shorter
            markupsNode = slicer.mrmlScene.GetNthNodeByClass(i, 'vtkMRMLMarkupsNode')
            slicer.mrmlScene.RemoveNode(markupsNode)

    return(markupsNodeExists)
