# ====================================================================================================
from __main__ import qt
from PelvicFloorLib import QMessageTitle, opacity, pathName, fileName
import os, json, slicer, vtk, numpy as np
from pathlib import Path
from vtk.util import numpy_support

# ====================================================================================================
# read model
def read(pathFileName, modelName, reader):
    """
    Read model from VTK (polydata or unstructured grid) or PC

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

    if reader == 'pc': # convert PC -> VTK
        nodeNumbers, nodeCoordinates, _, elementParts, elementNodes, _, _ = readMeshPC2Tria(pathFileName)
        ren = renumber(nodeNumbers, np.arange(len(nodeNumbers))) # renumber model
        elementNodes = ren[elementNodes].tolist() # tolist() necessary for pyvtk
        _ = renumberLandmarks(pathFileName, nodeNumbers) # renumber landmarks if exist
        
        # write VTK only if does not exist
        if not os.path.exists(str(pathName) + '/' + Path(pathFileName).stem + '.vtk'):
            writeMeshVTK(pathFileName, nodeCoordinates, elementParts, elementNodes) # write VTK
        
        # further work with VTK
        pathFileName = pathFileName[: -3] + '.vtk'
        reader = 'vtk'

    # volMesh = vtk.vtkUnstructuredGrid()
    if reader == 'vtk':
        modelReader = vtk.vtkUnstructuredGridReader()
    elif reader == 'vtu':
        modelReader = vtk.vtkXMLUnstructuredGridReader()
    elif reader ==  'vtp':
        modelReader = vtk.vtkPolyDataReader()
        
    try:
        modelReader.SetFileName(pathFileName)
        modelReader.Update()

        # make model active (necessary to access data)
        slicer.modules.models.logic().AddModel(modelReader.GetOutputPort())

        # number of model nodes
        noOfModelNode = slicer.mrmlScene.GetNumberOfNodesByClass('vtkMRMLModelNode')
            
        # last model node (last read)
        modelNode = slicer.mrmlScene.GetNthNodeByClass(noOfModelNode - 1, 'vtkMRMLModelNode')

    except:
        modelNode = []

    return(modelNode)

def readMeshPC2Tria(fileName):
    """
    Read PC file (with possible INCLUDE files recursively).
    Returns nodeNumbers, nodeCoordinates, elementNumbers, elementParts, elementNodes, boundaryNodes, boundaryConditions
    """
    fileIn = open(fileName, 'r')
    lines = fileIn.readlines()

    nodeNumbers = []
    nodeCoordinates = []

    elementNumbers = []
    elementParts = []
    elementNodes = []

    boundaryNodes = []
    boundaryConditions = []

    maximumElementNumber = 0

    # path to file
    path = fileName.rsplit('/', 1)[:-1]
    path = path[0] + '/' if path else ''

    # read by lines
    for line in lines:
        if len(line.split()) == 0:
            continue
        keyword = line.split()[0]

        if keyword == 'INCLU': # recursive include
            includeFile = path + line.split()[2]
            nodeNumbersi, nodeCoordinatesi, _, elementPartsi, elementNodesi, _, _ = readMeshPC2Tria(includeFile)

            nodeNumbers += nodeNumbersi
            nodeCoordinates += nodeCoordinatesi.tolist() if isinstance(nodeCoordinatesi, np.ndarray) else nodeCoordinatesi
            elementParts += elementPartsi
            elementNodes += elementNodesi.tolist() if isinstance(elementNodesi, np.ndarray) else elementNodesi

        elif keyword == 'NODE':
            nodeNumbers.append(int(line[8:16]))
            nodeCoordinates.append([float(line[16:32]),
                                    float(line[32:48]),
                                    float(line[48:64])])
        elif keyword in ('SHELL','MEMBR'):
            elemNum = int(line[8:16])
            elemPart = int(line[16:24])
            if maximumElementNumber < elemNum:
                maximumElementNumber = elemNum

            if len(line) <= 49: # triangle (line[49] = '\n')
                elementNumbers.append(elemNum)
                elementParts.append(elemPart)
                elementNodes.append([int(line[24:32]),
                                     int(line[32:40]),
                                     int(line[40:48])])
            elif line[55] == '0': # triangle (4 nodes)
                elementNumbers.append(elemNum)
                elementParts.append(elemPart)
                elementNodes.append([int(line[24:32]),
                                     int(line[32:40]),
                                     int(line[40:48])])
            else:  # quad -> 2 triangles
                # first
                elementNumbers.append(elemNum)
                elementParts.append(elemPart)
                elementNodes.append([int(line[24:32]),
                                     int(line[32:40]),
                                     int(line[40:48])])
                # second
                elementNumbers.append(elemNum)
                elementParts.append(elemPart)
                elementNodes.append([int(line[24:32]),
                                     int(line[40:48]),int(line[48:56])])
        elif keyword == 'BOUNC':
            boundaryNodes.append(int(line[8:16]))
            boundaryConditions.append(line[16:32])

    # renumber duplicated triangles after splitting quads
    for i in range(len(elementNumbers)-1):
        if elementNumbers[i] == elementNumbers[i+1]:
            elementNumbers[i+1] += maximumElementNumber

    fileIn.close()

    # convert to NumPy array only where needed
    nodeCoordinates = np.array(nodeCoordinates)
    elementNodes = np.array(elementNodes)

    return (nodeNumbers, nodeCoordinates,
            elementNumbers, elementParts, elementNodes,
            boundaryNodes, boundaryConditions)

# ====================================================================================================
# write mesh to VTK file
def writeMeshVTK(pathFileName, nodeCoordinates, elementParts, elementNodes):
    """
    Write mesh to VTK polydata and unstructured grid 'fileName'.

    Returns
    -------
    File saved to 'fileName' VTK, VTU and VTP.
    """     
    
    # import pyvtk as ptk
    # unstructuredGrid = ptk.UnstructuredGrid(points=nodeCoordinates, triangle=elementNodes)
    # polyData = ptk.PolyData(points=nodeCoordinates, polygons=elementNodes)
    #
    # celldata = ptk.CellData(ptk.Scalars(elementParts, name='Labels'))
    #
    # vtu = ptk.VtkData(unstructuredGrid, celldata, 'Pelvic floor model')
    # vtp = ptk.VtkData(polyData, celldata, 'Pelvic floor model')
    #
    # vtu.tofile(fileName, 'ascii') # binary not readable by meshio
    # vtu.tofile(fileName + '.vtu', 'ascii') # binary not readable by meshio
    # vtp.tofile(fileName + '.vtp', 'ascii') # binary not readable by meshio

    coordinates = np.asarray(nodeCoordinates) # (n, 3)
    triangles = np.asarray(elementNodes) # (m, 3)
    labels = np.asarray(elementParts) # (m,)

    points = vtk.vtkPoints() # body
    points.SetData(numpy_support.numpy_to_vtk(coordinates))
    n_cells = triangles.shape[0] # connectivity

    cells = np.hstack([np.full((n_cells, 1), 3), triangles]).astype(np.int64)
    cells_flat = cells.ravel()

    cells = vtk.vtkCellArray()
    cells.SetCells(n_cells, numpy_support.numpy_to_vtkIdTypeArray(cells_flat, deep=True))

    # unstructured grid
    grid = vtk.vtkUnstructuredGrid()
    grid.SetPoints(points)
    grid.SetCells(vtk.VTK_TRIANGLE, cells)

    # cell labels
    labels_vtk = numpy_support.numpy_to_vtk(labels)
    labels_vtk.SetName("Labels")
    grid.GetCellData().AddArray(labels_vtk)

    # writer to ASCII
    writer = vtk.vtkUnstructuredGridWriter()
    writer.SetFileName(pathFileName[: -3] + '.vtk') # from PC => [: -3]
    writer.SetInputData(grid)
    writer.SetFileTypeToASCII()
    writer.Write()

# ====================================================================================================
# write tetrahedral mesh mesh to PC file
def export(outputFileName, volMesh):
    """
    Write segmented tetrahedral mesh to 'fileName' in the structure of template file.

    Returns
    -------
    The file save to 'fileName' in the PC format.
    """     

    # get reference PC file from text node
    nodeText = slicer.util.getNode("pathFileName")
    pathFileName = nodeText.GetText()

    # nodes to be saved
    nodeCoordinates = numpy_support.vtk_to_numpy(volMesh.GetPoints().GetData())

    #writeMeshPC(pathFileName, outputFileName, nodeCoordinates)
    writeMeshPC(str(Path(pathFileName).parent / Path(pathFileName).stem) + '.pc', outputFileName, nodeCoordinates)

    #fileOut = open(outputFileName + '.pc', 'w')
    #writeMeshPC(pathFileName, fileOut, nodeCoordinates)
    #fileOut.close()

def writeMeshPC(inputFileName, outputFileName, nodeCoordinates, index = 0):
    """
    Write nodes in recursive includes.

    Returns
    -------
    The file save to 'fileName' with all includes in the PC format.
    """     
    
    # backup file if exists
    pathIn = Path(inputFileName)
    pathOut = Path(outputFileName)
    if pathIn.resolve() == pathOut.resolve():
        backupFileName = str(pathIn) + '.bak'
        backup = Path(backupFileName)
        if backup.exists():
            backup.unlink() # remove existing backup
        pathIn.rename(backup)
    else:
        backupFileName = str(pathIn)
    
    fileOut = open(outputFileName, 'w')
    fileIn = open(backupFileName, 'r')
    lines = fileIn.readlines() # read file

    for line in lines:
        keywords = line.split()
        if not keywords:
            fileOut.write(line)
            continue

        if keywords[0] == "INCLU":
            fileOut.write(line)
            includeIn = str(pathIn.parent / keywords[2])
            includeOut = str(pathOut.parent / keywords[2])
            index = writeMeshPC(includeIn, includeOut, nodeCoordinates, index)

        elif keywords[0] == "NODE":
            fileOut.write(line[0 : 16] + "%16.3f%16.3f%16.3f\n" % (nodeCoordinates[index][0], 
                                                                   nodeCoordinates[index][1], 
                                                                   nodeCoordinates[index][2]))
            index += 1 # next node

        else:
            fileOut.write(line)

    fileOut.close()
    fileIn.close()

    return index

# def writeMeshPC(inputFileName, fileOut, nodeCoordinates, index = 0):
#     """
#     Write nodes in recursive includes.

#     Returns
#     -------
#     The file save to 'fileName' with all includes in the PC format.
#     """     
    
#     # pathFileName include path
#     fileIn = open(inputFileName, 'r')
#     lines = fileIn.readlines() # read file

#     for line in lines:
#         keywords = line.split()
#         if not keywords:
#             fileOut.write(line)
#             continue

#         if keywords[0] == "INCLU":
#             path = inputFileName.rsplit('/', 1)[:-1]
#             path = path[0] + '/' if path else ''
#             includeIn = path + keywords[2]

#             #path = outputFileName.rsplit('/', 1)[:-1]
#             #path = path[0] + '/' if path else ''
#             #includeOut = path + keywords[2]
            
#             index = writeMeshPC(includeIn, fileOut, nodeCoordinates, index)

#         elif keywords[0] == "NODE":
#             fileOut.write(line[0 : 16] + "%16.3f%16.3f%16.3f\n" % (nodeCoordinates[index][0], 
#                                                                    nodeCoordinates[index][1], 
#                                                                    nodeCoordinates[index][2]))
#             index += 1 # next node

#         else:
#             fileOut.write(line)

#     fileIn.close()

#     return index

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
def renumberLandmarks(fileName, nodeNum):
    """
    Renumber landmarks.

    Returns
    -------
    landmarksDictionaryRenumbered : disctionary
        The landmarks renumbered according to node numbers.
    """
    
    # model landmark coordinates in renumbered model numbering
    landmarksDictionaryRenumbered = {}

    try: # landmark dictionary exists
        # model landmark coordinates in original model numbering
        with open(fileName + '.dic') as landmarks:
            landmarksDictionary = json.loads(landmarks.read())

        for i, key in enumerate(landmarksDictionary.keys()):
            if landmarksDictionary[key] in nodeNum:
                landmarksDictionaryRenumbered.update({key: nodeNum.index(landmarksDictionary[key])})

        # save landmarks (json to convert apostrophes to quotes)
        with open(fileName + '_landmarks.dic', 'w') as data:
            data.write(str(json.dumps(landmarksDictionaryRenumbered))) # dumps to exchange ' by "

    except: # landmark dictionary does not exists
        pass # nothing happens

    return landmarksDictionaryRenumbered

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
# transformation matrix
def landmarkPairs(modelNode):
    """
    Corresponding landmarks on template model and MRI data.
    
    Returns
    -------
    modelDictNodes : dictionary
        The dictionary of model landmark nodes.
    modelDictCoordinates : dictionary
        The dictionary of model landmark coordinates.
    imageDictCoordinates : dictionary
        The dictionary of MRI landmark coordinates.
    """

    # MRI landmarks coordinates (existing fiducial nodes in MRML screen are expected)
    imageDictCoordinates = {} # dictionary of MRI landmarks coordinates
    noOfMarkupsNode = slicer.mrmlScene.GetNumberOfNodesByClass('vtkMRMLMarkupsNode')

    # actual dictionary of landmark coordinates in active renumbered model
    with open(pathName + fileName + '_landmarks.dic') as landmarks:
        modelDictCoordinatesAll = json.loads(landmarks.read()) # renumbered and converted to dict

    modelDictNodes = {} # dictionary of model landmarks nodes
    modelDictCoordinates = {} # dictionary of model landmarks coordinates
    points = numpy_support.vtk_to_numpy(modelNode.GetMesh().GetPoints().GetData())

    # markupsNode remain in the database after deleting including all properties,
    # therefore reversed loop search from the newest one
    for i in reversed(range(noOfMarkupsNode)): # from the newest to find the newest

        markupsNode = slicer.mrmlScene.GetNthNodeByClass(i, 'vtkMRMLMarkupsNode')
        landmark = markupsNode.GetName()

        if landmark in modelDictCoordinatesAll.keys(): # landmark exists
            
            # add MRI and model landmarks coordinates and model landmarks nodes
            imageDictCoordinates.update({landmark: np.array(markupsNode.GetNthControlPointPosition(0))})
            modelDictCoordinates.update({landmark: points[modelDictCoordinatesAll[landmark]]})
            modelDictNodes.update({landmark: modelDictCoordinatesAll[landmark]})

        else: # not in model markpus, looking for pair

            for j in reversed(range(noOfMarkupsNode)): # from the newest to find the newest

                pairMarkupsNode = slicer.mrmlScene.GetNthNodeByClass(j, 'vtkMRMLMarkupsNode')
                pairLandmark = pairMarkupsNode.GetName()
                if  pairLandmark not in modelDictCoordinatesAll.keys() and \
                    len(pairLandmark) > 2 and pairLandmark[2 : ] == landmark and \
                        pairLandmark[0 : 2] == 'M_': # must begin with 'M_'

                    # add pair landmarks and landmark node number
                    imageDictCoordinates.update({pairLandmark[2 : ]: \
                                                 np.array(markupsNode.GetNthControlPointPosition(0))})
                    modelDictCoordinates.update({pairLandmark[2 : ]: \
                                                 np.array(pairMarkupsNode.GetNthControlPointPosition(0))})
                
                    break # pair found, out of loop

    return(modelDictCoordinates, imageDictCoordinates, modelDictNodes)

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

    # # updated landmarks
    # deleteMarkups()
    # with open(pathName + fileName + '_landmarks.dic') as landmarks:
    #     modelDictCoordinatesAll = json.loads(landmarks.read()) # renumbered and converted to dict

    # # show landmarks
    # points = numpy_support.vtk_to_numpy(volMesh.GetPoints().GetData())
    # for key, value in modelDictCoordinatesAll.items():
    #     markupsNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLMarkupsFiducialNode")
    #     markupsNode.SetName(key)
        
    #     display = markupsNode.GetDisplayNode()
    #     display.SetSelectedColor(1, 1, 1)
    #     display.SetColor(0, 1, 0)

    #     slicer.util.updateMarkupsControlPointsFromArray(markupsNode, points[value].reshape(1,3))

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
