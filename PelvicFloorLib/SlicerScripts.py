# ====================================================================================================
from __main__ import qt
from PelvicFloorLib import linePrint, timePrint, randomColor, \
    read, export, createMesh, separate, getLandmarks, landmarkPairs, \
        align, move, rotate, scale, update, \
            deleteModels, deleteSegmentations, deleteMarkups, RBF, \
                    rotationTransformMatrix, \
                        pelvicFloorVersion, QMessageTitle, epsilon, alpha, morphingType, \
                            pathName, fileName, extension
from scipy.interpolate import CubicSpline, splprep, splev
from scipy.optimize import least_squares
import os, json, slicer, time, vtk, numpy as np
from pathlib import Path
from scipy.interpolate import splprep, splev
from vtk import vtkTransform, vtkMatrix4x4
from vtk.util import numpy_support
from DICOMLib import DICOMUtils

# ====================================================================================================
# import VTK or PC
def importButton():

    if pelvicFloorVersion == 'development': # check for version
        currentTime = linePrint("Reading template model...")

    # open dialog
    fileDialog = qt.QFileDialog()
    pathFileName = fileDialog.getOpenFileName(fileDialog, '', pathName, '*.pc *.vtk *.vtu *.vtp')

    # update text node
    nodeText = slicer.util.getNode("pathFileName")
    nodeText.SetText(pathFileName)

    # read model
    extension = Path(pathFileName).suffix.lstrip('.')
    modelNode = read(pathFileName, 'Model', extension)

    if modelNode != []:
        volMesh = modelNode.GetMesh()
        points = numpy_support.vtk_to_numpy(volMesh.GetPoints().GetData())
        points[:, [0, 1]] *= -1 # -> convert LPS to (internal) RAS coordinate system
        update(modelNode)
        qt.QMessageBox.information(slicer.util.mainWindow(), 
                                QMessageTitle, 'Model read.')
    else:
        qt.QMessageBox.information(slicer.util.mainWindow(), 
                                QMessageTitle, 'Nothing read.')

    if pelvicFloorVersion == 'development': # check for version
        timePrint("... done", currentTime)


# ====================================================================================================
# read template model
def templateButton():
    """
    Read template model. Before reading delete all existing models.
    
    Returns
    -------
    The template model displayed in the main MRML screen.
    """

    if pelvicFloorVersion == 'development': # check for version
        currentTime = linePrint("Reading template model...")

    # delete existing models and segmentations
    deleteModels() # model nodes
    deleteSegmentations() # segmentation nodes

    # template file name
    pathFileName = str(pathName) + '\\' + fileName + extension

    # update text node
    nodeText = slicer.util.getNode("pathFileName")
    nodeText.SetText(pathFileName)

    # read model
    modelNode = read(pathFileName, 'Template model', 'vtk')

    # convert to LPS coordinate system (L, P, S) = (-z, x, -y)
    # A/P means anterior/posterior, L/R means left/right, S/I means superior/inferior
    # LPS coordinate system = (Left, Posterior, Superior)
    # RAS coordinate system = (Right, Anterior, Superior)
    # https://discourse.slicer.org/t/model-files-are-now-saved-in-lps-coordinate-system/10446
    # volMesh.GetPoints().GetPoint(i) # i-th point coordinates as tuple
    # volMesh = slicer.util.getNode(modelNode.GetName()).GetMesh()
    # points = numpy_support.vtk_to_numpy(volMesh.GetPoints().GetData())
    volMesh = modelNode.GetMesh()
    points = numpy_support.vtk_to_numpy(volMesh.GetPoints().GetData())
    
    # 'templateModel' in LPS
    # -> convert LPS to (internal) RAS coordinate system
    points[:, [0, 1]] *= -1

    # show template model
    update(modelNode)

    if pelvicFloorVersion == 'development': # check for version
        timePrint("... done", currentTime)

    qt.QMessageBox.information(slicer.util.mainWindow(), 
                               QMessageTitle, 'Template model read.')

# ====================================================================================================
# read MRI
def imagesButton():
    """
    Read MRI.
    
    Returns
    -------
    MRI displayed in the main MRML screen.
    """

    if pelvicFloorVersion == 'development': # check for version
        currentTime = linePrint("Reading MRI...")
    
    # open dialog
    fileDialog = qt.QFileDialog()
    # dicomDataDir = fileDialog.getOpenFileName(fileDialog, '', pathName, '*.dcm')
    dicomDataDir = fileDialog.getExistingDirectory(fileDialog, '', pathName, fileDialog.ShowDirsOnly)
    # dicomDataDir = 'C:/Users/hyncik/Work/HBM/Models/Pelvic_floor/MRI/Nulipary/Anonymized_to_UWA/MRI/P011'

    if not dicomDataDir == '':
        try: # always works, if 'dicomDataDir' is wrong, uses study instance
            DICOMUtils.importDicom(dicomDataDir)
            patientUIDs = slicer.dicomDatabase.patients()
            loadedNodeIDs = [] # this list will contain the list of all loaded node IDs
        
            for patientUID in patientUIDs:
                loadedNodeIDs.extend(DICOMUtils.loadPatientByUID(patientUID))  

            qt.QMessageBox.information(slicer.util.mainWindow(), 
                                    QMessageTitle, 'MRI read.')
        except:
            qt.QMessageBox.information(slicer.util.mainWindow(), 
                                    QMessageTitle, 'MRI do not exist.')
    else:
        qt.QMessageBox.information(slicer.util.mainWindow(), 
                                QMessageTitle, 'MRI do not exist.')

    if pelvicFloorVersion == 'development': # check for version
        timePrint("... done", currentTime)

# ====================================================================================================
# export model
def exportButton(file):
    """
    Export model to PC and VTK as unstructured grid.
    Rewritting the imported file name creates backup files, when PC exists,
    because it copies the files line by line and just updates nodal cooordinates.
    Rewriting any file name does not create any backup of the original file.
    
    Returns
    -------
    The model exported to PC and VTK.
    """

    if pelvicFloorVersion == 'development': # check for version
        currentTime = linePrint("Exporting model...")

    noOfModelNode = slicer.mrmlScene.GetNumberOfNodesByClass('vtkMRMLModelNode')
    if noOfModelNode > 0:
        
        # apply to last model only
        modelNode = slicer.mrmlScene.GetNthNodeByClass(noOfModelNode - 1, 'vtkMRMLModelNode')
        if modelNode.GetName().find('Model') == 0: # last name 'Model*' exist

            # save file to VTK and to PC (if PC exists)
            if file == '': # for export button
                fileDialog = qt.QFileDialog()
                outputFileName = fileDialog.getSaveFileName(fileDialog, '', pathName, '*.pc *.vtk')
                #if os.path.exists(outputFileName):
                #    outputFileName = ''
                #    qt.QMessageBox.warning(slicer.util.mainWindow(), 
                #                        QMessageTitle, 'File exists.')
            else: # for users scripts
                outputFileName = file # VTK

            path = Path(outputFileName)
            #print(path.parent)  # direcotry
            #print(path.stem)    # model
            #print(path.suffix)  # suffix

            if not outputFileName == '': # file name chosen => export
                # get reference PC file from text node
                nodeText = slicer.util.getNode("pathFileName")
                pathFileName = nodeText.GetText()

                # get points from model node
                points = numpy_support.vtk_to_numpy(modelNode.GetMesh().GetPoints().GetData())
                points[:, [0, 1]] = -points[:, [0, 1]] # RAS -> LPS for exporting

                # PolyData => UnstructuredGrid
                polyData = modelNode.GetPolyData()
                ugrid = vtk.vtkUnstructuredGrid()
                ugrid.SetPoints(polyData.GetPoints())

                # faces => Cells
                cells = polyData.GetPolys()
                ugrid.SetCells(vtk.VTK_TRIANGLE, cells)

                # export to VTK file (unstructured ASCII)
                writer = vtk.vtkUnstructuredGridWriter()
                writer.SetFileName(str(path.parent / path.stem) + '.vtk')
                writer.SetInputData(ugrid)
                #writer.SetFileTypeToBinary()
                writer.SetFileTypeToASCII()
                writer.Write()

                # export to PC file (RAS coordinate system) if reference PC exists
                if os.path.exists(str(Path(pathFileName).parent / Path(pathFileName).stem) + '.pc'):
                    export(str(path.parent / path.stem) + '.pc', modelNode.GetMesh())

                # return LPS -> RAS back for MRML
                points[:, [0, 1]] = -points[:, [0, 1]]

                qt.QMessageBox.information(slicer.util.mainWindow(), 
                                        QMessageTitle, 'Model exported.')
            else: # no file chosen
                pass
        else:
            qt.QMessageBox.warning(slicer.util.mainWindow(), 
                                   QMessageTitle, 'Nothing to export.')
    else:
        qt.QMessageBox.warning(slicer.util.mainWindow(), 
                               QMessageTitle, 'Nothing to export.')
    
    if pelvicFloorVersion == 'development': # check for version
        timePrint("... done", currentTime)

# ====================================================================================================
# rotate model around L
def rotateButton(axis):
    """
    Rotate model by 90 degrees around axis.
    
    Returns
    -------
    The rotated model displayed the main MRML screen.
    """

    if pelvicFloorVersion == 'development': # check for version
        currentTime = linePrint("Rotating model...")

    # check existence of model
    noOfModelNode = slicer.mrmlScene.GetNumberOfNodesByClass('vtkMRMLModelNode')
    if noOfModelNode > 0:
        
        # apply to last model only
        modelNode = slicer.mrmlScene.GetNthNodeByClass(noOfModelNode - 1, 'vtkMRMLModelNode')
        if modelNode.GetName().find('Model') == 0: # last name 'Model*' exist
            
            # get current model mesh points
            points = numpy_support.vtk_to_numpy(modelNode.GetMesh().GetPoints().GetData())

            if axis == 'L': # rotate around L
                points[:, [1, 2]] = points[:, [2, 1]]
                points[:, 2] = -points[:, 2]
            if axis == 'P': # rotate around P
                points[:, [0, 2]] = points[:, [2, 0]]
                points[:, 0] = -points[:, 0]
            if axis == 'S': # rotate around S
                points[:, [0, 1]] = points[:, [1, 0]]
                points[:, 0] = -points[:, 0]
            
            # update model node
            update(modelNode)
          
            qt.QMessageBox.information(slicer.util.mainWindow(), 
                                        QMessageTitle, 'Model rotated around L.')
        else:
            qt.QMessageBox.warning(slicer.util.mainWindow(), 
                                   QMessageTitle, 'Model does not exist.')
    else:
        qt.QMessageBox.warning(slicer.util.mainWindow(), 
                               QMessageTitle, 'Model does not exist.')
    
    if pelvicFloorVersion == 'development': # check for version
        timePrint("... done", currentTime)

# ====================================================================================================
# read model
def landmarksButton():
    """
    Read sample landmarks from CSV, JSON, FCSV or dictionary text files. 
    Dictionary text file is formated {'landmark: [L, P, S], ...}.
    
    Returns
    -------
    The landmarks displayed in the main MRML screen.
    """

    if pelvicFloorVersion == 'development': # check for version
        currentTime = linePrint("Reading landmarks...")

    # read file
    fileDialog = qt.QFileDialog()
    inputFileName = fileDialog.getOpenFileName(fileDialog, '', pathName, 
                                               '*.csv; *.json; *.fcsv; *.txt')
    fileType = inputFileName[-4 :] # postfix
    # fileType = os.path.splitext(inputFileName)[1].lower()

    try: # read merged markups
        if fileType == '.csv': # CSV (including labels) -> reads as a single node
            markupsNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLMarkupsFiducialNode")
            slicer.modules.markups.logic().ImportControlPointsFromCSV(markupsNode, inputFileName)
        elif fileType == 'json': # JSON (not labeled) -> does not read labels
            slicer.util.loadMarkups(inputFileName)
        elif fileType == 'fcsv': # FCSV (not labeled) -> does not read labels, reads in RAS
            #slicer.util.loadMarkupsFiducialList(inputFileName)

            with open(inputFileName, 'r') as lines:
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

        else: # dictionary (including labels) -> reads in LPS
            with open(inputFileName) as landmarks:
                landmarksDictXyz = json.loads(landmarks.read())

            for landmark in landmarksDictXyz: # separate markups
                node = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLMarkupsFiducialNode")
                node.SetName(landmark)
                display = node.GetDisplayNode()
                display.SetSelectedColor(1, 0, 0) # displayed in red
                display.SetColor(0, 1, 0) # green when selected

                # fiducial node can have more points 
                # -> [[x, y, z]]
                landmarkXyz = np.array([landmarksDictXyz[landmark]])
                landmarkXyz[0, [0, 1]] = -landmarkXyz[0, [0, 1]] # LSP -> RAS
                slicer.util.updateMarkupsControlPointsFromArray(node, landmarkXyz)

            qt.QMessageBox.information(slicer.util.mainWindow(), 
                                       QMessageTitle, 'Landmarks read.')

        # all except dictionary
        if fileType in ('.csv', 'json', 'fcsv'):

            # separate last (read) markups node
            noOfMarkupsNode = slicer.mrmlScene.GetNumberOfNodesByClass('vtkMRMLMarkupsFiducialNode')
            markupsNode = slicer.mrmlScene.GetNthNodeByClass(noOfMarkupsNode - 1, \
                                                            'vtkMRMLMarkupsFiducialNode')
            
            # number of fiducial points in markups node
            noOfMarkupsNodePoints = markupsNode.GetNumberOfControlPoints()
            if noOfMarkupsNodePoints > 0: # fiducial points exist

                # landmarksDictXyz = {} # JSON
                for landmark in range(noOfMarkupsNodePoints):
                    
                    # new node
                    node = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLMarkupsFiducialNode")
                    
                    if fileType == '.csv': # [: -2]) # avoid '-1'
                        node.SetName(markupsNode.GetNthControlPointLabel(landmark))
                    elif fileType == 'json':
                        # node.SetName(landmarks[landmark]) # JSON
                        # landmarksDictXyz.update({landmarks[landmark]: \ 
                        #         np.array(markupsNode.GetNthControlPointPosition(landmark))})
                        pass
                    elif fileType == 'fcsv':
                        node.SetName(landmarks[landmark])

                    # add point to new node
                    display = node.GetDisplayNode()
                    display.SetSelectedColor(1, 0, 0) # displayed in red
                    display.SetColor(0, 1, 0) # green when selected
                    slicer.util.updateMarkupsControlPointsFromArray(node, 
                            np.array([markupsNode.GetNthControlPointPosition(landmark)]))

                # remove merged markups
                slicer.mrmlScene.RemoveNode(markupsNode)

                qt.QMessageBox.information(slicer.util.mainWindow(), 
                                            QMessageTitle, 'Landmarks read.')

            else:
                qt.QMessageBox.information(slicer.util.mainWindow(), 
                                        QMessageTitle, 'Landmarks do not exist.')

    except:
        qt.QMessageBox.information(slicer.util.mainWindow(), 
                                QMessageTitle, 'Landmarks do not exist.')

    if pelvicFloorVersion == 'development': # check for version
        timePrint("... done", currentTime)

# ====================================================================================================
# align model
def alignButton(landmark):
    """
    Align model to 'landmark'.
    
    Returns
    -------
    The aligned model displayed in the main MRML screen.
    If landmark 'landmark' is not present, nothing happens and warning is printed.
    """

    if pelvicFloorVersion == 'development': # check for version
        currentTime = linePrint("Model aligning...")

    # check existence of model
    noOfModelNode = slicer.mrmlScene.GetNumberOfNodesByClass('vtkMRMLModelNode')
    if noOfModelNode > 0:

        # apply to last model only
        modelNode = slicer.mrmlScene.GetNthNodeByClass(noOfModelNode - 1, 'vtkMRMLModelNode')
        if modelNode.GetName().find('Model') == 0: # last name 'Model*' exist

            PSIExists = align(modelNode, 'PSI')
            
            if PSIExists:
                view = slicer.app.layoutManager().threeDWidget(0).threeDView()
                view.cornerAnnotation().SetText(vtk.vtkCornerAnnotation.UpperRight,
                                                "Aligned model") # set text
                qt.QMessageBox.information(slicer.util.mainWindow(), 
                                           QMessageTitle, 'Model aligned.')
            else:
                pass # message printed by get_landmarks
        else:
            qt.QMessageBox.warning(slicer.util.mainWindow(), 
                                QMessageTitle, 'Model does not exist.')
    else:
        qt.QMessageBox.warning(slicer.util.mainWindow(), 
                               QMessageTitle, 'Model does not exist.')
    
    if pelvicFloorVersion == 'development': # check for version
        timePrint("... done", currentTime)

# ====================================================================================================
# mesh scaling
def scaleButton():
    """
    Scale in two dimensions.
    
    Returns
    -------
    Scaled mesh based on existing landmarks.
    """

    if pelvicFloorVersion == 'development': # check for version
        currentTime = linePrint("Mesh scaling...")

    noOfModelNode = slicer.mrmlScene.GetNumberOfNodesByClass('vtkMRMLModelNode')
    if noOfModelNode > 0:

        # apply to last model
        modelNode = slicer.mrmlScene.GetNthNodeByClass(noOfModelNode - 1, 'vtkMRMLModelNode')
        if modelNode.GetName().find('Model') > -1:

            # check existence of landmarks PR, PSP, ILL, ILR
            landmarksDict = {'PSP': False, 'PR': False, 'IL_L': False, 'IL_R': False}
            allExist, modelDictXyz, imageDictXyz = getLandmarks(modelNode, landmarksDict)

            # all necessary landmarks exist
            if allExist:

                # model and image dimensions
                imageAPD = imageDictXyz['PR'] - imageDictXyz['PSP']
                modelAPD = modelDictXyz['PR'] - modelDictXyz['PSP']
                imageTD = imageDictXyz['IL_R'] - imageDictXyz['IL_L']
                modelTD = modelDictXyz['IL_R'] - modelDictXyz['IL_L']

                # scaling vector (expect APD perpendicular to TD)
                scaleAPD = np.linalg.norm(imageAPD) / \
                    np.linalg.norm(modelAPD) # np.divide(imageAPD, modelAPD)
                scaleTD = np.linalg.norm(imageTD) / \
                    np.linalg.norm(modelTD) # np.divide(imageTD, modelTD)

                # expecting APD lise in PS and TD lies in L
                # rotate model and image APD (rotate around PSP)
                angle = -np.arccos(np.dot(imageAPD, modelAPD) / \
                                   np.linalg.norm(imageAPD) / np.linalg.norm(modelAPD))
                modelNode = rotate(modelNode, 'PSP', [1, 0, 0], angle)

                # scale rotated model
                scaleVector = [scaleTD, scaleAPD, scaleAPD] # np.multiply(scaleAPD, scaleTD)
                modelNode = scale(modelNode, scaleVector)
                
                # align to PSP: same as align('PSP'), but without messages
                vector = imageDictXyz['PSP'] - modelDictXyz['PSP']
                move(modelNode, vector)

                view = slicer.app.layoutManager().threeDWidget(0).threeDView()
                view.cornerAnnotation().SetText(vtk.vtkCornerAnnotation.UpperRight,
                                                "Scaled model") # set text

                qt.QMessageBox.information(slicer.util.mainWindow(), 
                                           QMessageTitle, 'Model scaled.')
            else:
                pass # message printed by get_landmarks
        else:
            qt.QMessageBox.warning(slicer.util.mainWindow(), 
                                   QMessageTitle, 'Model does not exist.')
    else:
        qt.QMessageBox.warning(slicer.util.mainWindow(), 
                                QMessageTitle, 'Model does not exist.')

    if pelvicFloorVersion == 'development': # check for version
        timePrint("... done", currentTime)

# ====================================================================================================
# mesh morphing
def registrationButton(method):
    """
    Registration to MRI.
    
    Returns
    -------
    Registered mesh based on existing landmarks.
    """
    
    if pelvicFloorVersion == 'development': # check for version
        if method == 'centre':
            current_time = linePrint("Aligning landmarks...")
        elif method == 'rigid':
            current_time = linePrint("Mesh registering...")
        else: # method == 'morphing'
            current_time = linePrint("Mesh morphing...")

    noOfModelNode = slicer.mrmlScene.GetNumberOfNodesByClass('vtkMRMLModelNode')
    if noOfModelNode > 0:

        # applies to last model
        modelNode = slicer.mrmlScene.GetNthNodeByClass(noOfModelNode - 1, 'vtkMRMLModelNode')
        if modelNode.GetName().find('Model') > -1:
                        
            # corresponding (and additional by 'M_' marked) landmarks 
            # on model and MRI sorted from last to first
            modelDictCoordinates, imageDictCoordinates, modelDictNodes = landmarkPairs(modelNode)
            
            # len(modelDictCoordinates) = len(imageDictCoordinates)
            numberOfPairedLandmarks = len(modelDictCoordinates) 
            
            # modelDictCoordinates.landmarks() = imageDictCoordinates.landmarks()
            landmarks = list(modelDictCoordinates.keys())

            # check duplicate landmarks
            noOfCoincidentModelLandmarks = 0
            noOfCoincidentImageLandmarks = 0
            for i in range(numberOfPairedLandmarks - 1):
                for j in range(i + 1, numberOfPairedLandmarks):
                    pairModel = np.linalg.norm(np.r_[modelDictCoordinates[landmarks[i]]] - 
                                               np.r_[modelDictCoordinates[landmarks[j]]])
                    pairImage = np.linalg.norm(np.r_[imageDictCoordinates[landmarks[i]]] - 
                                               np.r_[imageDictCoordinates[landmarks[j]]])
                    if pairModel == 0: # model landmarks 'i' and 'j' coincident
                        noOfCoincidentModelLandmarks += 1
                    if pairImage == 0: # image landmarks 'i' and 'j' coincident
                        noOfCoincidentImageLandmarks += 1

            if noOfCoincidentModelLandmarks > 0:
                qt.QMessageBox.warning(slicer.util.mainWindow(), QMessageTitle, 
                                        'There are ' + str(noOfCoincidentModelLandmarks) + ' coincident model landmarks.')
            elif noOfCoincidentImageLandmarks > 0:
                qt.QMessageBox.warning(slicer.util.mainWindow(), QMessageTitle, 
                                        'There are ' + str(noOfCoincidentImageLandmarks) + ' coincident image landmarks.')
            else:
                # ----------------------------------------------------------------------------------------------------------------------------
                if method == 'centre': # aligning cetnres of gravity
                                       # -> moves model to align centres of gravity of template and target landmarks sets
                    if numberOfPairedLandmarks > 0:

                        modelX = np.array([modelDictCoordinates[landmark] for landmark in landmarks]) # model landmarks
                        imageX = np.array([imageDictCoordinates[landmark] for landmark in landmarks]) # target landmarks

                        # model mesh
                        AX = imageX.mean(axis=0) - modelX.mean(axis=0) # average point (centre of rotation)
                        volMesh = modelNode.GetMesh() # model coordinates (n, 3)
                        templatePoints = numpy_support.vtk_to_numpy(volMesh.GetPoints().GetData())
                        
                        # pre-registration (translation only)
                        X = templatePoints + AX # mesh
                        modelX = modelX + AX # landmarks

                        # pre-registered model
                        targetPoints = vtk.vtkPoints()
                        for i in range(len(X)):
                            targetPoints.InsertPoint(i, X[i])
                        volMesh.SetPoints(targetPoints)
                        # vtkArray = numpy_support.numpy_to_vtk(X, deep=True)
                        # volMesh.GetPoints().SetData(vtkArray)
                        # volMesh.Modified() # not necessary but better

                        # display preregistered landmarks in green
                        for i in range(len(modelX)):
                            node = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLMarkupsFiducialNode")
                            node.SetName('P_' + str(i))
                            display = node.GetDisplayNode()
                            display.SetSelectedColor(0, 1, 0) # displayed in green
                            display.SetColor(0, 1, 0) # green when selected
                            slicer.util.updateMarkupsControlPointsFromArray(node, np.array([modelX[i]]))

                        # modelNode = update(modelNode)
                        update(modelNode)

                        view = slicer.app.layoutManager().threeDWidget(0).threeDView()
                        view.cornerAnnotation().SetText(vtk.vtkCornerAnnotation.UpperRight,
                                                        "Landmarks aligned") # set text
                            
                        qt.QMessageBox.information(slicer.util.mainWindow(), 
                                                QMessageTitle, 'Landmarks aligned.')

                    else:
                        qt.QMessageBox.warning(slicer.util.mainWindow(), QMessageTitle, 
                                            'At least 1 landmark for pre-registration needed, ' +
                                            str(numberOfPairedLandmarks) + ' given.')            

                # ----------------------------------------------------------------------------------------------------------------------------
                elif method == 'rigid': # mesh rigid registraion by Kabsch algorithm
                                        # -> any two landmarks must not be coincident and 
                                        # -> at least one coincident landmark
                    if numberOfPairedLandmarks > 1:

                        modelX = np.array([modelDictCoordinates[landmark] for landmark in landmarks]) # model landmarks
                        imageX = np.array([imageDictCoordinates[landmark] for landmark in landmarks]) # target landmarks

                        # centroids
                        centroid_model = np.mean(modelX, axis=0)
                        centroid_image = np.mean(imageX, axis=0)

                        # model and target coordinates centered
                        Xm = modelX - centroid_model
                        Xi = imageX - centroid_image

                        # Kabsch algorithm
                        H = Xm.T @ Xi
                        U, S, Vt = np.linalg.svd(H)
                        R = Vt.T @ U.T

                        # correction of potential reflection
                        if np.linalg.det(R) < 0:
                            Vt[2, :] *= -1
                            R = Vt.T @ U.T

                        # model mesh
                        volMesh = modelNode.GetMesh() # model coordinates (n, 3)
                        templatePoints = numpy_support.vtk_to_numpy(volMesh.GetPoints().GetData())
                        # n = len(templatePoints) # number of nodes

                        # rigid registration
                        # X = rigid(res.x, templatePoints + AX)
                        t = centroid_image - R @ centroid_model # rigid translation
                        X = (R @ templatePoints.T).T + t # mesh rigid transform
                        modelX = (R @ modelX.T).T + t # landmarks rigid transform

                        # registered model
                        vtkArray = numpy_support.numpy_to_vtk(X, deep=True)
                        volMesh.GetPoints().SetData(vtkArray)
                        volMesh.Modified() # not necessary but better

                        # display registered landmarks in green
                        for i in range(len(modelX)):
                            node = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLMarkupsFiducialNode")
                            node.SetName('R_' + str(i))
                            display = node.GetDisplayNode()
                            display.SetSelectedColor(0, 1, 0) # displayed in green
                            display.SetColor(0, 1, 0) # green when selected
                            slicer.util.updateMarkupsControlPointsFromArray(node, np.array([modelX[i]]))

                        # modelNode = update(modelNode)
                        update(modelNode)

                        view = slicer.app.layoutManager().threeDWidget(0).threeDView()
                        view.cornerAnnotation().SetText(vtk.vtkCornerAnnotation.UpperRight,
                                                        "Registered model") # set text
                            
                        qt.QMessageBox.information(slicer.util.mainWindow(), 
                                                QMessageTitle, 'Model registered.')

                    else:
                        qt.QMessageBox.warning(slicer.util.mainWindow(), QMessageTitle, 
                                            'At least 2 landmarks for registration needed, ' +
                                            str(numberOfPairedLandmarks) + ' given.')            
                # ----------------------------------------------------------------------------------------------------------------------------
                else: # mesh morphing (method == 'morphing'):
                    # -> any two landmarks must not be coincident and 
                    # -> at least four landmarks not lying in a plane
                    if numberOfPairedLandmarks > 3:

                        # landmarks must not be in plane
                        modelX = np.array([modelDictCoordinates[landmark] for landmark in landmarks]) # model landmarks
                        imageX = np.array([imageDictCoordinates[landmark] for landmark in landmarks]) # target landmarks

                        # unit normal vector and plane coefficient
                        modelNormal = np.cross(modelX[1] - modelX[0], modelX[2] - modelX[0])
                        modelNormal /= np.linalg.norm(modelNormal)
                        imageNormal = np.cross(imageX[1] - imageX[0], imageX[2] - imageX[0])
                        imageNormal /= np.linalg.norm(imageNormal)

                        # check other landmarks not to be in plane
                        noOfModelLandmarksNotOnPlane = 0
                        noOfImageLandmarksNotOnPlane = 0
                        for i in range(3, numberOfPairedLandmarks):
                            if abs(np.dot(modelNormal, np.r_[modelDictCoordinates[landmarks[i]]]) - 
                                np.dot(modelNormal, modelX[0])) > epsilon:
                                noOfModelLandmarksNotOnPlane += 1
                            if abs(np.dot(imageNormal, np.r_[imageDictCoordinates[landmarks[i]]]) - 
                                np.dot(imageNormal, imageX[0])) > epsilon:
                                noOfImageLandmarksNotOnPlane += 1

                        if noOfModelLandmarksNotOnPlane > 0 and noOfImageLandmarksNotOnPlane > 0:

                            m = numberOfPairedLandmarks # number of landmarks

                            volMesh = modelNode.GetMesh() # model coordinates (n, 3)
                            X = numpy_support.vtk_to_numpy(volMesh.GetPoints().GetData())
                            n = len(X) # number of nodes

                            # baseline L (m, 3) and target T (m, 3) landmarks coordinates
                            # np.array([[0] * m] * 3) -> slower
                            L = np.zeros((m, 3))
                            T = np.zeros((m, 3))

                            # modelDictCoordinates = imageDictCoordinates
                            for j, landmark in enumerate(modelDictCoordinates): # for all present landmarks
                                L[j] = modelDictCoordinates[landmark] # baseline landmark L
                                T[j] = imageDictCoordinates[landmark] # target landmark T

                            # distance matrix landmarks
                            D = np.linalg.norm(L[:, None, :] - L[None, :, :], axis=2)
                            radius = np.max(D) * 0.4 # support radius
                            A = RBF(D, morphingType, radius) # RBF matrix
                            A += alpha * np.eye(m) # Tikhonov regularization

                            # baseline coordinates matrix [1, L] (affine part)
                            B = np.hstack((np.ones((m, 1)), L))

                            # transformation matrix
                            M = np.block([[A, B], 
                                        [B.T, np.zeros((4, 4))]])
                            R = np.vstack((T, np.zeros((4, 3))))
                            lambdac = np.linalg.solve(M, R)

                            # morph mesh nodes
                            D1 = np.linalg.norm(X[:,None,:] - L[None,:,:], axis=2)
                            A1 = RBF(D1, morphingType, radius)
                            B1 = np.hstack((np.ones((n, 1)), X))
                            X1 = np.dot(np.hstack((A1, B1)), lambdac)

                            # preserve landmark nodes
                            for landmark,node in modelDictNodes.items():
                                X1[node] = imageDictCoordinates[landmark]

                            # morphed model
                            vtkArray = numpy_support.numpy_to_vtk(X1, deep=True)
                            volMesh.GetPoints().SetData(vtkArray)
                            volMesh.Modified() # not necessary but better

                            # display morphed landmarks in blue
                            points = numpy_support.vtk_to_numpy(volMesh.GetPoints().GetData())
                            for landmark in modelDictNodes: # for all present landmarks
                                node = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLMarkupsFiducialNode")
                                node.SetName('M_' + landmark)
                                display = node.GetDisplayNode()
                                display.SetSelectedColor(0, 0, 1) # displayed in blue
                                display.SetColor(0, 1, 0) # green when selected
                                slicer.util.updateMarkupsControlPointsFromArray(node, np.array([points[modelDictNodes[landmark], :]]))

                            # modelNode = update(modelNode)
                            update(modelNode)

                            view = slicer.app.layoutManager().threeDWidget(0).threeDView()
                            view.cornerAnnotation().SetText(vtk.vtkCornerAnnotation.UpperRight,
                                                            "Morphed model") # set text
                            
                            if pelvicFloorVersion == 'development': # check for version
                                print(" condition number = " + str(np.linalg.cond(M)), end = '')
                            qt.QMessageBox.information(slicer.util.mainWindow(), 
                                                    QMessageTitle, 'Model morphed.')
                        else:
                            if noOfModelLandmarksNotOnPlane == 0:
                                qt.QMessageBox.warning(slicer.util.mainWindow(), 
                                                    QMessageTitle, 'All model landmarks are on single plane.')
                    else:
                        if numberOfPairedLandmarks <= 3:
                            qt.QMessageBox.warning(slicer.util.mainWindow(), QMessageTitle, 
                                                'At least 4 landmarks for morphing needed, ' +
                                                str(numberOfPairedLandmarks) + ' given.')
        else:
            qt.QMessageBox.warning(slicer.util.mainWindow(), 
                                   QMessageTitle, 'Nothing to register.')    
    else:
        qt.QMessageBox.warning(slicer.util.mainWindow(), 
                               QMessageTitle, 'Nothing to register.')    

    if pelvicFloorVersion == 'development': # check for version
        timePrint("... done", current_time)

# ====================================================================================================
# show birth
def birthButton():
    """
    Birth canal view.
    
    Returns
    -------
    The birth canal displayed in the main MRML screen.
    """

    if pelvicFloorVersion == 'development': # check for version
        currentTime = linePrint("Displaying birth canal and birth along curve of Carus...")

    # ------------------------------------------------------------
    # birth canal
    markupsNodes = slicer.util.getNodesByClass("vtkMRMLMarkupsNode")
    allExist = 0 # count on all existing landmarks

    for markupsNode in markupsNodes:

        # ------------------------------------------------------------
        # sequence of if is necessary, because PSI is twice,
        if markupsNode.GetName() == "IL_L": # there is always single landmark in the node
            pInletA = [0, 0, 0]
            markupsNode.GetNthControlPointPositionWorld(0, pInletA)
            allExist += 1

        # otherwise elif is better
        # inlet plane => PR, PSA, IL_L, IL_R (ellipse)
        if markupsNode.GetName() == "IL_R":
            pInletB = [0, 0, 0]
            markupsNode.GetNthControlPointPositionWorld(0, pInletB)
            allExist += 1
    
        if markupsNode.GetName() == "PSA":
            pInletC = [0, 0, 0]
            markupsNode.GetNthControlPointPositionWorld(0, pInletC)
            allExist += 1
    
        if markupsNode.GetName() == "PR":
            pInletD = [0, 0, 0]
            markupsNode.GetNthControlPointPositionWorld(0, pInletD)
            allExist += 1
    
        # ------------------------------------------------------------
        # greatest plane => S3, PSP, AC_L, AC_R (ellipse)
        if markupsNode.GetName() == "AC_L":
            pGreatestA = [0, 0, 0]
            markupsNode.GetNthControlPointPositionWorld(0, pGreatestA)
            allExist += 1
    
        if markupsNode.GetName() == "AC_R":
            pGreatestB = [0, 0, 0]
            markupsNode.GetNthControlPointPositionWorld(0, pGreatestB)
            allExist += 1
    
        if markupsNode.GetName() == "PSP":
            pGreatestC = [0, 0, 0]
            markupsNode.GetNthControlPointPositionWorld(0, pGreatestC)
            allExist += 1
    
        if markupsNode.GetName() == "S3":
            pGreatestD = [0, 0, 0]
            markupsNode.GetNthControlPointPositionWorld(0, pGreatestD)
            allExist += 1
    
        # ------------------------------------------------------------
        # least plane => S5, PSI, IS_L, IS_R (ellipse)
        if markupsNode.GetName() == "IS_L":
            pLeastA = [0, 0, 0]
            markupsNode.GetNthControlPointPositionWorld(0, pLeastA)
            allExist += 1
    
        if markupsNode.GetName() == "IS_R":
            pLeastB = [0, 0, 0]
            markupsNode.GetNthControlPointPositionWorld(0, pLeastB)
            allExist += 1
    
        if markupsNode.GetName() == "PSI":
            pLeastC = [0, 0, 0]
            markupsNode.GetNthControlPointPositionWorld(0, pLeastC)
            allExist += 1
    
        if markupsNode.GetName() == "S5":
            pLeastD = [0, 0, 0]
            markupsNode.GetNthControlPointPositionWorld(0, pLeastD)
            allExist += 1

        # ------------------------------------------------------------
        # outlet plane => PSI, TI_L, TI_R (circle)
        if markupsNode.GetName() == "PSI":
            pOutletA = [0, 0, 0]
            markupsNode.GetNthControlPointPositionWorld(0, pOutletA)
            allExist += 1
    
        if markupsNode.GetName() == "TI_L":
            pOutletB = [0, 0, 0]
            markupsNode.GetNthControlPointPositionWorld(0, pOutletB)
            allExist += 1
    
        if markupsNode.GetName() == "TI_R":
            pOutletC = [0, 0, 0]
            markupsNode.GetNthControlPointPositionWorld(0, pOutletC)
            allExist += 1

    # all landmarks defining pelvic planes exist
    if allExist == 15: # (4 points x 4 ellipses + 3 points x 1 triangle)

        # ellipse from 4 points
        # def ellipse_from_diameters(p1, p2, p3, p4, numberOfPoints=100):
        #     center = (p1 + p2 + p3 + p4) / 4
        #     axis_a = (p1 - p2) / 2
        #     axis_b = (p3 - p4) / 2
        #     theta = np.linspace(0, 2 * np.pi, numberOfPoints)
        #     ellipse = np.array([center + np.cos(t) * axis_a + np.sin(t) * axis_b for t in theta])
        #     return ellipse

        def spatial_ellipse_from_4_points(p1, p2, p3, p4, numberOfPoints=100):

            #points = np.array([p1, p2, p3, p4], dtype=float)
            points = np.array([p1, p3, p2, p4], dtype=float) # diameters (p1, p2) and (p3, p4)

            points_closed = np.vstack([points, points[0]]) # close the curve
            t = np.arange(len(points_closed)) # parameter
            spline, _ = splprep(points_closed.T, u=t, s=0,per=True, k=3) # periodic cubic spline
            tt = np.linspace(0, 4, numberOfPoints) # generate smooth curve
            curve = np.array(splev(tt, spline)).T

            return curve

        ellipses_points = [(np.array(pInletA), np.array(pInletB), np.array(pInletC), np.array(pInletD)), # pelvic inlet (ILL, ILR, PSA, PR)
                           (np.array(pGreatestA), np.array(pGreatestB), np.array(pGreatestC), np.array(pGreatestD)), # greatest plane (S3, PSP, AC_L, AC_R)
                           (np.array(pLeastA), np.array(pLeastB), np.array(pLeastC), np.array(pLeastD)), # least plane (S5, PSI, IS_L, IS_R)
                           None] # pelvic outlet (PSI, TI_L, TI_R)

        # ------------------------------------------------------------
        # diplay allipses
        major_ellipses = []
        for points in ellipses_points[:3]:
            ellipse = spatial_ellipse_from_4_points(points[0], points[1], points[2], points[3])
            major_ellipses.append(ellipse)

        for idx, ellipse in enumerate(major_ellipses):
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

            modelNodeEllipse = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLModelNode", f"Major Ellipse {idx + 1}")
            modelNodeEllipse.SetAndObservePolyData(polydataEllipse)
            displayNodeEllipse = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLModelDisplayNode")

            displayNodeEllipse.SetColor(1, 0, 0)
            displayNodeEllipse.SetLineWidth(3)
            displayNodeEllipse.SetOpacity(1.0)

            modelNodeEllipse.SetAndObserveDisplayNodeID(displayNodeEllipse.GetID())

        # ------------------------------------------------------------
        # parameter along pelvis
        t = np.array([0.0, 0.5, 1.0]) 
        numberOfLayers = 50
        ti = np.linspace(0.0, 1.0, numberOfLayers)
        ellipses_smooth = [] # interpolate the four anatomical points in 3D

        for tii in ti:

            interpolated_points = []
            for point in range(4):

                # coordinates of this anatomical point
                points = np.array([ellipses_points[0][point], ellipses_points[1][point], ellipses_points[2][point]])

                # interpolate X, Y, Z
                x = CubicSpline(t, points[:, 0])(tii)
                y = CubicSpline(t, points[:, 1])(tii)
                z = CubicSpline(t, points[:, 2])(tii)
                interpolated_points.append(np.array([x, y, z]))

            # construct ellipse from the four interpolated points
            ellipse = spatial_ellipse_from_4_points(interpolated_points[0], interpolated_points[1],
                                                    interpolated_points[2], interpolated_points[3])
            ellipses_smooth.append(ellipse)

        surface = np.array(ellipses_smooth)

        # ------------------------------------------------------------
        # create smooth surface
        points = vtk.vtkPoints()
        polys = vtk.vtkCellArray()

        numberOfLayers, numOfPoints, _ = surface.shape

        # add all surface points
        for layer in range(numberOfLayers):
            for point in range(numOfPoints):
                points.InsertNextPoint(surface[layer, point])

        # connect consecutive layers
        for layer in range(numberOfLayers - 1):

            for point in range(numOfPoints):

                nextPoint = (point + 1) % numOfPoints

                i0 = layer * numOfPoints + point
                i1 = layer * numOfPoints + nextPoint
                i2 = (layer + 1) * numOfPoints + point
                i3 = (layer + 1) * numOfPoints + nextPoint

                # triangle 1
                polys.InsertNextCell(3)
                polys.InsertCellPoint(i0)
                polys.InsertCellPoint(i2)
                polys.InsertCellPoint(i1)

                # triangle 2
                polys.InsertNextCell(3)
                polys.InsertCellPoint(i1)
                polys.InsertCellPoint(i2)
                polys.InsertCellPoint(i3)

        polydata = vtk.vtkPolyData()
        polydata.SetPoints(points)
        polydata.SetPolys(polys)

        # ------------------------------------------------------------
        # add smooth surface to MRML screen
        modelNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLModelNode", "Birth canal")
        modelNode.SetAndObservePolyData(polydata)
        displayNode = modelNode.GetDisplayNode()

        if displayNode is None:
            displayNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLModelDisplayNode")
            modelNode.SetAndObserveDisplayNodeID(displayNode.GetID())

        displayNode.SetColor(0.3, 0.7, 1.0)
        displayNode.SetOpacity(0.7)
        displayNode.SetBackfaceCulling(False)

        # ------------------------------------------------------------
        # pelvic outlet triangle
        points = vtk.vtkPoints()
        points.InsertNextPoint(float(pOutletA[0]), float(pOutletA[1]), float(pOutletA[2]))
        points.InsertNextPoint(float(pOutletB[0]), float(pOutletB[1]), float(pOutletB[2]) )
        points.InsertNextPoint(float(pOutletC[0]), float(pOutletC[1]), float(pOutletC[2]))

        # triangle edges
        lines = vtk.vtkCellArray()

        for i, j in [(0, 1), (1, 2), (2, 0)]:

            line = vtk.vtkLine()
            line.GetPointIds().SetId(0, i)
            line.GetPointIds().SetId(1, j)

            lines.InsertNextCell(line)

        # filled triangle
        polys = vtk.vtkCellArray()

        triangle = vtk.vtkPolygon()
        triangle.GetPointIds().SetNumberOfIds(3)

        triangle.GetPointIds().SetId(0, 0)
        triangle.GetPointIds().SetId(1, 1)
        triangle.GetPointIds().SetId(2, 2)

        polys.InsertNextCell(triangle)

        # polyData
        polydata = vtk.vtkPolyData()
        polydata.SetPoints(points)
        polydata.SetLines(lines)
        polydata.SetPolys(polys)

        # view
        outletNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLModelNode", "Pelvic outlet")
        outletNode.CreateDefaultDisplayNodes()
        outletNode.SetAndObservePolyData(polydata)

        displayNode = outletNode.GetDisplayNode()

        displayNode.SetColor(1.0, 0.0, 0.0)
        displayNode.SetOpacity(0.5)
        displayNode.SetLineWidth(3.0)
        displayNode.SetBackfaceCulling(False)

        outletNode.SetDisplayVisibility(True)

    else:

        qt.QMessageBox.warning(slicer.util.mainWindow(), QMessageTitle, 
                               'At least 1 landmark for defining birth canal is missing.')            

    # ------------------------------------------------------------
    # curve of carus
    noOfModelNode = slicer.mrmlScene.GetNumberOfNodesByClass('vtkMRMLModelNode')
    if noOfModelNode > 0:

        # apply to last model (head)
        modelNodeHead = slicer.mrmlScene.GetNthNodeByClass(noOfModelNode - 1, 'vtkMRMLModelNode')
        if modelNodeHead.GetName().find('Model') > -1: # pelvic floow (-1) + head (0)

            # new transform node
            #transformNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLLinearTransformNode", "MotionTransform")
            transformNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLLinearTransformNode", "StepTransform")

            # add transformation to model
            modelNodeHead.SetAndObserveTransformNodeID(transformNode.GetID())

            # model points
            volMeshModelNodeHead = modelNodeHead.GetMesh()
            pointsModelNodeHead = numpy_support.vtk_to_numpy(volMeshModelNodeHead.GetPoints().GetData())

            # head COG
            v9083001 = pointsModelNodeHead[7]
            COG = v9083001

            # mentovertical axis
            v803001 = pointsModelNodeHead[0] # point at chin (down) on mentovertical axis
            v803002 = pointsModelNodeHead[1] # point at top (up) on mentovertical axis
            MVD = v803001
            MVU = v803002
            mentovertical = MVU - MVD

            # parietal axis
            v803003 = pointsModelNodeHead[2] # left point on parietal axis
            v803004 = pointsModelNodeHead[3] # right point on parietal axis
            BPL = v803003
            BPR = v803004
            DPR = np.linalg.norm(BPR - BPL) / 2.0 # biparietal diameter
            biparietal = BPR - BPL # equal to asis L between P3 and P5

            v803005 = pointsModelNodeHead[4] # rear point on suboccipitobregmatic diameter
            v803006 = pointsModelNodeHead[5] # frontal point on suboccipitobregmatic diameter
            SOBR = v803005
            SOBF = v803006
            #suboccipit = SOBF - SOBR

            # curve of Carus
            modelNodePelvicFloor = slicer.mrmlScene.GetNthNodeByClass(noOfModelNode - 2, 'vtkMRMLModelNode')
            volMeshmodelNodePelvicFloor = modelNodePelvicFloor.GetMesh()
            pointsModelNodePelvicFloor = numpy_support.vtk_to_numpy(volMeshmodelNodePelvicFloor.GetPoints().GetData())
            
            # landmarks defining curve of Carus
            # PSA = pointsModelNodePelvicFloor[35484]
            # PSP = pointsModelNodePelvicFloor[35485]
            # PSI = pointsModelNodePelvicFloor[35473]
            # PR = pointsModelNodePelvicFloor[35487]
            # S3 = pointsModelNodePelvicFloor[35488]
            # S5 = pointsModelNodePelvicFloor[35490] # S5 = SJ

            # markupsNodes = slicer.util.getNodesByClass("vtkMRMLMarkupsNode")
            # allExist = 0 # count on all existing landmarks

            # for markupsNode in markupsNodes:

            #     # ------------------------------------------------------------
            #     # sequence of if is necessary, because PSI is twice,
            #     if markupsNode.GetName() == "PSA":
            #         PSA = [0, 0, 0]
            #         markupsNode.GetNthControlPointPositionWorld(0, PSA)
            #         allExist += 1
            
            #     elif markupsNode.GetName() == "PSP": # there is always single landmark in the node
            #         PSP = [0, 0, 0]
            #         markupsNode.GetNthControlPointPositionWorld(0, pInletA)
            #         allExist += 1

            #     elif markupsNode.GetName() == "PSI": # there is always single landmark in the node
            #         PSI = [0, 0, 0]
            #         markupsNode.GetNthControlPointPositionWorld(0, pInletA)
            #         allExist += 1

            #     elif markupsNode.GetName() == "PR": # there is always single landmark in the node
            #         PR = [0, 0, 0]
            #         markupsNode.GetNthControlPointPositionWorld(0, pInletA)
            #         allExist += 1

            #     elif markupsNode.GetName() == "S3": # there is always single landmark in the node
            #         S3 = [0, 0, 0]
            #         markupsNode.GetNthControlPointPositionWorld(0, pInletA)
            #         allExist += 1

            #     elif markupsNode.GetName() == "S5": # there is always single landmark in the node
            #         S5 = [0, 0, 0]
            #         markupsNode.GetNthControlPointPositionWorld(0, pInletA)
            #         allExist += 1

            # # all landmarks defining pelvic planes exist
            # if allExist == 6: # (6 landmarks)
            #     pass
            # else:
            #     qt.QMessageBox.warning(slicer.util.mainWindow(), QMessageTitle, 
            #                            'At least 1 landmark for defining birth canal is missing.')            

            PSA = pInletC
            PSP = pGreatestC
            PSI = pLeastC
            PR = pInletD
            S3 = pGreatestD
            S5 = pLeastD

            # ----------------------------------------------------------------------------------------------------------------------------
            P1 = (PSA + PR) / 2.0 # middle point of pelvic inlet (from the proximal end of pubic symphysis to proximal end of sacrum)
            P2 = (PSP + S3) / 2.0 # middle point of the midpelvic cavity (from middle of pubic symphysis to middle of sacrum (third sacral vertebrae))
            P3 = (PSI + S5) / 2.0 # middle point of pelvic outlet (from distal end of pubic symphysis to distal end of sacrum)

            # coordinates of P4 (closest point to pubic symphysis on plane of pubic arch that allowed the passage of fetal head)
            A = PSI
            B = pointsModelNodePelvicFloor[7]
            C = pointsModelNodePelvicFloor[6]

            E = (B - A) + (C - A) # auxiliary point E lies on line AP4
            alpha = np.arcsin(np.linalg.norm(np.cross((B - A), (E - A))) / (np.linalg.norm(B - A) * np.linalg.norm(E - A)))
            dP4A = DPR / np.sin(alpha)
            p = 0.95; # p is the approximation constant to symphysis
            P4 = (E - A) / np.linalg.norm(E - A) * dP4A * p

            # coordinates of P5 (point where fetal head is fully delivered)
            dP5 = 60.0 # estimate of sufficient distance between P4 and P5
            #k = 1.5 # k tilts P5 back to pelvis
            #P5 = np.array([0.0, k * dP5 + P4[2] - k * P4[0], -dP5])
            P5 = P4 + np.array([0.0, dP5, -dP5])

            # points defining curve of Carus
            P = np.array([P1, P2, P3, P4, P5])

            # display points defining curve of Carus
            markups = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLMarkupsFiducialNode")
            for i, p in enumerate(P):
                markups.AddControlPoint(*p, f"P{i+1}")

            # number of points along curve of Carus
            numberOfPoints = 20
            t = np.linspace(0, 1, numberOfPoints)
            
            # parametric B-spline (tck = tuple with parameters)
            tck, Pi = splprep([P[:, 1], P[:, 2]], s=0) # s=0 means accurately through points
            Pi = np.round(Pi * numberOfPoints)
            yi, zi = splev(t, tck)
            curveOfCarus = np.column_stack((np.zeros_like(yi), yi, zi))

            # ----------------------------------------------------------------------------------------------------------------------------
            # display curve of Carus
            points = vtk.vtkPoints()
            lines = vtk.vtkCellArray()

            for p in curveOfCarus:
                points.InsertNextPoint(p)

            lines.InsertNextCell(len(curveOfCarus))
            for i in range(len(curveOfCarus)):
                lines.InsertCellPoint(i)

            polyData = vtk.vtkPolyData()
            polyData.SetPoints(points)
            polyData.SetLines(lines)

            curve = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLModelNode", "Curve of Carus")
            curve.SetAndObservePolyData(polyData)

            curve.CreateDefaultDisplayNodes()
            displayNode = curve.GetDisplayNode()
            displayNode.SetColor(0, 1, 0)
            displayNode.SetLineWidth(1.0)
            # ----------------------------------------------------------------------------------------------------------------------------

            angleMV13 = np.radians(90.0) # mentovertical angle between P1 and P3
            angleBPextension12 = np.arcsin(np.linalg.norm(np.cross(PR - PSA, S3 - PSP)) / (np.linalg.norm(PR - PSA) * np.linalg.norm(S3 - PSP)))
            angleBPextension23 = np.arcsin(np.linalg.norm(np.cross(S3 - PSP, S5 - PSI)) / (np.linalg.norm(S3 - PSP) * np.linalg.norm(S5 - PSI)))
            angleBPextension34 = np.arcsin(np.linalg.norm(np.cross(S5 - PSI, P4 - PSI)) / (np.linalg.norm(S5 - PSI) * np.linalg.norm(P4 - PSI)))
            angleBP45 = np.radians(70.0) # final angle

            transform = vtkTransform() # transformation initialisation
            transform.Identity() # for cumulative transformation during loop

            T_global = np.eye(4) # initialisation of global transformation
            dtheta31 = angleMV13 / Pi[2] # mentovertical angle increment between P1 and P2
            #mentoverticalP1 = mentovertical # at P1

            # ------------------------------------------------------------------------------
            # def rotm2globXYZ(R):
            #     """
            #     Rozklad rotační matice na globální rotace:
            #     R = Rz(gamma) @ Ry(beta) @ Rx(alpha)

            #     Parametry:
            #         R : (3,3) numpy array
            #             rotační matice

            #     Návrat:
            #         alpha  ... rotace kolem globální osy X [rad]
            #         beta   ... rotace kolem globální osy Y [rad]
            #         gamma  ... rotace kolem globální osy Z [rad]
            #     """

            #     R = np.asarray(R, dtype=float)

            #     if R.shape != (3, 3):
            #         raise ValueError("R musí být 3x3 matice")

            #     # volitelná kontrola ortogonality
            #     if np.linalg.norm(R @ R.T - np.eye(3)) > 1e-6:
            #         print("Warning: R není přesně ortogonální")

            #     beta = np.arcsin(-R[2, 0])

            #     # ošetření gimbal locku
            #     if abs(np.cos(beta)) > 1e-6:
            #         alpha = np.arctan2(R[2, 1], R[2, 2])
            #         gamma = np.arctan2(R[1, 0], R[0, 0])
            #     else:
            #         # gimbal lock (beta = +-pi/2)
            #         alpha = 0.0
            #         gamma = np.arctan2(-R[0, 1], R[1, 1])

            #     return alpha, beta, gamma
            # ------------------------------------------------------------------------------

            # ------------------------------------------------------------------------------
            # dragon
            polyPoints = vtk.vtkPoints()
            polyPoints.InsertNextPoint(SOBR)
            polyPoints.InsertNextPoint(MVD)
            polyPoints.InsertNextPoint(SOBF)
            polyPoints.InsertNextPoint(MVU)

            polygon = vtk.vtkPolygon()
            polygon.GetPointIds().SetNumberOfIds(4)
            for i in range(4):
                polygon.GetPointIds().SetId(i, i)

            polys = vtk.vtkCellArray()
            polys.InsertNextCell(polygon)

            polyData = vtk.vtkPolyData()
            polyData.SetPoints(polyPoints)
            polyData.SetPolys(polys)

            polygonModel = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLModelNode", "MVBP_Polygon")
            polygonModel.SetAndObservePolyData(polyData)

            displayNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLModelDisplayNode")
            displayNode.SetColor(1, 0, 0)
            displayNode.SetOpacity(0.8)
            displayNode.SetBackfaceCulling(False)

            polygonModel.SetAndObserveDisplayNodeID(displayNode.GetID())
            # ------------------------------------------------------------------------------

            # ------------------------------------------------------------------------------
            # planes axes
            points = vtk.vtkPoints()
            points.InsertNextPoint(PSA.tolist())  # 0
            points.InsertNextPoint(PR.tolist())   # 1
            points.InsertNextPoint(PSP.tolist())  # 2
            points.InsertNextPoint(S3.tolist())   # 3
            points.InsertNextPoint(PSI.tolist())  # 4
            points.InsertNextPoint(S5.tolist())   # 5

            lines = vtk.vtkCellArray()

            def addLine(i, j):
                line = vtk.vtkLine()
                line.GetPointIds().SetId(0, i)
                line.GetPointIds().SetId(1, j)
                lines.InsertNextCell(line)

            addLine(0, 1)  # PSA - PR
            addLine(2, 3)  # PSP - S3
            addLine(4, 5)  # PSI - S5

            linePolyData = vtk.vtkPolyData()
            linePolyData.SetPoints(points)
            linePolyData.SetLines(lines)

            psLinesModel = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLModelNode", "PS_Lines")
            psLinesModel.SetAndObservePolyData(linePolyData)

            display = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLModelDisplayNode")
            display.SetColor(0, 0, 1)
            display.SetLineWidth(3)
            display.SetOpacity(1.0)
            display.SetBackfaceCulling(False)

            psLinesModel.SetAndObserveDisplayNodeID(display.GetID())
            # ------------------------------------------------------------------------------

            # ------------------------------------------------------------------------------   
            # biparietal line
            bpAxisPoints = vtk.vtkPoints()
            bpAxisPoints.InsertNextPoint(0, 0, 0)  # BPL
            bpAxisPoints.InsertNextPoint(0, 0, 0)  # BPR

            bpLine = vtk.vtkLine()
            bpLine.GetPointIds().SetId(0, 0)
            bpLine.GetPointIds().SetId(1, 1)

            bpLines = vtk.vtkCellArray()
            bpLines.InsertNextCell(bpLine)

            bpAxisPolyData = vtk.vtkPolyData()
            bpAxisPolyData.SetPoints(bpAxisPoints)
            bpAxisPolyData.SetLines(bpLines)

            bpAxisModel = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLModelNode", "Biparietal_Axis")
            bpAxisModel.SetAndObservePolyData(bpAxisPolyData)

            bpDisplay = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLModelDisplayNode")
            bpDisplay.SetColor(0, 1, 0)
            bpDisplay.SetLineWidth(5)
            bpDisplay.SetOpacity(1.0)
            bpDisplay.SetBackfaceCulling(False)

            bpAxisModel.SetAndObserveDisplayNodeID(bpDisplay.GetID())
            # ------------------------------------------------------------------------------

            # ------------------------------------------------------------------------------
            #  (P1, P2, P3, P4, P5) = (Pi[0], Pi[1], Pi[2], Pi[3], Pi[4])
            for tindex in range(1, numberOfPoints):

                # COG always on curve of Carus
                COG = curveOfCarus[tindex - 1, :3]

                if tindex <= Pi[1]: # around mentovertical and biparietal between P1 and P2
                    T21MV = rotationTransformMatrix(mentovertical, COG, dtheta31) # mentovertical = mentovertical at P1
                    extension21 = angleBPextension12 * tindex / (Pi[1] - Pi[0])
                    #T21BP = rotationTransformMatrix(biparietal, COG, extension21)
                    T21BP = rotationTransformMatrix([1.0, 0.0, 0.0], COG, extension21)
                    T_rotation = T21BP @ T21MV # increment rotation matrix

                elif tindex > Pi[1] and tindex <= Pi[2]: # around mentovertical and biparietal between P2 and P3
                    T32MV = rotationTransformMatrix(mentovertical, COG, dtheta31) # mentovertical rotates -> original mentovertical
                    extension31 = angleBPextension23 * tindex / (Pi[2] - Pi[1])
                    #T32BP = rotationTransformMatrix(biparietal, COG, extension31) # biparietal at P2
                    T32BP = rotationTransformMatrix([1.0, 0.0, 0.0], COG, extension31) # biparietal at P2
                    T_rotation = T32BP @ T32MV # increment rotation matrix

                elif tindex > Pi[2] and tindex <= Pi[3]: # around biparietal between P3 and P4
                    extension43 = angleBPextension34 / (Pi[3] - Pi[2]) # extension angle increment between P4 and P5
                    T43BP = rotationTransformMatrix(biparietal, COG, extension43) # extension matrix between P4 and P5
                    T_rotation = T43BP # increment rotation matrix

                else:  # around biparietal between P4 and P5
                    extension54 = angleBP45 / (Pi[4] - Pi[3]) # extension angle increment between P4 and P5
                    #T54BP = rotationTransformMatrix(biparietal, COG, extension54) # extension matrix between P4 and P5
                    T54BP = rotationTransformMatrix(biparietal, PSI, extension54) # extension matrix between P4 and P5
                    T_rotation = T54BP # increment rotation matrix

                # translation
                T_translation = np.eye(4)
                T_translation[:3, 3] = curveOfCarus[tindex, :3] - COG

                # increment transformation
                T_step = T_translation @ T_rotation # rotation followed by translation
                T_global = T_step @ T_global # accummulation to global transformation

                matrix = vtkMatrix4x4()
                for i in range(4):
                    for j in range(4):
                        matrix.SetElement(i, j, T_global[i, j])

                transformNode.SetMatrixTransformToParent(matrix)
                slicer.app.processEvents()

                #MVD = (T_global @ np.append(v803001, 1))[:3]
                #MVU = (T_global @ np.append(v803002, 1))[:3]
                MVD = (T_step @ np.append(MVD, 1))[:3]
                MVU = (T_step @ np.append(MVU, 1))[:3]
                mentovertical = MVU - MVD

                #BPL = (T_global @ np.append(v803003, 1))[:3]
                #BPR = (T_global @ np.append(v803004, 1))[:3]
                BPL = (T_step @ np.append(BPL, 1))[:3]
                BPR = (T_step @ np.append(BPR, 1))[:3]
                biparietal = BPR - BPL
                
                SOBF = (T_step @ np.append(SOBF, 1))[:3]
                SOBR = (T_step @ np.append(SOBR, 1))[:3]
                #suboccipit = SOBF - SOBR
                
                # ------------------------------------------------------------------------------
                # dragon
                polyPoints.SetPoint(0, SOBR.tolist())
                polyPoints.SetPoint(1, MVD.tolist())
                polyPoints.SetPoint(2, SOBF.tolist())
                polyPoints.SetPoint(3, MVU.tolist())

                polyPoints.Modified()
                polyData.Modified()
                polygonModel.Modified()
                # ------------------------------------------------------------------------------

                # ------------------------------------------------------------------------------
                # biparietal line
                bpAxisPoints.SetPoint(0, BPL.tolist())
                bpAxisPoints.SetPoint(1, BPR.tolist())

                bpAxisPoints.Modified()
                bpAxisPolyData.Modified()
                bpAxisModel.Modified()
                # ------------------------------------------------------------------------------

                time.sleep(0.5)

        else:
            qt.QMessageBox.warning(slicer.util.mainWindow(), 
                                   QMessageTitle, 'Model does not exist.')
            
    else:
        qt.QMessageBox.warning(slicer.util.mainWindow(), 
                                QMessageTitle, 'Model does not exist.')

    if pelvicFloorVersion == 'development': # check for version
        timePrint("... done", currentTime)

# ====================================================================================================
# mesh segmentation
def segmentationButton():
    """
    Segment model.
    
    Returns
    -------
    The segmented model displayed the main MRML screen.
    """

    if pelvicFloorVersion == 'development': # check for version
        currentTime = linePrint("Segmenting...")

    noOfModelNode = slicer.mrmlScene.GetNumberOfNodesByClass('vtkMRMLModelNode')
    if noOfModelNode > 0:

        # apply to last model
        modelNode = slicer.mrmlScene.GetNthNodeByClass(noOfModelNode - 1, 'vtkMRMLModelNode')
        if modelNode.GetName().find('Model') > -1:

            # delete current and create new segmentation
            deleteSegmentations() # delete current segmentation
            segmentationNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLSegmentationNode")
            slicer.modules.segmentations.logic().ImportModelToSegmentationNode(modelNode, segmentationNode)

            # # get closed surface representation of the segment
            # shellThickness = 3.0 # mm
            # segmentationNode = slicer.mrmlScene.GetNthNodeByClass(0, 'vtkMRMLSegmentationNode')
            # # segmentationNode.GetSegmentation().GetNumberOfSegments()
            # segmentationNode.CreateClosedSurfaceRepresentation()
            # polyData = segmentationNode.GetClosedSurfaceInternalRepresentation(modelNode.GetName())
            # #polyData = slicer.util.getNode(modelNode.GetName()).GetMesh()

            # seg = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLSegmentationNode")
            # seg.SetReferenceImageGeometryParameterFromVolumeNode(polyData)
            # seg.SetMasterRepresentationToClosedSurface()
            # slicer.modules.segmentations.logic().ImportModelToSegmentationNode(modelNode, seg)

            # # create shell
            # extrude = vtk.vtkLinearExtrusionFilter()
            # extrude.SetInputData(polyData)
            # extrude.SetExtrusionTypeToNormalExtrusion()
            # extrude.SetScaleFactor(shellThickness)

            # # compute consistent surface normals
            # triangle_filter = vtk.vtkTriangleFilter()
            # triangle_filter.SetInputConnection(extrude.GetOutputPort())
            # normals = vtk.vtkPolyDataNormals()
            # normals.SetInputConnection(triangle_filter.GetOutputPort())
            # normals.FlipNormalsOn()

            # # save result into new model node
            # slicer.modules.models.logic().AddModel(normals.GetOutputPort())

            view = slicer.app.layoutManager().threeDWidget(0).threeDView()
            view.cornerAnnotation().SetText(vtk.vtkCornerAnnotation.UpperRight,
                                            "Segmented model") # set text
            
            qt.QMessageBox.information(slicer.util.mainWindow(), 
                                       QMessageTitle, 'Model segmented.')
        else:
            qt.QMessageBox.warning(slicer.util.mainWindow(), 
                                   QMessageTitle, 'Nothing to segment.')
    else:
        qt.QMessageBox.warning(slicer.util.mainWindow(), 
                               QMessageTitle, 'Nothing to segment.')    

    if pelvicFloorVersion == 'development': # check for version
        timePrint("... done", currentTime)

# ====================================================================================================
# mesh labeled segmentation
def labelsButton():
    """
    Segment model by labeled slices.
    
    Returns
    -------
    The segmented model displayed the main MRML screen.
    """

    if pelvicFloorVersion == 'development': # check for version
        currentTime = linePrint("Labeling...")

    noOfModelNode = slicer.mrmlScene.GetNumberOfNodesByClass('vtkMRMLModelNode')
    if noOfModelNode > 0:

        # apply to last model
        modelNode = slicer.mrmlScene.GetNthNodeByClass(noOfModelNode - 1, 'vtkMRMLModelNode')
        if modelNode.GetName().find('Model') > -1:

            # try:
            #     slicer.app.pauseRender()
            #     ...
            # finally:
            #     slicer.app.resumeRender()
            slicer.app.pauseRender() # render just once at end
            segments, segmentVertices, segmentTriangles = separate(modelNode.GetMesh()) # separate mesh by lables
            deleteSegmentations() # delete current segmentations

            segmentationNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLSegmentationNode") # create new segmentation
            segmentationNode.CreateDefaultDisplayNodes()

            segmentation = segmentationNode.GetSegmentation()
            closedSurface = slicer.vtkSegmentationConverter.GetSegmentationClosedSurfaceRepresentationName()

            wasModified = segmentationNode.StartModify()

            for segment in segments:

                segmentMesh = createMesh(segmentVertices[segment], [], segmentTriangles[segment], [], [])

                extractSurface = vtk.vtkGeometryFilter()
                extractSurface.SetInputData(segmentMesh)
                extractSurface.Update()

                polydata = extractSurface.GetOutput()

                seg = slicer.vtkSegment()
                seg.SetName(str(segment))

                r, g, b = randomColor()
                seg.SetColor(r, g, b)

                seg.AddRepresentation(closedSurface, polydata)

                segmentation.AddSegment(seg)

            segmentationNode.EndModify(wasModified)
            slicer.app.resumeRender()

            view = slicer.app.layoutManager().threeDWidget(0).threeDView()
            view.cornerAnnotation().SetText(vtk.vtkCornerAnnotation.UpperRight,
                                            "Labeled model") # set text

            qt.QMessageBox.information(slicer.util.mainWindow(), 
                                       QMessageTitle, 'Model labeled.')
        else:
            qt.QMessageBox.warning(slicer.util.mainWindow(), 
                                   QMessageTitle, 'Nothing to be labeled.')
    else:
        qt.QMessageBox.warning(slicer.util.mainWindow(), 
                               QMessageTitle, 'Nothing to be labeled.')

    if pelvicFloorVersion == 'development': # check for version
        timePrint("... done", currentTime)

# ====================================================================================================
# delete all
def deleteButton():
    """
    Delete all markups, segmentations and models.
    
    Returns
    -------
    All markups, segmentations and models in the main MRML screen deleted.
    """

    if pelvicFloorVersion == 'development': # check for version
        currentTime = linePrint("Deleting models...")
    
    # renderer = slicer.app.layoutManager().threeDWidget(0).threeDView(). \
    #     renderWindow().GetRenderers().GetFirstRenderer()
    # renderer.RemoveAllViewProps()

    modelExists = deleteModels() # model nodes
    segmentationExists = deleteSegmentations() # segmentation nodes
    markupExists = deleteMarkups() # segmentation nodes

    if modelExists or segmentationExists or markupExists:
        qt.QMessageBox.information(slicer.util.mainWindow(), 
                                   QMessageTitle, 'All deleted.')
    else:
        qt.QMessageBox.warning(slicer.util.mainWindow(), 
                               QMessageTitle, 'Nothing to delete.')
    
    if pelvicFloorVersion == 'development': # check for version
        timePrint("... done", currentTime)

# ====================================================================================================
# delete markups
def deleteMarkupsButton():
    """
    Delete all models.
    
    Returns
    -------
    All markups in the main MRML screen deleted.
    """

    if pelvicFloorVersion == 'development': # check for version
        currentTime = linePrint("Deleting models...")
    
    markupExists = deleteMarkups() # segmentation nodes

    if markupExists:
        qt.QMessageBox.information(slicer.util.mainWindow(), 
                                   QMessageTitle, 'Markups deleted.')
    else:
        qt.QMessageBox.warning(slicer.util.mainWindow(), 
                               QMessageTitle, 'Nothing to delete.')
    
    if pelvicFloorVersion == 'development': # check for version
        timePrint("... done", currentTime)

# ====================================================================================================
# delete segmentations
def deleteSegmentationsButton():
    """
    Delete all segmentations.
    
    Returns
    -------
    All segmentations in the main MRML screen deleted.
    """

    if pelvicFloorVersion == 'development': # check for version
        currentTime = linePrint("Deleting models...")
    
    segmentationExists = deleteSegmentations() # segmentation nodes

    if segmentationExists:
        qt.QMessageBox.information(slicer.util.mainWindow(), 
                                   QMessageTitle, 'All deleted.')
    else:
        qt.QMessageBox.warning(slicer.util.mainWindow(), 
                               QMessageTitle, 'Nothing to delete.')
    
    if pelvicFloorVersion == 'development': # check for version
        timePrint("... done", currentTime)

# ====================================================================================================
# delete models
def deleteModelsButton():
    """
    Delete all models.
    
    Returns
    -------
    All models in the main MRML screen deleted.
    """

    if pelvicFloorVersion == 'development': # check for version
        currentTime = linePrint("Deleting models...")
    
    modelExists = deleteModels() # model nodes

    if modelExists:
        qt.QMessageBox.information(slicer.util.mainWindow(), 
                                   QMessageTitle, 'All deleted.')
    else:
        qt.QMessageBox.warning(slicer.util.mainWindow(), 
                               QMessageTitle, 'Nothing to delete.')
    
    if pelvicFloorVersion == 'development': # check for version
        timePrint("... done", currentTime)
