# ====================================================================================================
from __main__ import qt
from PelvicFloorLib import linePrint, timePrint, randomColor, \
    createMesh, separate, getLandmarks, \
        read, export, align, move, rotate, scale, update, \
            deleteModels, deleteSegmentations, deleteMarkups, \
                pelvicFloorVersion, pathName, fileName, extension, QMessageTitle, opacity
import json, slicer, vtk, numpy as np
from vtk.util import numpy_support
from DICOMLib import DICOMUtils

# ====================================================================================================
# open VTK
def openButton():

    # open dialog
    fileDialog = qt.QFileDialog()
    inputFileName = fileDialog.getOpenFileName(fileDialog, '', pathName, '*.vtk')

    # read model
    modelNode = read(inputFileName, 'Template model', 'vtp')
    volMesh = modelNode.GetMesh()
    points = numpy_support.vtk_to_numpy(volMesh.GetPoints().GetData())
    
    # show template model
    update(modelNode)

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

    # # read file (https://www.tutorialspoint.com/pyqt/pyqt_qfiledialog_widget.htm)
    # fileDialog = qt.QFileDialog()
    # inputFileName = fileDialog.getOpenFileName(fileDialog, '', pathName, '*.vtk; *.vtp; *.vtu')
    
    # # default file name
    # if inputFileName == '': # no file chosen
    #     pathFileName = pathName + fileName + extension
    # else:
    #     pathFileName = inputFileName
    pathFileName = pathName + fileName + extension

    # read model
    modelNode = read(pathFileName, 'Template model', 'vtu')

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
    
    # convert to LPS (internally RAS) coordinate system
    points[:, [0, 1, 2]] = -points[:, [2, 0, 1]]
    
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
def exportButton():
    """
    Export model to VTK as unstructured grid.
    
    Returns
    -------
    The model exported to VTK.
    """

    if pelvicFloorVersion == 'development': # check for version
        currentTime = linePrint("Exporting model...")

    noOfModelNode = slicer.mrmlScene.GetNumberOfNodesByClass('vtkMRMLModelNode')
    if noOfModelNode > 0:
        
        # apply to last model only
        modelNode = slicer.mrmlScene.GetNthNodeByClass(noOfModelNode - 1, 'vtkMRMLModelNode')
        if modelNode.GetName().find('Model') == 0: # last name 'Model*' exist

            # save file
            fileDialog = qt.QFileDialog()
            outputFileName = fileDialog.getSaveFileName(fileDialog, '', pathName, '*.vtk')

            if not outputFileName == '': # file name chosen
                
                # LPS -> RAS
                points = numpy_support.vtk_to_numpy(modelNode.GetMesh().GetPoints().GetData())
                points[:, 1] = -points[:, 1]
                # update(modelNode)

                # file exported in RAS coordinate system
                # #modelWriter = vtk.vtkUnstructuredGridWriter()
                modelWriter = vtk.vtkPolyDataWriter()
                modelWriter.SetFileName(outputFileName)
                modelWriter.SetInputData(modelNode.GetMesh())
                #modelWriter.SetFileTypeToBinary()
                #modelWriter.Update()
                modelWriter.Write()

                # export to PC file (RAS coordinate system)
                export(outputFileName[ : -4], modelNode.GetMesh())

                # RAS -> LPS back for MRML
                points[:, 1] = -points[:, 1]
                # update(modelNode)

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

    try: # read merged markups
        if fileType == '.csv': # CSV (including labels)
            markupsNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLMarkupsFiducialNode")
            slicer.modules.markups.logic().ImportControlPointsFromCSV(markupsNode, inputFileName)
        elif fileType == 'json': # JSON (not labeled)
            slicer.util.loadMarkups(inputFileName)
            # landmarks = ('PSI', 'PSP', 'PR', 'IL_L', 'IL_R')
        elif fileType == 'fcsv': # FCSV (not labeled)
            slicer.util.loadMarkupsFiducialList(inputFileName)
            landmarks = ('PSA',       # index 1 (pubic syphysis anterior, PSA)
                         'PSP',       # index 2 (pubic syphysis posterior, PSP)
                         'PSI',       # index 3 (pubic syphysis inferior, PSI)
                         'L5',        # index 4 (L5 spinous process)
                         'PR',        # index 5 (sacral promontory)
                         'S3',        # index 6 (S3 inferior surface)
                         'S45',       # index 7 (border between S4 and S5)
                         'S5',        # index 8 (S5 inferior surface)
                         'S6',        # index 9 (S6 inferior surface)
                         'IT_L',      # index 10 (left ischial tuberosity)
                         'IT_R',      # index 11 (right ischial tuberosity)
                         'IS_L',      # index 12 (left ischial spine)
                         'IS_R',      # index 13 (right ischial spine)
                         'BS_L',      # index 14 (left iliac spine)
                         'BS_R',      # index 15 (right iliac spine)
                         'BC_L',      # index 16 (left iliac crest)
                         'BC_R',      # index 17 (right iliac crest)
                         'TR_L',      # index 18 (left trochanter)
                         'TR_R',      # index 19 (right trochanter)
                         'IL_L',      # index 20 (left transverse diameter)
                         'IL_R',      # index 21 (right transverse diameter)
                         'PSI_t0',    # index 22 (PSI at t0)
                         'S5_t0',     # index 23 (S5 at t0)
                         'Cervix_t0', # index 24 (cervix at t0)
                         'PSI_t1',    # index 25 (PSI at t1)
                         'S5_t1',     # index 26 (S5 at t1)
                         'Cervix_t1', # index 27 (cervix at t1)
                         'age',       # index 28 (years)
                         'height',    # index 29 (cm)
                         'weight',    # index 30 (kg)
                         'SB',        # index 31 (back point on SJ of S6)
                         'TI_L',      # index 32 (left tuber ischiadicum)
                         'TI_R',      # index 33 (right tuber ischiadicum)
                         'PSIS_L',    # index 34 (left PSIS)
                         'PSIS_R',    # index 35 (right PSIS)
                         'SA_L',      # index 36 (left sacrum)
                         'SA_R',      # index 37 (right sacrum)
                         'SB2',       # index 38 (back point on SJ of S5)
                         'Aa',        # index 39 (at midline of anterior vaginal wall)
                         'Ba',        # index 40 (most superior location of front vaginal wall)
                         'C',         # index 41 (lowest edge of cervix or the vaginal cuff)
                         'D',         # index 42 (topmost point of posterior vaginal wall)
                         'Ap',        # index 43 (midline of posterior vaginal wall)
                         'Bp',        # index 44 (uppermost point of posterior vaginal wall)
                         'E',         # index 45 (additional point to measure GH (genital hiatus))
                         'F',         # index 46 (additional point to measure GH/PB)
                         'G');        # index 47 (additional point to measure PB (perineal body))
        else: # dictionary (including labels)
            with open(inputFileName) as landmarks:
                landmarksDictXyz = json.loads(landmarks.read())

            for landmark in landmarksDictXyz: # separate markups
                node = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLMarkupsFiducialNode")
                node.SetName(landmark)
                slicer.util.updateMarkupsControlPointsFromArray(node, 
                                                np.array([landmarksDictXyz[landmark]]))
                
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
    If landmark PSI is not present, nothing happens and warning is printed.
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
            deleteSegmentations()
            segmentationNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLSegmentationNode")
            segmentationNode.CreateDefaultDisplayNodes() # only needed for display

            # import model into segmentation node
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

            # separate mesh by lables
            segments, segmentVertices, segmentTriangles = separate(modelNode.GetMesh())

            # delete old and create new segmentation
            deleteSegmentations()
            segmentationNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLSegmentationNode")
            segmentationNode.CreateDefaultDisplayNodes() # only needed for display

            for segment in segments: # over all segments in cell data
                segmentMesh = createMesh(segmentVertices[segment], [], 
                                         segmentTriangles[segment], [], [])
                
                extractSurface = vtk.vtkGeometryFilter()
                extractSurface.SetInputData(segmentMesh)
                # subModelNode.SetAndObservePolyData(segmentMesh)
                extractSurface.Update()

                segmentModelNode = slicer.modules.models.logic().AddModel(extractSurface.GetOutput())
                segmentModelNode.SetDisplayVisibility(0) # segmentModelNode.GetDisplayNode().SetVisibility(False)
                segmentModelNode.SetName(str(segment))
                r, g, b = randomColor()
                segmentModelNode.GetDisplayNode().SetColor(r, g, b)
                segmentModelNode.GetDisplayNode().SetOpacity(opacity)

                # add to segmentation
                slicer.modules.segmentations.logic().ImportModelToSegmentationNode(segmentModelNode, 
                                                                                segmentationNode)
                slicer.mrmlScene.RemoveNode(segmentModelNode)

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
# mesh scaling
def scaleButton():
    """
    Scale in two dimensions.
    
    Returns
    -------
    Scaled mesh based on existing landmarks.
    """

    if pelvicFloorVersion == 'development': # check for version
        current_time = linePrint("Mesh scaling...")

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
        timePrint("... done", current_time)

# ====================================================================================================
# show birth canal
def canalButton():
    """
    Birth canal view.
    
    Returns
    -------
    The birth canal displayed in the main MRML screen.
    """

    if pelvicFloorVersion == 'development': # check for version
        currentTime = linePrint("Displaying birth canal...")

    noOfModelNode = slicer.mrmlScene.GetNumberOfNodesByClass('vtkMRMLModelNode')
    if noOfModelNode > 0:

        # apply to last model
        modelNode = slicer.mrmlScene.GetNthNodeByClass(noOfModelNode - 1, 'vtkMRMLModelNode')
        if modelNode.GetName().find('Model') > -1:

            # check existence of landmarks PR, PSP, S4, SJ, ILL, ILR, ISL, ISR, ITL, ITR
            landmarksDict = {'PSP': False, 'PR': False, 'S4': False, 'SJ': False, 
                             'IL_L': False, 'IL_R': False, 'IS_L': False, 'IS_R': False, 
                             'IT_L': False, 'IT_R': False}
            allExist, modelDictXyz, imageDictXyz = getLandmarks(modelNode, landmarksDict)

            # all necessary landmarks exist
            if allExist:

                # draw birth canal
                pass
                
                # update model node
                update(modelNode)

                qt.QMessageBox.information(slicer.util.mainWindow(), 
                                           QMessageTitle, 'Birth canal displayed.')
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
