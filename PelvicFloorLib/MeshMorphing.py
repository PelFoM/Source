# ====================================================================================================
from __main__ import qt
from PelvicFloorLib import linePrint, timePrint, update, \
    pelvicFloorVersion, pathName, fileName, QMessageTitle, epsilon
import json, slicer, vtk, numpy as np
from vtk.util import numpy_support

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

    # MRI landmarks coordinates
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
# radial basis function
def RBF(rij, method):
    """
    Radial Basis Function mesh morphing.
    
    Returns
    -------
    phi : float
        The kernel value.
    """
    
    if method == 'linear':
        phi = rij
    elif method == 'quadratic':
        phi = (1 - rij) ** 2
    elif method == 'cubic':
        phi = rij ** 3
    elif method == 'thin plate spline':
        if rij > 0:
            phi = (rij ** 2) * np.log(rij)
        else:
            phi = 0.0
    elif method == 'Gaussian':
        phi = np.exp(-rij * rij / 100)
    else: # nothing happens
        phi = 0.0

    return(phi)

# ====================================================================================================
# mesh morphing
def morphButton():
    """
    Radial Basis Function mesh morphing.
    
    Returns
    -------
    Morphed mesh based on existing landmarks.
    """
    
    if pelvicFloorVersion == 'development': # check for version
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

            # any two landmarks must not be coincident and
            # at least four landmarks not lying in a plane
            if numberOfPairedLandmarks > 3 and \
                noOfCoincidentModelLandmarks == 0 and noOfCoincidentImageLandmarks == 0:

                # image landmarks must not be in plane
                modelX = [] # first three landmarks
                imageX = [] # first three landmarks
                for i in range(3):
                    # modelX.append(np.array(list(modelDictCoordinates.values())[landmark]))
                    # imageX.append(np.array(list(imageDictCoordinates.values())[landmark]))
                    modelX.append(np.r_[modelDictCoordinates[landmarks[i]]])
                    imageX.append(np.r_[imageDictCoordinates[landmarks[i]]])

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
                    alpha = 0 # morphing coefficient
                    method = 'thin plate spline' # radial basis function method ('linear', 'quadratic', 'cubic', 'thin plate spline', 'Gaussian')
                    #method = 'Gaussian'

                    volMesh = modelNode.GetMesh() # model coordinates (n, 3)
                    X = numpy_support.vtk_to_numpy(volMesh.GetPoints().GetData())
                    n = len(X) # number of nodes

                    # baseline L (m, 3) and target T (m, 3) landmarks coordinates
                    # np.array([[0] * m] * 3) -> slower
                    L = np.zeros((m, 3))
                    T = np.zeros((m, 3))

                    # modelDictCoordinates = imageDictCoordinates
                    for j, landmark in enumerate(modelDictCoordinates): # for all present  landmarks
                        L[j] = np.r_[modelDictCoordinates[landmark]] # baseline landmark L
                        T[j] = np.r_[imageDictCoordinates[landmark]] # target landmark T

                    # ----------------------------------------------------------------------------------------------------------------------------
                    # original_control_points (numpy.ndarray) – it is an (n_control_points, 3) array with the coordinates of the original interpolation control points
                    # deformed_control_points (numpy.ndarray) – it is an (n_control_points, 3) array with the coordinates of the interpolation control points after deformation
                    # func – the basis function to use in the transformation
                    # radius (float) – the scaling parameter r that affects the shape of the basis functions. For details see the class RBF. The default value is 0.5.
                    # extra_parameter (dict) – the additional parameters that may be passed to the kernel function. Default is None.                   

                    # from pygem import RBF
                    # rbf = RBF(func='gaussian_spline',original_control_points=L,deformed_control_points=T)
                    # new_mesh = rbf(X)
                    # ----------------------------------------------------------------------------------------------------------------------------

                    A = np.zeros((m, m)) # landmarks RBF matrix A (m, m)
                    for i in range(m): # for all landmarks
                        for j in range(m): # distance between landmarks
                            A[i][j] = RBF(np.linalg.norm(T[i] - L[j]), method)
                    
                    # baseline coordinates matrix [1, L]
                    B = np.concatenate((np.ones((m, 1)), L), axis=1)

                    # transformation matrix
                    M = np.concatenate(\
                        (np.concatenate((A + alpha * np.identity(m), B), axis=1), 
                         np.concatenate((np.transpose(B), np.zeros((4, 4))), axis=1)), 
                         axis=0)
                    R = np.concatenate((T, np.zeros((4, 3))), axis=0)
                    lambdac = np.linalg.solve(M, R)

                    A1 = np.zeros(m) # nodes RBF matrix A (n, m) by rows to save memory
                    targetPoints = vtk.vtkPoints() # target model node points

                    for i in range(n): # for all nodes
                        if i in modelDictNodes.values(): # node = baseline landmark
                            index = list(modelDictNodes.values()).index(i)
                            landmark = list(modelDictNodes.keys())[index]
                            X1 = imageDictCoordinates[landmark] # target node = target landmark
                        else:
                            for j in range(m): # for all landmarks
                                A1[j] = RBF(np.linalg.norm(X[i] - L[j]), method)
                            
                            B1 = np.r_[1, X[i]]
                            X1 = np.dot(np.append(A1, B1), lambdac)

                        targetPoints.InsertPoint(i, X1)

                    modelNode.GetMesh().SetPoints(targetPoints)
                    # modelNode = update(modelNode)
                    update(modelNode)

                    view = slicer.app.layoutManager().threeDWidget(0).threeDView()
                    view.cornerAnnotation().SetText(vtk.vtkCornerAnnotation.UpperRight,
                                                    "Morphed model") # set text

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
                if noOfCoincidentModelLandmarks > 0:
                    qt.QMessageBox.warning(slicer.util.mainWindow(), 
                                           QMessageTitle, 'There are ' + 
                                           str(noOfCoincidentModelLandmarks) + 
                                           ' coincident model landmarks.')
        else:
            qt.QMessageBox.warning(slicer.util.mainWindow(), 
                                   QMessageTitle, 'Nothing to morph.')    
    else:
        qt.QMessageBox.warning(slicer.util.mainWindow(), 
                               QMessageTitle, 'Nothing to morph.')    

    if pelvicFloorVersion == 'development': # check for version
        timePrint("... done", current_time)
