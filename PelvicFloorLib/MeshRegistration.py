# ====================================================================================================
# from __main__ import qt
# from PelvicFloorLib import linePrint, timePrint, update, \
#     pelvicFloorVersion, pathName, fileName, QMessageTitle, epsilon
# import json, slicer, vtk
import numpy as np
# from vtk.util import numpy_support

# ====================================================================================================
def rigid_cost(c, X, Y):
    """
    Cost function to find transformation to fit landmarks 'X' to landmarks 'Y' by least squares.
    
    Returns
    -------
    Cost function values.
    """

    return(np.linalg.norm(rigid(c, X) - Y))

# ====================================================================================================
# parameters for rigid registration (rotation and translation)
def rigid(c, X):
    """
    Function to transform landmarks 'X' to landmarks 'Y' by rotation matrix (angles c[0], c[1], c[2])
    and translation vector (shifts c[3], c[4], c[5]).
    
    Returns
    -------
    Transformed landmarks.
    """

    # rotation matrices
    Rx = [[1.0,           0.0,          0.0],
          [0.0,  np.cos(c[0]), np.sin(c[0])], 
          [0.0, -np.sin(c[0]), np.cos(c[0])]]
 
    Ry = [[ np.cos(c[1]), 0.0, np.sin(c[1])], 
          [          0.0, 1.0,          0.0], 
          [-np.sin(c[1]), 0.0, np.cos(c[1])]]
    
    Rz = [[ np.cos(c[2]), np.sin(c[2]), 0.0],
          [-np.sin(c[2]), np.cos(c[2]), 0.0],
          [          0.0,          0.0, 1.0]]

    R = np.dot(np.dot(Rz, Ry), Rx)
    # R = np.eye(3)
    # T = np.vstack([np.hstack([R, c[3 : 6].reshape(3, 1)]), np.array([0.0, 0.0, 0.0, 1.0])])
    # Y = np.dot(T, np.vstack([np.transpose(X), np.ones([1, np.size(X, 0)])]))

    # in case X is (n, 3) -> np.transpose(R * X) returns (n, 3) and c is added to each row
    # in case X is (1, 3) -> np.transpose(R * X) returns (1, 3) and c is added to row
    #return(X + c[3 : 6])

    return(np.transpose(np.dot(R, np.transpose(X))) + c[3 : 6])
    # return(np.transpose(Y[:-1]))

# ====================================================================================================
# radial basis function
def RBF(rij, type, radius):
    """
    Radial Basis Function mesh morphing.
    
    Returns
    -------
    phi : float
        The kernel value.
    """
    
    if type == 'linear':
        phi = rij
    elif type == 'quadratic':
        phi = (1 - rij) ** 2
    elif type == 'cubic':
        phi = rij ** 3
    elif type == 'Gaussian':
        phi = np.exp(-rij * rij / 100)
    elif type == 'thin plate spline':
        phi = np.zeros_like(rij, dtype=float)
        mask = rij > 0
        phi[mask] = rij[mask] ** 2 * np.log(rij[mask])
    elif type == 'Wendland':
        q = rij / radius
        phi = np.zeros_like(q)
        mask = q < 1
        qm = q[mask]
        phi[mask] = (1 - qm) ** 4 * (4 * qm + 1)
    else: # nothing happens
        phi = 0.0

    return(phi)

# ====================================================================================================
# Wendland C2 compact support RBF
# def RBFW(rij, radius):
#     """
#     Wendland C2 compact support Radial Basis Function mesh morphing.
#    
#     Returns
#     -------
#     phi : float
#         The kernel value.
#     """
#    
#     q = rij / radius
#     phi = np.zeros_like(q)
#
#     mask = q < 1
#     qm = q[mask]
#    
#     phi[mask] = (1 - qm) ** 4 * (4 * qm + 1)
#    
#     return(phi)
