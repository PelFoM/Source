# ====================================================================================================
import numpy as np
import random, sys, time
from slicer import mrmlScene
from pathlib import Path

# ----------------------------------------------------------------------------------------------------
# variables
# 'development' prints to CLI
# 'release' prints nothing
pelvicFloorVersion = 'development' # 'release'
QMessageTitle = 'Pelvic Floor Model'
opacity = 0.5
epsilon = 1e-3
alpha = 1e-6 # morphing coefficient
morphingType = 'Wendland' # radial basis function method ('linear', 'quadratic', 'cubic', 'thin plate spline', 'Gaussian')

# ----------------------------------------------------------------------------------------------------
# home directory in 'slicer.app.slicerHome/..'
PelvicFloorDirectory = Path(__file__).resolve().parent # this script directory

#localDirectory = '/../Work/HBM/Models/PFD/Morphing/slicer/pelvic_floor/TemplateModel/'
#pathName = mrmlScene.GetRootDirectory() + localDirectory # mrmlScene.GetRootDirectory() = C:/Users/hyncik/Documents

# tempate model relative path to this script directory
localDirectory = '../TemplateModel'
pathName = (PelvicFloorDirectory / localDirectory).resolve()

# template file name
fileName = 'template'
#fileName = 'templateModelPelvis'
extension = '.vtk'

#pathFileName = pathName + fileName + extension
pathFileName = str(pathName / f"{fileName}{extension}")

nodeText = mrmlScene.AddNewNodeByClass("vtkMRMLTextNode", "pathFileName")
nodeText.HideFromEditorsOn() # if not hidden => collides with other nodes
nodeText.SetText(pathFileName)
#nodeText = slicer.util.getNode("pathFileName")
#fileName = nodeText.GetText()

# ====================================================================================================
# print text and return actual time
def linePrint(text):
    """
    Print 'text' without EOL and continues.

    Returns
    -------
    The 'text' printed to stdout without EOL.
    Current time : float
        Current time is returend.
    """
    
    # starting with Slicer versions 5.5+ it is no longer possible to do this directly
    #slicer.app.pythonConsole().clear() # collapse in older versions
    print("\n" * 100)
    print(text, end = '')
    sys.stdout.flush()
    return time.time()

# ====================================================================================================
# print text and actual time and return actual time
def timePrint(text, initialTime):
    """
    Print 'text' at 'initialTime'.

    Returns
    -------
    Current time : float
        Current time is returend.
    """     
    
    print(text + " in %.3f seconds." % (time.time() - initialTime))
    return time.time()

# ====================================================================================================
# generate random color
# def randomColor():
#     levels = range(32,256,32)
#     return tuple(random.choice(levels) for _ in range(3))
def randomColor():
    """
    Return random color.

    Returns
    -------
    Random color : tuple
        Random color in (red, green, blue) is returned.
    """     

    return tuple(random.choice(np.linspace(0.0, 1.0, num=11)) for _ in range(3))

# ====================================================================================================
def rotationMatrix(axis, theta):
    """Rodrigues rotation matrix 3x3"""
    axis = np.array(axis, dtype=float)
    axis = axis / np.linalg.norm(axis)
    ux, uy, uz = axis
    c = np.cos(theta)
    s = np.sin(theta)
    C = 1 - c
    R = np.array([
        [c + ux**2*C,     ux*uy*C - uz*s,  ux*uz*C + uy*s],
        [uy*ux*C + uz*s,  c + uy**2*C,     uy*uz*C - ux*s],
        [uz*ux*C - uy*s,  uz*uy*C + ux*s,  c + uz**2*C]
    ])
    return R

# ====================================================================================================
def rotationTransformMatrix(axis, point, theta):
    """Rotation matrix 4x4 around axis at point"""
    R = rotationMatrix(axis, theta)
    point = np.array(point)
    t = point - R @ point
    T = np.eye(4)
    T[:3, :3] = R
    T[:3, 3] = t
    return T

# ====================================================================================================
def ellipseFromDiameters(p1, p2, p3, p4, numberOfPoints=100):
    center = (p1 + p2 + p3 + p4) / 4
    axis_a = (p1 - p2) / 2
    axis_b = (p3 - p4) / 2
    theta = np.linspace(0, 2 * np.pi, numberOfPoints)
    ellipse = np.array([center + np.cos(t) * axis_a + np.sin(t) * axis_b for t in theta])
    return ellipse
