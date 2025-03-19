# ====================================================================================================
import numpy as np
import random, slicer, sys, time
from slicer import mrmlScene

# ----------------------------------------------------------------------------------------------------
# variables
# 'development' prints to CLI
# 'release' prints nothing
pelvicFloorVersion = 'development' # 'release'
QMessageTitle = 'Pelvic Floor Model'
opacity = 0.5
epsilon = 1e-3

# ----------------------------------------------------------------------------------------------------
# home directory in 'slicer.app.slicerHome/..'
# localDirectory = '/../Work/HBM/Models/PFD/Partners/UWA/Work/pelvic_floor/TemplateModel/'
localDirectory = '/../Work/HBM/Models/PFD/Morphing/slicer/pelvic_floor/TemplateModel/'
pathName = mrmlScene.GetRootDirectory() + localDirectory
# fileName = 'bony_pelvis'
# fileName = 'bony_pelvis_template'
fileName = 'bony_pelvis_template_single'
# fileName = 'bony_pelvis_template_single_cut'
# fileName = 'bony_pelvis_template_single_cut_2'
# fileName = '221221_Gynecoid'
# fileName = 'PelvicFloor'
# fileName = 'simpleBall'
extension = '.vtk'

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
    
    slicer.app.pythonConsole().clear()
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
