
import sys
import os

# Import modules and functions

sys.path.append(os.getcwd())

from preprocessing.generate_case import default_settings
from postprocessing.routines import write_settings
from postprocessing.plot_contours import main

import numpy as np

def preprocess(casename, nsteps = 5000, cfl = 0.4, sfac = 0.5, meshf = 1.0):
    av = default_settings(casename)
    av['cfl'] = cfl;
    av['sfac'] = sfac;

    av['nsteps'] = nsteps
    av['ni'] = np.round(meshf*53)
    av['nj'] = np.round(meshf*37)

    write_settings(av)



preprocess('bump')
res = os.system('./build/cmdApp.exe')

