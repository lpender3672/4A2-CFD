import numpy as np
import os
import subprocess
from timeit import default_timer as timer
import matplotlib.pyplot as plt

from postprocessing.routines import (
    read_settings,
    write_settings,
    parse_output
)

global rel_folder
global outfname

rel_folder = 'cases/'

cases = [
    'bump',
    'bend'
]

# plot of non dimensionalized time per iteration against other parameters would be good
# its difficult to compare cases 

# the problem we have is that the solver doesnt recognize when it has converged properly
# and is wasting compute time on the last few iterations
# this means that python cannot actually measure the 
wdir = os.getcwd()
exe_path = wdir + '/build/solverApp.exe'
outfname = wdir + '/report/output.txt'

def timed_run_solver(casename):
    casepath = rel_folder + casename + '/input_' + casename + '.txt'
    start_time = timer()
    with open(outfname, 'w') as f:
        subprocess.run([exe_path, '--path', casepath], check=True, stdout=f)
    elapsed_time = timer() - start_time

    return elapsed_time

def apply_settings(casename, **kwargs):
    casepath = rel_folder + casename + '/input_' + casename + '.txt'
    av = read_settings(casepath)
    av['casename'] = casename
    
    valid_keys = ['cfl', 'sfac', 'nsteps', 'ni', 'nj']
    for k, v in kwargs.items():
        if k not in valid_keys:
            raise ValueError(f'Invalid key {k} provided. Valid keys are {valid_keys}')
        av[k] = v
    
    for key,val in av.items():
        if isinstance(val, list):
            av[key] = val[0]

    write_settings(av)

default_times = {}
cfls = [0.1, 0.2, 0.3, 0.4, 0.5]

for casename in cases:
    # loop to change a variable
    fig, ax = plt.subplots()
    times = []

    for cfl in cfls:
        apply_settings(casename, cfl=cfl)
        elapsed_time = timed_run_solver(casename)

        # read output
        with open(outfname, 'r') as f:
            lines = f.readlines()
        
        result = parse_output(lines)
        
        times.append(elapsed_time)
    
    ax.plot(cfls, times, label=casename)


plt.show()