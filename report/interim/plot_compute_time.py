import numpy as np
import os
import subprocess
from timeit import default_timer as timer
import matplotlib.pyplot as plt
import concurrent.futures

from postprocessing.routines import (
    default_settings,
    read_settings,
    write_settings,
    parse_run_output,
    check_run_converged,
    read_conv
)

global rel_folder
global outfname
global appver

rel_folder = 'cases/'
appver = '1.0'

cases = [
    'bump',
    #'bend'
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

    return av


def cfl_sfac_grid_run(casename):
    
    try:
        history = np.loadtxt(f'report/data/{casename}_runs_{appver}.txt')
    except FileNotFoundError:
        history = np.zeros((0, 9))

    cfls = np.logspace(-2, np.log10(0.5), 10, endpoint = True)
    sfacs = np.arange(0.05, 0.65, 0.05)
    ni = 53
    nj = 37

    print(f'Running {casename} with {len(cfls) * len(sfacs)} combinations')
    print("cfls:", cfls)
    print("sfacs:", sfacs)
    

    for cfl in cfls:
        for sfac in sfacs:
            # we want 3 unique runs for each cfl and sfac to average
            # decided to just do this in another nested loop
            duplicates =  np.sum((history[:, 0] == cfl) * (history[:, 1] == sfac) * (history[:, 7] == ni))
            if duplicates >= 1:
                continue

            av = apply_settings(casename, cfl=cfl, sfac=sfac, ni=ni, nj=nj)
            time = 0
            n = 3
            for _ in range(n):
                time += timed_run_solver(casename)
            time /= n

            with open(outfname, 'r') as f:
                lines = f.readlines()
            
            convpath = rel_folder + casename + '/conv_' + casename + '.csv'
            conv_hist = read_conv(convpath)
            
            if conv_hist['nstep'].shape[0] > 0:
                iterations = conv_hist['nstep'][-1]
                dro_max = conv_hist['dro_max'][-1]
                dro_avg = conv_hist['dro_avg'][-1]
            else:
                iterations = 1
                dro_max = np.inf
                dro_avg = np.inf

            converged = check_run_converged(lines)

            newrow = np.array([
                cfl, sfac, time, converged, iterations, dro_max, dro_avg, av['ni'], av['nj']
            ])
            history = np.vstack((history, newrow))

            print(f'cfl: {cfl}, sfac: {sfac}, time: {time}, converged: {converged}, iterations: {iterations}, dro_max: {dro_max}, dro_avg: {dro_avg}')

    np.savetxt(f'report/data/{casename}_runs_{appver}.txt', history)


def plot_cfl_sfac_scatter(casename):

    cfl_sfac_grid_run(casename)
    
    fig, ax = plt.subplots(figsize=[12,9])

    fig.subplots_adjust(
        top=0.872,
        bottom=0.087,
        left=0.073,
        right=0.983,
        hspace=0.2,
        wspace=0.2
    )

    try:
        history = np.loadtxt(f'report/data/{casename}_runs_{appver}.txt')
    except FileNotFoundError:
        return ax

    # remove rows where the run did not converge
    converged_within = history[history[:, 3] == 0]
    converged_outside = history[history[:, 3] == 1]
    diverged = history[history[:, 3] == 2]
    maxiter = history[history[:, 3] == 3]
    
    #cols = np.log10(passed_runs[:, 2] / np.min(passed_runs[:, 2]))
    if converged_within.shape[0] > 0:
        colsin = converged_within[:, 2] / np.min(converged_within[:, 2])
        #colsin = np.log10(colsin)
    else:
        colsin = []
    conin = ax.scatter(converged_within[:, 0], 
                        converged_within[:, 1], 
                        c=colsin, 
                        marker='o', 
                        label='Converged within bounds',
                        cmap='jet',
                        s=100)

    if converged_outside.shape[0] > 0:
        colsout = converged_outside[:, 2] / np.min(converged_outside[:, 2])
        #colsout = np.log10(colsout)
    else:
        colsout = []
    conout = ax.scatter(converged_outside[:, 0], 
                        converged_outside[:, 1], 
                        c=colsout, 
                        marker='d', 
                        label='Converged outside bounds',
                        cmap='jet',
                        s=100)

    div = ax.scatter(diverged[:, 0], 
                        diverged[:, 1], 
                        c='r', 
                        marker='x', 
                        label='Diverged',
                        s=100)
    
    maxit = ax.scatter(maxiter[:, 0],
                        maxiter[:, 1],
                        c='k',
                        marker='s',
                        label='Max iterations reached',
                        s=100)
    # add failed runs to plot
    
    ax.set_xlabel('CFL')
    ax.set_ylabel('SFAC')

    # add colourbar
    cbar = plt.colorbar(conout)
    cbar.set_label('Minimum normalised run time')


    ax.grid(which='both', linestyle='--')
    ax.set_xscale('log')
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15),
          ncol=2, fancybox=True, shadow=True)
    
    return ax


def d_avg_cfl_run(casename):

    sfac = 0.5 # higher sfac allows larger range of cfls
    cfls = np.logspace(-2, np.log10(0.5), 10, endpoint = True)
    ni = 53
    nj = 37

    try:
        history = np.loadtxt(f'report/data/{casename}_runs_{appver}.txt')
    except FileNotFoundError:
        history = np.zeros((0, 9))

    for cfl in cfls:
        duplicates =  np.sum((history[:, 0] == cfl) * (history[:, 1] == sfac) * (history[:, 7] == ni))
        if duplicates >= 1:
            continue

        av = apply_settings(casename, cfl=cfl, sfac=sfac, ni=ni, nj=nj)
        time = 0
        n = 3
        for _ in range(n):
            time += timed_run_solver(casename)
        time /= n

        with open(outfname, 'r') as f:
            lines = f.readlines()
        
        convpath = rel_folder + casename + '/conv_' + casename + '.csv'
        conv_hist = read_conv(convpath)
        
        if conv_hist['nstep'].shape[0] > 0:
            iterations = conv_hist['nstep'][-1]
            dro_max = conv_hist['dro_max'][-1]
            dro_avg = conv_hist['dro_avg'][-1]
        else:
            iterations = 1
            dro_max = np.inf
            dro_avg = np.inf

        converged = check_run_converged(lines)

        newrow = np.array([
            cfl, sfac, time, converged, iterations, dro_max, dro_avg, av['ni'], av['nj']
        ])
        history = np.vstack((history, newrow))

        print(f'cfl: {cfl}, sfac: {sfac}, time: {time}, converged: {converged}, iterations: {iterations}, dro_max: {dro_max}, dro_avg: {dro_avg}')

    np.savetxt(f'report/data/{casename}_runs_{appver}.txt', history)


def d_avg_ni_run(casename):

    sfac = 0.5 # higher sfac allows larger range of cfls
    cfl = 0.13572088082974532
    nis = np.logspace(1, 3, 10, endpoint = True).astype(int)
    nj = 37

    try:
        history = np.loadtxt(f'report/data/{casename}_runs_{appver}.txt')
    except FileNotFoundError:
        history = np.zeros((0, 9))

    for ni in nis:
        duplicates =  np.sum((history[:, 0] == cfl) * (history[:, 1] == sfac) * (history[:, 7] == ni))
        if duplicates >= 1:
            continue

        av = apply_settings(casename, cfl=cfl, sfac=sfac, ni=ni, nj=nj)
        time = 0
        n = 3
        for _ in range(n):
            time += timed_run_solver(casename)
        time /= n

        with open(outfname, 'r') as f:
            lines = f.readlines()
        
        convpath = rel_folder + casename + '/conv_' + casename + '.csv'
        conv_hist = read_conv(convpath)

        if conv_hist['nstep'].shape[0] > 0:
            iterations = conv_hist['nstep'][-1]
            dro_max = conv_hist['dro_max'][-1]
            dro_avg = conv_hist['dro_avg'][-1]
        else:
            iterations = 1
            d_max = np.inf
            d_avg = np.inf

        converged = check_run_converged(lines)

        newrow = np.array([
            cfl, sfac, time, converged, iterations, dro_max, dro_avg, ni, nj
        ])
        history = np.vstack((history, newrow))

        print(f'cfl: {cfl}, sfac: {sfac}, time: {time}, converged: {converged}, iterations: {iterations}, dro_max: {d_max}, dro_avg: {d_avg}')

    np.savetxt(f'report/data/{casename}_runs_{appver}.txt', history)

def reshape_saved_data(casename):
    
    history = np.loadtxt(f'report/data/{casename}_runs_{appver}.txt')
    # add two columns to end with ni and nj
    ni = 53
    nj = 37
    newdata = np.zeros((history.shape[0], 9))
    newdata[:, :-2] = history
    newdata[:, -2] = ni
    newdata[:, -1] = nj

    np.savetxt(f'report/data/{casename}_runs_{appver}.txt', newdata)

def plot_d_avg_cfl(casename):

    fig,ax = plt.subplots()

    try:
        history = np.loadtxt(f'report/data/{casename}_runs_{appver}.txt')
    except FileNotFoundError:
        return ax
    
    sfac = 0.5 # higher sfac allows larger range of cfls
    
    # filter by sfac and cfl
    history = history[(history[:, 1] == sfac)]
    # sort by cfl
    history = history[history[:, 0].argsort()]


    x = np.linspace(np.min(history[:, 0]), np.max(history[:, 0]), 100)
    ax.loglog(x, 1e-4 * x**1, 'k--', label=r'$O(CFL^1)$')

    ax.loglog(history[:, 0], history[:, 6], 'o-', label = 'Average density residual error')

    ax.set_xlabel('CFL')
    ax.set_ylabel('Average density residual error')

    ax.grid(which='both', linestyle='--')
    ax.legend()

    return ax


if __name__ == '__main__':
    
    
    d_avg_cfl_run('bump')
    plot_d_avg_cfl('bump')

    plt.tight_layout()
    plt.show()