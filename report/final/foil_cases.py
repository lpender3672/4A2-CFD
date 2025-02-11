
import os
from pathlib import Path
from postprocessing.routines import (
    read_settings,
    default_settings,
)
from runtime_enviroment import CFDManager

import numpy as np
import pandas as pd
import copy

import matplotlib.pyplot as plt
import mat73

import scipy.io as sio


# load the .mat file
#naca0012_sweep = sio.loadmat('Data/0012_swp.mat')

# print the contents of the .mat file
# print(naca0012_sweep.keys())


def create_cfd_env(avs, name):

    n_workers = 8
    base_folder = Path(os.getcwd()) / "case_env"
    data_file = Path(os.getcwd()) / "report" / "final" / "data" / name
    cfd_executable = 'build/solverApp.exe'
    #av = read_settings('cases/bump/input_bump.txt')
    manager = CFDManager(n_workers, base_folder, data_file, cfd_executable, avs)

    return manager

def plot_clcd_alpha():

    av_0012 = read_settings('cases/naca0012/input_naca0012.txt')
    av_2412 = read_settings('cases/naca2412/input_naca2412.txt')

    av_0012['cfl'] = 0.2
    av_2412['cfl'] = 0.2

    av_templates = [av_0012, av_2412]
    alphas = np.linspace(-10, 20, 13)
    print(alphas)

    avs = []
    for av in av_templates:
        for alpha in alphas:
            av['alpha'] = alpha
            avs.append(av.copy())

    manager = create_cfd_env(avs, 'cl_alpha.csv')
    manager.clear_worker_folders()
    manager.start_workers(1)

    lift_file = manager.shared_file.parent / "cl_alpha.csv"
    df = pd.read_csv(lift_file)

    # filter by converged
    df = df[df['converged'] < 2]
    # aggregate dt column, merging rows with the same alpha
    dt_mean = df.groupby('alpha')['dt'].mean().reset_index()
    df = df.merge(dt_mean, on='alpha', how='left')

    # separate naca0012 and naca2412
    df_naca0012 = df[df['casename'] == 'naca0012']
    df_naca2412 = df[df['casename'] == 'naca2412']

    # sort by alpha
    df_naca0012 = df_naca0012.sort_values('alpha')
    df_naca2412 = df_naca2412.sort_values('alpha')

    # load BEM data
    naca0012_sweep = sio.loadmat('report/final/data/0012_swp.mat')

    alpha0012 = naca0012_sweep['alpha'][0,:]
    cl0012 = naca0012_sweep['clswp'][0,:]
    cd0012 = naca0012_sweep['cdswp'][0,:]

    naca2412_sweep = sio.loadmat('report/final/data/2412_swp.mat')

    alpha2412 = naca2412_sweep['alpha'][0,:]
    cl2412 = naca2412_sweep['clswp'][0,:]
    cd2412 = naca2412_sweep['cdswp'][0,:]

    fig, ax = plt.subplots()

    ax.plot(df_naca0012['alpha'].to_numpy(), 
            df_naca0012['cl'].to_numpy(), 
            '-o',
            label='CFD NACA0012')
    ax.plot(df_naca2412['alpha'].to_numpy(), 
            df_naca2412['cl'].to_numpy(),
            '-o',
            label='CFD NACA2412')
    
    ax.plot(alpha0012, cl0012, label='PM NACA0012')
    ax.plot(alpha2412, cl2412, label='PM NACA2412')

    cl1 = df_naca2412['cl'][df_naca2412['alpha'] == 0]
    cl2 = cl2412[np.where(alpha2412 == 0)[0][0]]
    print(f'cl1 {cl1}, cl2 {cl2}, diff {np.abs(cl1 - cl2)/cl2}')

    ax.set_ylim([
        df_naca0012['cl'].min() - 0.2,
        df_naca2412['cl'].max() + 0.2
    ])

    ax.grid( linestyle='--', linewidth=0.5)
    ax.legend()

    ax.set_xlabel('Angle of attack [deg]')
    ax.set_ylabel('Lift coefficient')

    fig.tight_layout()
    fig.savefig('report/final/figures/cl_alpha.png', dpi=300)

    fig, ax = plt.subplots()


    ax.plot(df_naca0012['alpha'].to_numpy(),
            df_naca0012['cd'].to_numpy(),
            '-o',
            label='NACA0012')
    ax.plot(df_naca2412['alpha'].to_numpy(),
            df_naca2412['cd'].to_numpy(),
            '-o',
            label='NACA2412')
        
    theoretical_gradient = 2 / np.pi
    gradient_0012 =  np.diff(df_naca0012['cl']) / np.diff(np.deg2rad(df_naca0012['alpha']))
    gradient_0012 = np.mean(gradient_0012[:-4])
    print(f'raw gradient {gradient_0012}, difference {(gradient_0012 - theoretical_gradient) / theoretical_gradient}')
    gradient_2412 =  np.diff(df_naca2412['cl']) / np.diff(np.deg2rad(df_naca2412['alpha']))
    gradient_2412 = np.mean(gradient_2412[:-4])
    print(f'raw gradient {gradient_2412}, difference {(gradient_2412 - theoretical_gradient) / theoretical_gradient}')

    ax.grid( linestyle='--', linewidth=0.5)
    ax.legend()

    ax.set_xlabel('Angle of attack [deg]')
    ax.set_ylabel('Drag coefficient')

    fig.tight_layout()
    fig.savefig('report/final/figures/cd_alpha.png', dpi=300)


if __name__ == "__main__":

    plot_clcd_alpha()

    plt.show()

    