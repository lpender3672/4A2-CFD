
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

    av_templates = [av_0012, av_2412]
    alphas = np.linspace(-5, 15, 9)

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

    fig, ax = plt.subplots()

    ax.plot(df_naca0012['alpha'], df_naca0012['cl'], label='NACA0012')
    ax.plot(df_naca2412['alpha'], df_naca2412['cl'], label='NACA2412')

    ax.grid( linestyle='--', linewidth=0.5)
    ax.legend()

    ax.set_xlabel('Angle of attack [deg]')
    ax.set_ylabel('Lift coefficient')

    fig.tight_layout()
    fig.savefig('report/final/figures/cl_alpha.png', dpi=300)


if __name__ == "__main__":

    plot_clcd_alpha()

    plt.show()

    