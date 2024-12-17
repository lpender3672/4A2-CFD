
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

def plot_cl_alpha():

    av_default = read_settings('cases/naca0012/input_naca0012.txt')

    alphas = np.linspace(-5, 15, 3)

    avs = []
    for alpha in alphas:
        av_default['alpha'] = alpha
        avs.append(av_default.copy())

    manager = create_cfd_env(avs, 'cl_alpha.csv')
    manager.clear_worker_folders()
    manager.start_workers()

    lift_file = manager.shared_file.parent / "cl_alpha.csv"
    df = pd.read_csv(lift_file)

    print(df)



if __name__ == "__main__":

    plot_cl_alpha()
    