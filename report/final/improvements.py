
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


def test_rk4():
    # test improvements on bump case
    av = read_settings('cases/bump/input_bump.txt')

    av['facsec'] = 0
    av['fcorr'] = 0
    av['nrkuts'] = 1
    av['guess_method'] = 2

    # create 8 workers
    avs = []
    nkruts = [1, 2, 3, 4, 5]
    cfls = np.logspace(-2, np.log10(0.5), 10, endpoint = True)
    cfls = np.append(cfls, 0.8)
    cfls = np.append(cfls, 1.0)
    cfls = np.append(cfls, 1.1)
    cfls = np.append(cfls, 1.3)

    for cfl in cfls:
        av['cfl'] = cfl
        for nrkut in nkruts:
            av['nrkuts'] = nrkut
            avs.append(av.copy())

    manager = create_cfd_env(avs, "rkN_bump.csv")
    manager.clear_worker_folders()
    manager.start_workers()

    return manager



# now load and plot residuals

def plot_rkn_accuracy():

    manager = test_rk4()
    df = pd.read_csv(manager.shared_file)

    nkruts = [1, 2, 3, 4, 5]

    for nrkut in nkruts:
        grouped_df = df.groupby('nrkuts').get_group(nrkut)
        plt.loglog(
            grouped_df['cfl'],
            grouped_df['dro_avg'],
            label = 'Average residual'
        )

    print(df['dro_avg'])
    plt.show()

def plot_rkn_time():

    manager = test_rk4()
    df = pd.read_csv(manager.shared_file)
    df = df[df['converged'] < 2]

    group_headers = list(manager.avs[0].keys())

    fig, ax = plt.subplots()

    nkruts = [1, 2, 3, 4, 5]
    averaged_df = df.groupby(group_headers, as_index=False).apply(
        lambda g : pd.Series({
            'dt' : g['dt'].mean()
        })
    )

    print(len(averaged_df), len(df))

    for nrkut in nkruts:
        grouped_df = averaged_df.groupby('nrkuts').get_group(nrkut)
        ax.loglog(
            grouped_df['cfl'],
            grouped_df['dt'],
            label = f'nrkuts = {nrkut}'
        )

    ax.grid( which='both', linestyle='--', linewidth=0.5)
    ax.set_xlabel('CFL')
    ax.set_ylabel('Run time (s)')
    ax.legend()

    fig.tight_layout()
    fig.savefig('report/final/figures/rk4_time.png', dpi=300)
    plt.show()


if __name__ == "__main__":

    plot_rkn_time()