
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
    av['tstep_method'] = 1

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


def get_improvement_setting_templates():

    av_default = read_settings('cases/bump/input_bump.txt')

    av_default['facsec'] = 0
    av_default['fcorr'] = 0
    av_default['nrkuts'] = 1
    av_default['guess_method'] = 2
    av_default['tstep_method'] = 1

    av_rk4 = read_settings('cases/bump/input_bump.txt')

    av_rk4['facsec'] = 0
    av_rk4['fcorr'] = 0
    av_rk4['nrkuts'] = 4
    av_rk4['guess_method'] = 2
    av_rk4['tstep_method'] = 1

    av_facsec = read_settings('cases/bump/input_bump.txt')
    
    av_facsec['facsec'] = 0.5
    av_facsec['fcorr'] = 0
    av_facsec['nrkuts'] = 1
    av_facsec['guess_method'] = 2
    av_facsec['tstep_method'] = 1

    av_fcorr = read_settings('cases/bump/input_bump.txt')

    av_fcorr['facsec'] = 0
    av_fcorr['fcorr'] = 0.8
    av_fcorr['nrkuts'] = 1
    av_fcorr['guess_method'] = 2
    av_fcorr['tstep_method'] = 1

    av_tinac = read_settings('cases/bump/input_bump.txt')

    av_tinac['facsec'] = 0
    av_tinac['fcorr'] = 0
    av_tinac['nrkuts'] = 1
    av_tinac['guess_method'] = 2
    av_tinac['tstep_method'] = 2

    templates = [av_default, av_rk4, av_facsec, av_fcorr, av_tinac]
    return templates

def plot_improvement_cfl():

    templates = get_improvement_setting_templates()

    cfls = np.logspace(-2, np.log10(0.5), 10, endpoint = True)
    avs = []
    for avt in templates:
        for cfl in cfls:
            avt['cfl'] = cfl
            avs.append(avt.copy())

    manager = create_cfd_env(avs, 'improve_comparison.csv')
    manager.clear_worker_folders()
    manager.start_workers()

    df = pd.read_csv(manager.shared_file)

    fig, ax = plt.subplots()

    df = df[df['converged'] < 2]
    
    df_rk4 = df[df['nrkuts'] == 4].sort_values('cfl')
    df_facsec = df[df['facsec'] == 0.5].sort_values('cfl')
    df_fcorr = df[df['fcorr'] == 0.8].sort_values('cfl')
    df_tinac = df[df['tstep_method'] == 2].sort_values('cfl')

    ax.loglog(
        df_rk4['cfl'].to_numpy(),
        df_rk4['dro_avg'].to_numpy(),
        label = 'RK4'
    )
    """ax.loglog(
        df_facsec['cfl'].to_numpy(),
        df_facsec['dro_avg'].to_numpy(),
        label = 'Facsec'
    )"""
    ax.loglog(
        df_fcorr['cfl'].to_numpy(),
        df_fcorr['dro_avg'].to_numpy(),
        label = '$f_{corr}$'
    )
    ax.loglog(
        df_tinac['cfl'].to_numpy(),
        df_tinac['dro_avg'].to_numpy(),
        label = 'Variable $\Delta t$'
    )

    default_headers = ['facsec', 'fcorr', 'nrkuts', 'guess_method', 'tstep_method']
    default_items = [0, 0, 1, 2, 1]
    default_df = df.groupby(default_headers, as_index=False).get_group(tuple(default_items))
    default_df = default_df.sort_values('cfl')
    ax.loglog(
        default_df['cfl'].to_numpy(),
        default_df['dro_avg'].to_numpy(),
        label = 'Default'
    )

    # perform averaging
    default_df = default_df.groupby('cfl', as_index=False).agg({'dt': 'mean'})
    df_rk4 = df_rk4.groupby('cfl', as_index=False).agg({'dt': 'mean'})
    df_facsec = df_facsec.groupby('cfl', as_index=False).agg({'dt': 'mean'})
    df_fcorr = df_fcorr.groupby('cfl', as_index=False).agg({'dt': 'mean'})
    df_tinac = df_tinac.groupby('cfl', as_index=False).agg({'dt': 'mean'})

    ax.grid( which='both', linestyle='--', linewidth=0.5)
    ax.set_xlabel('CFL')
    ax.set_ylabel('Average residual')

    ax.legend()
    fig.tight_layout()

    fig.savefig('report/final/figures/improvements_cfl_residual.png', dpi=300)

    fig, ax = plt.subplots()

    # perform averaging

    ax.loglog(
        df_rk4['cfl'].to_numpy(),
        df_rk4['dt'].to_numpy(),
        label = 'RK4'
    )
    """ax.loglog(
        df_facsec['cfl'].to_numpy(),
        df_facsec['dt'].to_numpy(),
        label = 'Facsec'
    )"""
    ax.loglog(
        df_fcorr['cfl'].to_numpy(),
        df_fcorr['dt'].to_numpy(),
        label = '$f_{corr}$'
    )
    ax.loglog(
        df_tinac['cfl'].to_numpy(),
        df_tinac['dt'].to_numpy(),
        label = '$\Delta t_{i,j}$'
    )
    ax.loglog(
        default_df['cfl'].to_numpy(),
        default_df['dt'].to_numpy(),
        label = 'Default'
    )

    ax.grid( which='both', linestyle='--', linewidth=0.5)
    ax.set_xlabel('CFL')
    ax.set_ylabel('Run time (s)')
    ax.legend()
    fig.tight_layout()

    fig.savefig('report/final/figures/improvements_cfl_time.png', dpi=300)

    plt.show()

if __name__ == "__main__":

    plot_improvement_cfl()