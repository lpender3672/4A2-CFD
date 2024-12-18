
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


def get_improvement_setting_templates():

    av_default = read_settings('cases/bump/input_bump.txt')
    av_default['cfl'] = 0.2
    av_default['sfac'] = 0.5

    av_default['sfac_res'] = 0
    av_default['facsec'] = 0
    av_default['fcorr'] = 0
    av_default['nrkuts'] = 1
    av_default['guess_method'] = 2
    av_default['tstep_method'] = 1

    av_rk4 = av_default.copy()

    av_rk4['sfac_res'] = 0
    av_rk4['facsec'] = 0
    av_rk4['fcorr'] = 0
    av_rk4['nrkuts'] = 4
    av_rk4['guess_method'] = 2
    av_rk4['tstep_method'] = 1

    av_rsfac = av_default.copy()
    
    av_rsfac['sfac_res'] = 0.5
    av_rsfac['facsec'] = 0
    av_rsfac['fcorr'] = 0
    av_rsfac['nrkuts'] = 1
    av_rsfac['guess_method'] = 2
    av_rsfac['tstep_method'] = 1

    av_fcorr = av_default.copy()

    av_fcorr['sfac_res'] = 0
    av_fcorr['facsec'] = 0
    av_fcorr['fcorr'] = 0.8
    av_fcorr['nrkuts'] = 1
    av_fcorr['guess_method'] = 2
    av_fcorr['tstep_method'] = 1

    av_tinac = av_default.copy()

    av_tinac['sfac_res'] = 0
    av_tinac['facsec'] = 0
    av_tinac['fcorr'] = 0
    av_tinac['nrkuts'] = 1
    av_tinac['guess_method'] = 2
    av_tinac['tstep_method'] = 2

    templates = [av_default, av_rk4, av_rsfac, av_fcorr, av_tinac]
    return templates

def plot_improvement(manager, name):
    df = pd.read_csv(manager.shared_file)

    fig, ax = plt.subplots()

    df = df[df['converged'] < 2]

    default_headers = ['facsec', 'fcorr', 'nrkuts', 'guess_method', 'tstep_method']
    default_items = [0, 0, 1, 2, 1]
    default_df = df.groupby(default_headers, as_index=False).get_group(tuple(default_items))
    default_df = default_df.sort_values(name)
    
    df_rk4 = df[df['nrkuts'] == 4].sort_values(name)
    df_rsfac = df[df['sfac_res'] == 0.5].sort_values(name)
    df_fcorr = df[df['fcorr'] == 0.8].sort_values(name)
    df_tinac = df[df['tstep_method'] == 2].sort_values(name)

    default_df_r = default_df.groupby(name, as_index=False).agg({'dro_avg': 'mean'})
    df_rk4_r = df_rk4.groupby(name, as_index=False).agg({'dro_avg': 'mean'})
    df_rsfac_r = df_rsfac.groupby(name, as_index=False).agg({'dro_avg': 'mean'})
    df_fcorr_r = df_fcorr.groupby(name, as_index=False).agg({'dro_avg': 'mean'})
    df_tinac_r = df_tinac.groupby(name, as_index=False).agg({'dro_avg': 'mean'})

    ax.loglog(
        df_rk4_r[name].to_numpy(),
        df_rk4_r['dro_avg'].to_numpy(),
        label = 'RK4'
    )
    ax.loglog(
        df_rsfac_r[name].to_numpy(),
        df_rsfac_r['dro_avg'].to_numpy(),
        label = 'sfac_res'
    )
    ax.loglog(
        df_fcorr_r[name].to_numpy(),
        df_fcorr_r['dro_avg'].to_numpy(),
        label = '$f_{corr}$'
    )
    ax.loglog(
        df_tinac_r[name].to_numpy(),
        df_tinac_r['dro_avg'].to_numpy(),
        label = '$\Delta t_{i,j}$'
    )
    ax.loglog(
        default_df_r[name].to_numpy(),
        default_df_r['dro_avg'].to_numpy(),
        label = 'Default'
    )

    # perform averaging
    default_df_t = default_df.groupby(name, as_index=False).agg({'dt': 'mean'})
    df_rk4_t = df_rk4.groupby(name, as_index=False).agg({'dt': 'mean'})
    df_rsfac_t = df_rsfac.groupby(name, as_index=False).agg({'dt': 'mean'})
    df_fcorr_t = df_fcorr.groupby(name, as_index=False).agg({'dt': 'mean'})
    df_tinac_t = df_tinac.groupby(name, as_index=False).agg({'dt': 'mean'})

    ax.grid( which='both', linestyle='--', linewidth=0.5)
    ax.set_xlabel(name)
    ax.set_ylabel('Average residual')

    ax.legend()
    fig.tight_layout()

    fig.savefig(f'report/final/figures/improvements_{name}_residual.png', dpi=300)

    fig, ax = plt.subplots()

    # perform averaging

    ax.loglog(
        df_rk4_t[name].to_numpy(),
        df_rk4_t['dt'].to_numpy(),
        label = 'RK4'
    )
    ax.loglog(
        df_rsfac_t[name].to_numpy(),
        df_rsfac_t['dt'].to_numpy(),
        label = 'sfac_res'
    )
    ax.loglog(
        df_fcorr_t[name].to_numpy(),
        df_fcorr_t['dt'].to_numpy(),
        label = '$f_{corr}$'
    )
    ax.loglog(
        df_tinac_t[name].to_numpy(),
        df_tinac_t['dt'].to_numpy(),
        label = '$\Delta t_{i,j}$'
    )
    ax.loglog(
        default_df_t[name].to_numpy(),
        default_df_t['dt'].to_numpy(),
        label = 'Default'
    )

    ax.grid( which='both', linestyle='--', linewidth=0.5)
    ax.set_xlabel(name)
    ax.set_ylabel('Run time (s)')
    ax.legend()
    fig.tight_layout()

    fig.savefig(f'report/final/figures/improvements_{name}_time.png', dpi=300)

def plot_improvement_cfl():

    templates = get_improvement_setting_templates()

    cfls = np.logspace(-2, np.log10(0.5), 10, endpoint = True)
    avs = []
    for avt in templates:
        for cfl in cfls:
            avt['cfl'] = cfl
            avs.append(avt.copy())

    manager = create_cfd_env(avs, 'improve_comparison_cfl.csv')
    manager.clear_worker_folders()
    manager.start_workers()

    plot_improvement(manager, 'cfl')

def plot_improvement_ni():

    templates = get_improvement_setting_templates()

    nis = np.logspace(1, 3, 10, endpoint = True).astype(int)
    nis = np.insert(nis, 0, [5, 8])

    avs = []
    for avt in templates:
        for ni in nis:
            avt['ni'] = ni
            avs.append(avt.copy())

    manager = create_cfd_env(avs, 'improve_comparison_ni.csv')
    manager.clear_worker_folders()
    manager.start_workers()

    plot_improvement(manager, 'ni')

def plot_smoothing_cfl(av_template, data):
    
    # keep sfac = 0.8

    sfacs = np.linspace(0.05, 0.8, 10, endpoint = True)
    #cfls = np.logspace(-2, np.log10(1.5), 10, endpoint = True)
    # rewrite with arrange
    cfls = 10**np.arange(-2, np.log10(1.5), 0.2)
    print(cfls)

    avs = []
    for cfl in cfls:
        for sfac in sfacs:
            av_template['cfl'] = cfl
            av_template['sfac'] = sfac
            avs.append(av_template.copy())

    manager = create_cfd_env(avs, 'smoothing_cfl.csv')
    manager.clear_worker_folders()
    manager.start_workers()

    plot_scatter(manager, 'cfl', 'sfac', data)


def plot_smoothing_fcorr(av_template, data):
    
    # keep sfac = 0.8
    av_template['cfl'] = 0.3

    sfacs = np.linspace(0.05, 0.8, 10, endpoint = True)
    #cfls = np.logspace(-2, np.log10(1.5), 10, endpoint = True)
    # rewrite with arrange
    fcorrs = np.linspace(0.1, 0.9, 10, endpoint = True)

    avs = []
    for fcorr in fcorrs:
        for sfac in sfacs:
            av_template['fcorr'] = fcorr
            av_template['sfac'] = sfac
            avs.append(av_template.copy())

    manager = create_cfd_env(avs, 'smoothing_fcorr.csv')
    manager.clear_worker_folders()
    manager.start_workers()

    plot_scatter(manager, 'sfac', 'fcorr', data, logx=False)


def plot_smoothing_sfac_res(av_template, data):
    
    # keep sfac = 0.8
    av_template['cfl'] = 0.3

    sfacs = np.linspace(0.05, 0.95, 10, endpoint = True)
    #cfls = np.logspace(-2, np.log10(1.5), 10, endpoint = True)
    # rewrite with arrange
    sfac_ress = np.linspace(0.3, 0.95, 10, endpoint = True)

    avs = []
    for sfac_res in sfac_ress:
        for sfac in sfacs:
            av_template['sfac_res'] = sfac_res
            av_template['sfac'] = sfac
            avs.append(av_template.copy())

    manager = create_cfd_env(avs, 'smoothing_sfac_res.csv')
    manager.clear_worker_folders()
    manager.start_workers()

    plot_scatter(manager, 'sfac', 'sfac_res', data, logx=False)

def plot_scatter(manager, xlabel, ylabel, clabel, logx=True):

    df = pd.read_csv(manager.shared_file)

    df_converged_within = df[df['converged'] == 0]
    df_converged_outside = df[df['converged'] == 1]
    df_diverged = df[df['converged'] == 2]
    df_max_iter = df[df['converged'] == 3]

    # scatter
    fig, ax = plt.subplots(figsize=[12,9])

    fig.subplots_adjust(
        top=0.872,
        bottom=0.087,
        left=0.073,
        right=0.983,
        hspace=0.2,
        wspace=0.2
    )

    if clabel == 'FM':
        cval_in = - np.log10(df_converged_within['dro_avg']) / df_converged_within['dt']
        cval_out = - np.log10(df_converged_outside['dro_avg']) / df_converged_outside['dt']
    else:
        cval_in = np.log10(df_converged_within[clabel])
        cval_out = np.log10(df_converged_outside[clabel])

    conin = ax.scatter(
        df_converged_within[xlabel].to_numpy(),
        df_converged_within[ylabel].to_numpy(),
        c =  cval_in,
        s = 100,
        label = 'Converged within',
        marker = 'o',
        cmap = 'jet'
    )
    conout = ax.scatter(
        df_converged_outside[xlabel].to_numpy(),
        df_converged_outside[ylabel].to_numpy(),
        c = cval_out,
        s = 100,
        label = 'Converged outside',
        marker = 'd',
        cmap = 'jet'
    )
    ax.scatter(
        df_diverged[xlabel].to_numpy(),
        df_diverged[ylabel].to_numpy(),
        s = 100,
        color = 'red',
        label = 'Diverged',
        marker = 'x'
    )
    ax.scatter(
        df_max_iter[xlabel].to_numpy(),
        df_max_iter[ylabel].to_numpy(),
        s = 100,
        color = 'black',
        label = 'Max iterations',
        marker = 's'
    )

    if len(df_converged_within) > len(df_converged_outside):
        cbar = plt.colorbar(conin)
    else:
        cbar = plt.colorbar(conout)

    if clabel == 'dro_avg':
        cbar.set_label('Log10 average residual density error')
    elif clabel == 'dt':
        cbar.set_label('Log10 run time')
    elif clabel == 'FM':
        cbar.set_label('FM')

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if logx:
        ax.set_xscale('log')

    ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15),
          ncol=2, fancybox=True)
    fig.tight_layout()

    fig.savefig(f'report/final/figures/{xlabel}_{ylabel}_{clabel}.png', dpi=300)


def effort_vs_accuracy_fcorr():

    df = pd.read_csv('report/final/data/smoothing_fcorr.csv')

    df = df[df['converged'] < 2]

    fig, ax = plt.subplots(figsize = [6, 4])

    scat = ax.scatter(
        np.log10(df['dro_avg']),
        np.log10(df['dt']),
        c = df['sfac'],
        s = 60*df['fcorr']
    )

    cbar = plt.colorbar(scat)

    cbar.set_label('sfac')
    ax.set_xlabel('Log10 average residual density error')
    ax.set_ylabel('Log10 run time')

    ax.grid( which='both', linestyle='--', linewidth=0.5)

    fig.tight_layout()

if __name__ == "__main__":

    #plot_improvement_cfl()
    #plot_improvement_ni()

    av_template = read_settings('cases/bump/input_bump.txt')
    av_template['fcorr'] = 0.8
    print(av_template)
    #plot_smoothing_cfl(av_template, 'dro_avg')
    #plot_smoothing_cfl(av_template, 'dt')
    #plot_smoothing_fcorr(av_template, 'dt')
    #plot_smoothing_sfac_res(av_template, 'dt')
    #plot_smoothing_cfl_residual(av_template, 'dt')

    effort_vs_accuracy_fcorr()

    plt.show()

