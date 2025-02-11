
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

    n_workers = 2
    base_folder = Path(os.getcwd()) / "case_env"
    data_file = Path(os.getcwd()) / "report" / "final" / "data" / name
    cfd_executable = 'build/solverApp.exe'

    #av = read_settings('cases/bump/input_bump.txt')

    manager = CFDManager(n_workers, base_folder, data_file, cfd_executable, avs)

    return manager


def plot_effort_accuracy():

    av_template = read_settings('cases/bump/input_bump.txt')

    av_template['facsec'] = 0
    av_template['fcorr'] = 0.8
    av_template['nrkuts'] = 4
    av_template['guess_method'] = 2
    av_template['tstep_method'] = 2


