
import os
from pathlib import Path
from postprocessing.routines import (
    read_settings,
    default_settings,
)
from runtime_enviroment import CFDManager

import numpy as np
import pandas as pd


def create_cfd_env():

    n_workers = 2
    base_folder = Path(os.getcwd()) / "case_env"
    data_file = Path(os.getcwd()) / "report" / "final" / "data" / "results.csv"
    cfd_executable = 'build/solverApp.exe'

    #av = read_settings('cases/bump/input_bump.txt')
    av = default_settings('bump')
    av['nsteps'] = 500
    avs = [av for _ in range(4)]

    manager = CFDManager(n_workers, base_folder, data_file, cfd_executable, avs)

    return manager
    

manager = create_cfd_env()
manager.clear_worker_folders()
manager.start_workers()