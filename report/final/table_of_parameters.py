

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




# basically read all cases and format latex table

def create_table():

    cases = [
        "bend",
        "bump",
        "naca0012",
        "naca2412",
        "turbine_c",
        "turbine_h"
    ]

    #create empty df
    df = pd.DataFrame(columns=default_settings("").keys())

    for casename in cases:
        av = read_settings(f'cases/{casename}/input_{casename}.txt')
        df.loc[len(df)] = av

    print(df.columns)
    # drop some columns
    df = df.drop(columns=['rgas', 'gam', 'facsec', 'nrkuts', 'tstep_method', 'guess_method', 'd_max', 'd_var'])
    # set dtypes correctly
    df = df.astype({
        'alpha': 'int',
        'cfl': 'float',
        'sfac_res' : 'float',
        'fcorr' : 'float',
        'nsteps' : 'int'
    })
    # round sig figs not decimals
    
    latex_output = (
    df.style
    .format({
        'alpha': '{:.1f}', 
        'cfl': '{:.2f}', 
        'sfac' : '{:.2f}',
        'sfac_res': '{:.2f}', 
        'fcorr': '{:.2f}', 
        'nsteps': '{:.0e}', 
        'pstag' : '{:.2e}',
        'tstag' : '{:.2f}',
        'p' : '{:.2e}',
        'rfin' : '{:.2f}',

    })
    .hide(axis='index')  # Removes the index
    .to_latex()
    )

    latex_output = latex_output.replace('-1', '-')
    
    for header in df.columns:
        if header == 'p' or header == 'sfac':
            continue
        texttt = r'\texttt{' + header + '}'
        latex_output = latex_output.replace(
            header, texttt)
        
    latex_output = latex_output.replace(
        '_', r'\_'
    )

    print(latex_output)

if __name__ == "__main__":
    create_table()