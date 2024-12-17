
#
#   plot_contours
#                               
#   Script to plot a converged flowfield from the 4A2 solver
#
#   Change to the directory you want to execute the script within and execute 
#   with "python path_to_script/plot_contours.py casename"
import os
# Import modules and functions
from postprocessing.routines import *


def separate(arr):

    npts = arr.shape[0] // 2
    up_var = np.flip(arr[:npts+1], axis=0)
    low_var = arr[npts:]
    return up_var, low_var

def calc_lift(av, gs):

    for i in range(len(gs)):
        gs[i] = calc_secondary(av,gs[i])

    cut = cut_i(gs[0], 0)
    pstag_ref = mass_av(cut, 'pstag')[0]
    p_ref = area_av(cut, 'p')[0]

    for i in range(len(gs)):
        gs[i]['cp'] = (gs[i]['p'] - p_ref) / (pstag_ref - p_ref)

    cut = cut_j(gs[0], 0) 

    xs_u, xs_l = separate(cut['x'])
    ys_u, ys_l = separate(cut['y'])

    cpup, cplo = separate(cut['cp'])

    alpha = av['alpha'] * np.pi / 180
    dtheta_u = np.arctan2(np.diff(ys_u), np.diff(xs_u))
    dtheta_l = np.arctan2(np.diff(ys_l), np.diff(xs_l))

    cl_upper = -cpup[:-1] * np.cos(dtheta_u - alpha) * np.diff(xs_u)
    cl_lower = -cplo[:-1] * np.cos(dtheta_l - alpha) * np.diff(xs_l)

    Cl = np.sum(cl_lower) - np.sum(cl_upper)

    return Cl

def main():

    # Construct full filenames to read the run data
    inname = 'cases/' + sys.argv[-1] + '/input_' + sys.argv[-1] + '.txt'
    outname = 'cases/' + sys.argv[-1] + '/out_final_' + sys.argv[-1] + '.bin'

    # Read the settings and the case from file
    av = read_settings(inname)
    gs = read_case(outname)

    # When presenting results all values should be non-dimensionalised. Two
    # variables of interest might be:
    #    1. Static pressure coefficient, (p - p_ref) / (pstag_ref - p_ref)
    #    2. Mach number, v / (ga * rgas * t)**0.5

    # First complete the "calc_secondary" function within "routines.py" to
    # calculate static pressure and Mach number, and any others you want!

    for i in range(len(gs)):
        gs[i] = calc_secondary(av,gs[i])    

    # Use the "cut_i", "mass_av" AND "area_av" functions to calculate the
    # reference pressures at the inlet plane and therefore the static pressure
    # coefficient
    cut = cut_i(gs[0], 0)
    pstag_ref = mass_av(cut, 'pstag')[0]
    p_ref = area_av(cut, 'p')[0]

    for i in range(len(gs)):
        gs[i]['cp'] = (gs[i]['p'] - p_ref) / (pstag_ref - p_ref)

    cut = cut_j(gs[0], 0)

    fig = plt.figure(figsize=[9.6,3.8]); ax = plt.axes();

    xs_u, xs_l = separate(cut['x'])
    ys_u, ys_l = separate(cut['y'])

    lens_u = np.cumsum(np.sqrt(np.diff(xs_u)**2 + np.diff(ys_u)**2))
    lens_l = np.cumsum(np.sqrt(np.diff(xs_l)**2 + np.diff(ys_l)**2))
    lens_u = np.insert(lens_u, 0, 0)
    lens_l = np.insert(lens_l, 0, 0)

    cpup, cplo = separate(cut['cp'])

    Cl = calc_lift(av, gs)
    print(f'Cl: {Cl}')

    ax.plot(lens_u, cpup)
    ax.plot(lens_l, cplo)

    

    # flip y
    ax.invert_yaxis()
    ax.grid()

    ax.set_xlabel('Upper surface [m]')
    ax.set_ylabel('Cp [-]')

    fig.tight_layout()

#    fig.savefig(f'report/final/figures/{sys.argv[-1]}_{name}.png', dpi=300)

    # Show all the plots
    plt.show()

if __name__ == '__main__':
    main()


