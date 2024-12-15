
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

    fig = plt.figure(figsize=[9.6,3.8]); ax = plt.axes();

    cut = cut_j(gs[0], 0)
    slen = np.sqrt(cut['lx']**2 + cut['ly']**2)
    n = slen.shape[0] // 2
    slo = np.cumsum(slen[:n][::-1])
    sup = np.cumsum(slen[n:])
    ax.plot(slo, cut['cp'][:n][::-1])
    ax.plot(sup, cut['cp'][n+1:])
    #ax.plot(cut['x'][:10], cut['y'][:10])
    

    fig.tight_layout()

#    fig.savefig(f'report/final/figures/{sys.argv[-1]}_{name}.png', dpi=300)

    # Show all the plots
    plt.show()

if __name__ == '__main__':
    main()


