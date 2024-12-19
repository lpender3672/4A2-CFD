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
    if "turbine" in av['casename']:
        cut = cut_i(gs[0], -1)
    else:
        cut = cut_i(gs[0], 0)
        
    pstag_ref = mass_av(cut, 'pstag')[0]
    p_ref = area_av(cut, 'p')[0]

    for i in range(len(gs)):
        gs[i]['cp'] = (gs[i]['p'] - p_ref) / (pstag_ref - p_ref)

    # Specify the parameters to plot
    fieldnames = ['cp', 'mach']; 
    colnames = ['Static pressure coefficient','Mach number']

    sizedict = {
        "bump" : [9.6,3.8],
        "bend" : [9.6,3.8],
        "tube" : [9.6,3.8],
        "tunnel" : [9.6,3.8],
        "waves" : [9.6,3.8],
        "turbine_h" : [6,6],
        "turbine_c" : [6,6],
    }
    if "naca" in sys.argv[-1]:
        figsize = [8,5.4]
    else:
        figsize = sizedict[sys.argv[-1]]

    # Plot the calculated non-dimensional parameters to show the flow solution
    for n,name in enumerate(fieldnames):

        # Open figure window
        fig = plt.figure(figsize=figsize); ax = plt.axes();
    
        # Set aspect ratio as equal and remove axes labels
        ax.set_aspect('equal',adjustable='box'); ax.axis('off')

        min_col = np.inf
        max_col = -np.inf

        for g in gs:
            min_col = min(min_col, np.min(g[name]))
            max_col = max(max_col, np.max(g[name]))
 
        # Plot filled contour levels
        for g in gs:
            hc = ax.pcolormesh(g['x'],g['y'],g[name],shading='gouraud', vmax = max_col, vmin = min_col)        

        # Add colorbar with variable name
        colorbar(hc,colnames[n])

        # Add Mach = 1 contours
        if name == 'mach':
            for g in gs:
                ax.contour(g['x'],g['y'],g['mach'],[1.0],colors='w',
                    linewidths=0.5)

        # Draw the walls of the block
        for g in gs:
            plot_wall(ax,g)

        if "naca" in sys.argv[-1]:
            # zoom in on the airfoil
            ax.set_xlim([-0.5, 1.5])
            ax.set_ylim([-0.6, 0.6])

        fig.tight_layout()

        fig.savefig(f'report/final/figures/{sys.argv[-1]}_{name}.png', dpi=300)

    # Show all the plots
    plt.show()

if __name__ == '__main__':
    main()


