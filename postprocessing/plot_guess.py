#
#   plot_guess                       
#                               
#   Script to plot an initial flowfield guess created using the 4A2 solver
#
#   Change to the directory you want to execute the script within and execute 
#   with "python path_to_script/plot_guess.py casename"

# Import modules and functions
from postprocessing.routines import *

def main():

    # Construct full filenames to read the guess data
    filename = 'cases/' + sys.argv[-1] + '/out_guess_' + sys.argv[-1] + '.bin'

    # Read the case from file
    gs = read_case(filename)

    # Open figure window and open four subplots
    fig,ax = plt.subplots(2,2,sharex=True,sharey=True,figsize=[14.4,7.2]); 
    fig.tight_layout()

    # Set subplot aspect ratios as equal and remove axes labels
    ax = ax.flatten()
    for a in ax:
        a.set_aspect('equal',adjustable='box'); a.axis('off')

    # Plot the primary flow variables to show the guess
    fieldnames = ['ro','roe','rovx','rovy']
    for n,name in enumerate(fieldnames):

        min_col = np.inf
        max_col = -np.inf

        for g in gs:
            min_col = min(min_col, np.min(g[name]))
            max_col = max(max_col, np.max(g[name]))
 
        # Plot filled contour levels
        for g in gs:
            hc = ax[n].pcolormesh(g['x'],g['y'],g[name],shading='gouraud', vmax = max_col, vmin = min_col)

  	# Add colorbar with variable name
        colorbar(hc,name)

        # Draw the walls of the block
        for g in gs:
            plot_wall(ax[n],g)

    # Show all the plots
    plt.show()

    
main()


