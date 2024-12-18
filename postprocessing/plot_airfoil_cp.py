
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

from scipy.signal import find_peaks
import mat73

def separate(arr):

    npts = arr.shape[0] // 2
    up_var = np.flip(arr[:npts+1], axis=0)
    low_var = arr[npts:]
    return up_var, low_var

def gradient_seperate(arr, gradarr):

    # take middle array between highest 2nd gradient points
    ddarr = np.abs(np.diff(gradarr, 3))
    # get indicies of two highest points
    idxs = np.argsort(ddarr)[-2:]
    idxs = np.sort(idxs)
    slic = arr[idxs[0]:idxs[1]]

    return slic

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

    cd_upper = -cpup[:-1] * np.sin(dtheta_u - alpha) * np.diff(xs_u)
    cd_lower = -cplo[:-1] * np.sin(dtheta_l - alpha) * np.diff(xs_l)

    Cl = np.sum(cl_lower) - np.sum(cl_upper)

    Cd = np.sum(cd_lower) + np.sum(cd_upper)

    return Cl, Cd


def plot_blades(av, gs):

    fig, ax = plt.subplots()

    # plot turbine blade
    data_dir = '../MEng-CW/IIB/4A3/turbine_cascade/data/'
    filename = '4A3_cascade_lwp26_jrt66_30-Oct-2024_1.mat'
    e = mat73.loadmat(data_dir + filename); e = e['e'];

    # Plot the suction surface
    xys = np.array(e['xy_ss'])
    xyp =  np.array(e['xy_ps'])
    # concate the two arrays but reverse the second one
    all_xy = np.concatenate((xys,xyp[::-1,:]),axis=0)
    all_xy = np.append(all_xy, [xys[0]], axis=0) # close the loop
    all_xy[:,1] *= -1 # flip y axis
    # offset all points by pos
    Cx = 1 # mm
    all_xy = Cx * all_xy
    
    # plot the blade as closed plot
    ax.plot(all_xy[:,0],all_xy[:,1], label = 'Experimental blade')

    # plot CFD blade
    assert av['casename'] == 'turbine_c'

    cut = cut_j(gs[0], 0)
    max_x = np.max(cut['x'])
    scaled_x = cut['x'] / max_x
    scaled_y =  cut['y'] / max_x
    _,yup = separate(scaled_y)
    plt.plot(scaled_x,scaled_y - yup[0], label = 'CFD blade')

    ax.legend()

def channels(N):
    # Convert all channel numbers to python indices

    # Loop over all keys, minus 1 and convert to integers
    for instr in N.keys():
        for chan in N[instr]:
            N[instr][chan] = N[instr][chan] - 1
            N[instr][chan] = N[instr][chan].astype(int)

    return(N)

def plot_experimental_cp(ax):
    # Use manometer data for detailed calculations on pressure distribution

    data_dir = '../MEng-CW/IIB/4A3/turbine_cascade/data/'
    filename = '4A3_cascade_lwp26_jrt66_30-Oct-2024_1.mat'

    e = mat73.loadmat(data_dir + filename); e = e['e'];
    b = mat73.loadmat(data_dir + 'cascade.mat'); b = b['b']['cfd'];
    N = channels(e['N'])

    # Initialise dictionary for processed data
    s = {}

    # Figure window for lift distribution

    cols = gen_cols()

    # Plot the CFD velocity distribution
    b['Vs_V2s'] = b['v'] * np.cos(np.deg2rad(b['alpha_2']))
    b['x_cx'] = (b['x'] - np.min(b['x'])) / (np.max(b['x']) - np.min(b['x'])) 

    # Concatenate the SS and PS measurements in a single complete loop
    iss = N['h2o']['Pss']; ips = N['h2o']['Pps'];
    leng = e['xy_ss'].shape[0]
    P = np.concatenate((e['P_h2o'][iss],np.flip(e['P_h2o'][ips]),
        np.expand_dims(e['P_h2o'][iss][0],axis=0)))

    # Join coordinates from both sides together
    s['xy'] = np.concatenate((e['xy_ss'],np.flip(e['xy_ps'],axis=0),
        np.expand_dims(e['xy_ss'][0],axis=0)))

    # Blade surface pressure coefficients
    
    Po_1 = e['P_h2o'][N['h2o']['Po_1']]; P_2 = 0;
    s['Cp'] = (Po_1 - P) / (Po_1 - P_2)

    ax.plot(s['xy'][:,0], s['Cp'], label='Experimental blade')

    return ax

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

    print(f'pstag_ref: {pstag_ref}')
    print(f'p_ref: {p_ref}')

    for i in range(len(gs)):
        gs[i]['cp'] = (gs[i]['p'] - p_ref) / (pstag_ref - p_ref)


    if av['casename'] == 'turbine_h':
        topcut = cut_j(gs[0], 0)
        botcut = cut_j(gs[0], -1)
        xs_u = gradient_seperate(topcut['x'], topcut['y'])
        ys_u = gradient_seperate(topcut['y'], topcut['y'])
        xs_l = gradient_seperate(botcut['x'], botcut['x'])
        ys_l = gradient_seperate(botcut['y'], botcut['x'])
        cpup = gradient_seperate(topcut['cp'], topcut['y'])
        cplo = gradient_seperate(botcut['cp'], botcut['x'])
    
    else:
        cut = cut_j(gs[0], 0)
        xs_u, xs_l = separate(cut['x'])
        ys_u, ys_l = separate(cut['y'])
        cpup, cplo = separate(cut['cp'])

    lens_u = np.cumsum(np.sqrt(np.diff(xs_u)**2 + np.diff(ys_u)**2))
    lens_l = np.cumsum(np.sqrt(np.diff(xs_l)**2 + np.diff(ys_l)**2))
    lens_u = np.insert(lens_u, 0, 0) / np.max(lens_u)
    lens_l = np.insert(lens_l, 0, 0) / np.max(lens_l)

    Cl, _ = calc_lift(av, gs)
    print(f'Cl: {Cl}')

    if av['casename'] == 'turbine_c':
        plot_blades(av, gs)

    fig = plt.figure(figsize=[8,5.4]); ax = plt.axes();

    ax.plot(lens_u, cpup, label='Upper surface')
    ax.plot(lens_l, cplo, label='Lower surface')

    # these are uncomparable
    #if av['casename'] == 'turbine_c':
    #    ax = plot_experimental_cp(ax)

    # flip y
    ax.invert_yaxis()
    ax.grid()
    ax.legend()

    ax.set_xlabel('Normalised surface path length [-]')
    ax.set_ylabel('Cp [-]')

    fig.tight_layout()

    fig.savefig(f'report/final/figures/{sys.argv[-1]}_surface_cp.png', dpi=300)

    # Show all the plots
    plt.show()

if __name__ == '__main__':
    main()


