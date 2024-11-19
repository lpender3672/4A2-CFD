
#include "routines.h"

#define M_PI 3.14159265358979323846

void read_settings(const std::string& fpath, t_appvars& av, t_bconds& bcs) {
    std::ifstream infile(fpath);
    if (!infile.is_open()) {
        std::cerr << "Error opening file: " << fpath << std::endl;
        return;
    }

    // Read and trim the case name
    std::string tempname;
    std::getline(infile, tempname);
    std::strncpy(av.casename, tempname.c_str(), 128);
    av.casename[128 - 1] = '\0';

    size_t last_slash_index = fpath.find_last_of('/');
    if (last_slash_index != std::string::npos) {
        std::string folderpath = fpath.substr(0, last_slash_index);
        std::strncpy(av.casefolder, folderpath.c_str(), 128);
        av.casefolder[128 - 1] = '\0';
    } else {
        av.casefolder[0] = '\0'; 
    }

    av.crashed = false;

    infile >> av.rgas >> av.gam;
    infile >> av.cfl >> av.sfac >> av.d_max >> av.d_var >> av.facsec >> av.fcorr;
    infile >> av.nsteps;
    infile >> av.ni >> av.nj;

    av.cp = av.rgas * av.gam / (av.gam - 1.0);
    av.cv = av.cp / av.gam;
    av.fgam = (av.gam - 1.0) / av.gam;

    av.sfac *= av.cfl;
    av.d_max *= av.cfl;
    av.d_avg = 0.5 * av.d_max;

    infile >> bcs.pstag >> bcs.tstag >> bcs.alpha >> bcs.rfin;

    bcs.alpha = bcs.alpha * M_PI / 180.0;

    bcs.rostag = bcs.pstag / (av.rgas * bcs.tstag);

    infile >> bcs.p_out;

    infile.close();
}
