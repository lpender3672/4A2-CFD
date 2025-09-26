
#include "routines.h"
#include <cstring>

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

    av.crashed = 0;

    infile >> av.rgas >> av.gam;
    infile >> av.cfl >> av.sfac >> av.sfac_res >> av.d_max >> av.d_var >> av.facsec >> av.fcorr;
    infile >> av.nsteps >> av.nkruts >> av.guess_method >> av.tstep_method;
    infile >> av.ni >> av.nj;

    av.cp = av.rgas * av.gam / (av.gam - 1.0);
    av.cv = av.cp / av.gam;
    av.fgam = (av.gam - 1.0) / av.gam;

    infile >> bcs.pstag >> bcs.tstag >> bcs.alpha >> bcs.rfin;

    infile >> bcs.p_out;

    infile.close();
}

void write_settings(const std::string& fpath, const t_appvars& av, const t_bconds& bcs) {
    std::ofstream outfile(fpath);

    outfile << av.casename << std::endl;
    outfile << av.rgas << " " << av.gam << std::endl;
    outfile << av.cfl << " " << av.sfac << " " << av.sfac_res << " " << av.d_max << " " << av.d_var << " " << av.facsec << " " << av.fcorr << std::endl;
    outfile << av.nsteps << " " << av.nkruts << " " << av.guess_method << " " << av.tstep_method << std::endl;
    outfile << av.ni << " " << av.nj << std::endl;
    outfile << bcs.pstag << " " << bcs.tstag << " " << bcs.alpha << " " << bcs.rfin << std::endl;
    outfile << bcs.p_out << std::endl;

    outfile.close();

}
