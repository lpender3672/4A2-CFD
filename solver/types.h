#include <string>

#pragma once

#ifndef TYPES_H
#define TYPES_H

extern "C" {
    struct t_grid {
        int ni, nj;

        // Pointers for 2D array data
        float* x;
        float* y;
        float* area;
        float* lx_i;
        float* ly_i;
        float* lx_j;
        float* ly_j;
        float l_min;

        // Primary variables at nodes
        float* ro;
        float* roe;
        float* rovx;
        float* rovy;

        // Variables to hold cell increments
        float* dro;
        float* droe;
        float* drovx;
        float* drovy;

        // Secondary variables at nodes
        float* p;
        float* hstag;
        float* vx;
        float* vy;

        int* wall;
    };

    struct t_appvars {
        char casename[128];
        char casefolder[128];
        
        float rgas, gam, cp, cv, fgam;
        float cfl, sfac, dt, d_max, d_avg;
        float d_var, facsec, fcorr;
        int nsteps, nstep, nkruts, guess_method;
        float ro_ref, roe_ref, rov_ref;
        int ni, nj;
        int nn, nm;
        bool crashed;

    };
    

    struct t_bconds {
        float pstag, tstag, alpha, rfin, rostag;

        float* p;
        float* ro;

        float p_out;
        int n_in, n_out;
    };

    struct t_conv_point {
        int iter;
        float d_max, d_avg;
    };
}

#endif // TYPES_H