#include <string>
#include <cstdint>

#pragma once

#ifndef TYPES_H
#define TYPES_H

extern "C" {

    // block mesh grid structs
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
        float cfl, sfac, sfac_res, dt, d_max, d_avg;
        float d_var, facsec, fcorr;
        int nsteps, nstep, nkruts, guess_method, tstep_method;
        float ro_ref, roe_ref, rov_ref;
        int ni, nj;
        int nn, nm;
        uint8_t crashed;

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

    // curve fill mesh structs

    struct cell2d {
        double xmin;
        double xmax;
        double ymin;
        double ymax;
        int    level;
        int    id;
        int    neigh_offset;
        uint8_t neigh_count;
        uint8_t sc1, sc2, sc3, sc4;
    };

    struct lod_mesh {
        int length;
        cell2d* cells;
        int* neigh_indices;
    };
}

#endif // TYPES_H