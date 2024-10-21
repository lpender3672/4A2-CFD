
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

        // Logical array for wall locations
        bool* wall;
    };
}

#endif // TYPES_H