#ifndef _DIM3_SWEEP_CUDA_H
#define _DIM3_SWEEP_CUDA_H

#include "snap_data.h"

int test();

void dim3_sweep_cuda ( input_data *input_vars, 
    bool firsty, bool lasty, bool firstz, bool lastz,
    geom_data *geom_vars, sn_data *sn_vars,
    data_data *data_vars, control_data *control_vars,
    solvar_data *solvar_vars, dim_sweep_data *dim_sweep_vars,
    int ich, int i_dir, int d1, int d2, int d3, int d4, int j_dir,
    int k_dir, int j_low, int k_low, int j_high, int k_high, int j_step, int k_step,
    int i1, int i2, int oct, int g, int *ierr );

#endif