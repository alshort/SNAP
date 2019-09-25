#include <cuda.h>
#include <cuda_runtime.h>
#include <iostream>
extern "C" {
  #include "dim3_sweep_cuda.h"
}

// Local variable array macro
#define PSI_1D(ANG)   psi[ANG]
#define PC_1D(ANG)    pc[ANG]
#define DEN_1D(ANG)   den[ANG]

#ifdef ROWORDER
#define HV_2D(ANG, X) hv[ ANG*4                 \
                          + X ]
#else
#define HV_2D(ANG, X) hv[ X*NANG                \
                          + ANG ]
#endif

#ifdef ROWORDER
#define FXHV_2D(ANG, X) fxhv[ ANG*4             \
                              + X ]
#else
#define FXHV_2D(ANG, X) fxhv[ X*NANG            \
                              + ANG ]
#endif

// Simplify array indexing when certain values constant throughout module
#define PSII_3D(ANG, Y, Z)       PSII_4D(ANG, Y, Z, (g-1))
#define PSIJ_3D(ANG, CHUNK, Z)   PSIJ_4D(ANG, CHUNK, Z, (g-1))
#define PSIK_3D(ANG, CHUNK, Y)   PSIK_4D(ANG, CHUNK, Y, (g-1))
#define QTOT_4D(MOM1, X, Y, Z)   QTOT_5D(MOM1, X, Y, Z, (g-1))
#define EC_2D(ANG, MOM1)         EC_3D(ANG, MOM1, (oct-1))
#define VDELT_CONST              VDELT_1D(g-1)
#define PTR_IN_4D(ANG, X, Y, Z)  PTR_IN_6D(ANG, X, Y, Z, (i1-1), (i2-1))
#define PTR_OUT_4D(ANG, X, Y, Z) PTR_OUT_6D(ANG, X, Y, Z, (i1-1), (i2-1))
#define DINV_4D(ANG, X, Y, Z)    DINV_5D(ANG, X, Y, Z, (g-1))
#define FLUX_3D(X, Y, Z)         FLUX_4D(X, Y, Z, (g-1))
#define FLUXM_4D(MOM1, X, Y, Z)  FLUXM_5D(MOM1, X, Y, Z, (g-1))
#define JB_IN_3D(ANG, CHUNK, Z)  JB_IN_4D(ANG, CHUNK, Z, (g-1))
#define JB_OUT_3D(ANG, CHUNK, Z) JB_OUT_4D(ANG, CHUNK, Z, (g-1))
#define KB_IN_3D(ANG, CHUNK, Y)  KB_IN_4D(ANG, CHUNK, Y, (g-1))
#define KB_OUT_3D(ANG, CHUNK, Y) KB_OUT_4D(ANG, CHUNK, Y, (g-1))
#define FLKX_3D(X, Y, Z)         FLKX_4D(X, Y, Z, (g-1))
#define FLKY_3D(X, Y, Z)         FLKY_4D(X, Y, Z, (g-1))
#define FLKZ_3D(X, Y, Z)         FLKZ_4D(X, Y, Z, (g-1))
#define T_XS_3D(X, Y, Z)         T_XS_4D(X, Y, Z, (g-1))


// CUDA vars
#define N   10



__global__ void add( int *a, int *b, int *c )
{
    int tid = blockIdx.x;    // this thread handles the data at its thread id
    if (tid < N)
        c[tid] = a[tid] + b[tid];
}

int test( void )
{
    int a[N], b[N], c[N];
    int *dev_a, *dev_b, *dev_c;

    // allocate the memory on the GPU
    cudaMalloc( (void**)&dev_a, N * sizeof(int) );
    cudaMalloc( (void**)&dev_b, N * sizeof(int) );
    cudaMalloc( (void**)&dev_c, N * sizeof(int) );

    // fill the arrays 'a' and 'b' on the CPU
    for (int i=0; i<N; i++) {
        a[i] = -i;
        b[i] = i * i;
    }

    // copy the arrays 'a' and 'b' to the GPU
    cudaMemcpy( dev_a, a, N * sizeof(int), cudaMemcpyHostToDevice );
    cudaMemcpy( dev_b, b, N * sizeof(int), cudaMemcpyHostToDevice );

    add<<<N,1>>>( dev_a, dev_b, dev_c );

    // copy the array 'c' back from the GPU to the CPU
    cudaMemcpy( c, dev_c, N * sizeof(int), cudaMemcpyDeviceToHost );

    // display the results
    for (int i=0; i<N; i++) {
        printf( "%d + %d = %d\n", a[i], b[i], c[i] );
    }

    // free the memory allocated on the GPU
    cudaFree( dev_a );
    cudaFree( dev_b );
    cudaFree( dev_c );

    return 0;
}

__global__
void diagonal_loop( input_data *input_vars, 
    bool firsty, bool lasty, bool firstz, bool lastz,
    geom_data *geom_vars, sn_data *sn_vars,
    data_data *data_vars, control_data *control_vars,
    solvar_data *solvar_vars, dim_sweep_data *dim_sweep_vars,

    int ich, int i_dir, int d1, int d2, int d3, int d4, int j_dir,
    int k_dir, int j_low, int k_low, int j_high, int k_high, int j_step, int k_step,
    int i1, int i2, int oct, int g, int *ierr,

    int nang, 
    double *psi, double *pc, double *den, double *hv, double *fxhv,

    int d)
{
    // Local variables
    int n = threadIdx.x;
    int i_step, ic, i, j, k, l, ibl, ibb, ibf;

    int ang, indx1 = 4;

    double sum_hv = 0, sum_hv_tmp = 0, sum_wpsi = 0, sum_ecwpsi = 0,
        sum_wmupsii = 0, sum_wetapsij = 0, sum_wxipsik = 0;

    // Set up the sweep order in the i-direction.
    i_step = -1;
    if ( i_dir == 2 ) i_step = 1;

    // Loop over cells along the diagonals. When only 1 diagonal, it's
    // normal sweep order. Otherwise, nested threading performs mini-KBA.
    ic = DIAG_1D(d-1).cell_id_vars[n-1].ic;

    if ( i_step < 0 )
    {
        i = ich*ICHUNK - ic + 1;
    }
    else
    {
        i = (ich-1)*ICHUNK + ic;
    }

    if ( i <= NX )
    {
        j = DIAG_1D(d-1).cell_id_vars[n-1].jc;

        if ( j_step < 0 )
        {
            j = NY - j + 1;
        }

        k = DIAG_1D(d-1).cell_id_vars[n-1].kc;

        if ( k_step < 0 )
        {
            k = NZ - k + 1;
        }

        // Left/right boundary conditions, always vacuum.
        ibl = 0;

        if ( (i == NX) && (i_step == -1) )
        {
            for ( ang = 0; ang < nang; ang++ )
            {
                PSII_3D(ang,(j-1),(k-1)) = 0;
            }
        }
        else if ( i == 1 && i_step == 1 )
        {
            switch ( ibl )
            {
            case 0:
            case 1:
                for ( ang = 0; ang < nang; ang++ )
                {
                    PSII_3D(ang,(j-1),(k-1)) = 0;
                }
            }
        }

        // Top/bottom boundary condtions. Vacuum at global boundaries,
        // but set to some incoming flux from neighboring proc.
        ibb = 0;
        
        if ( j == j_low )
        {
            if ( j_dir == 1 && lasty )
            {
                for ( ang = 0; ang < nang; ang++ )
                {
                    PSIJ_3D(ang,(ic-1),(k-1)) = 0;
                }
            }
            else if ( j_dir == 2 && firsty )
            {
                switch ( ibb )
                {
                case 0:
                case 1:
                    for ( ang = 0; ang < nang; ang++ )
                    {
                        PSIJ_3D(ang,(ic-1),(k-1)) = 0;
                    }
                }
            }
            else
            {
                for ( ang = 0; ang < nang; ang++ )
                {
                    PSIJ_3D(ang,(ic-1),(k-1))
                        = JB_IN_3D(ang,(ic-1),(k-1));
                }
            }
        }

        // Front/back boundary condtions. Vacuum at global boundaries, 
        // but set to some incoming flux from neighboring proc.
        ibf = 0;
        
        if ( k == k_low )
        {
            if ( (k_dir == 1 && lastz) || NDIMEN < 3 )
            {
                for ( ang = 0; ang < nang; ang++ )
                {
                    PSIK_3D(ang,(ic-1),(j-1)) = 0;
                }
            }
            else if ( k_dir == 2 && firstz )
            {
                switch ( ibf )
                {
                case 0:
                case 1:
                    for ( ang = 0; ang < nang; ang++ )
                    {
                        PSIK_3D(ang,(ic-1),(j-1)) = 0;
                    }
                }
            }
            else
            {
                for ( ang = 0; ang < nang; ang++ )
                {
                    PSIK_3D(ang,(ic-1),(j-1))
                        = KB_IN_3D(ang,(ic-1),(j-1));
                }
            }
        }

        // Compute the angular source
        for ( ang = 0; ang < nang; ang++ )
        {
            PSI_1D(ang) = QTOT_4D(0,(i-1),(j-1),(k-1));

            if ( SRC_OPT == 3 )
            {
                PSI_1D(ang) +=
                    QIM_6D(ang,(i-1),(j-1),(k-1),(oct-1),(g-1));
            }
        }

        for ( l = 2; l <= CMOM; l++ )
        {
            for ( ang = 0; ang < nang; ang++ )
            {
                PSI_1D(ang) +=
                    EC_2D(ang,(l-1))
                    *QTOT_4D((l-1),(i-1),(j-1),(k-1));
            }
        }

        // Compute the numerator for the update formula
        for ( ang = 0; ang < nang; ang++ )
        {
            PC_1D(ang) = PSI_1D(ang)
                + PSII_3D(ang,(j-1),(k-1)) *MU_1D(ang)*HI
                + PSIJ_3D(ang,(ic-1),(k-1))*HJ_1D(ang)
                + PSIK_3D(ang,(ic-1),(j-1))*HK_1D(ang);

            if ( VDELT_CONST != 0 )
            {
                PC_1D(ang) += VDELT_CONST
                    *PTR_IN_6D(ang,(i-1),(j-1),(k-1),(i1-1),(i2-1));
            }
        }

        // Compute the solution of the center. Use DD for edges. 
        // Use fixup if requested.
        if ( FIXUP == 0 )
        {
            for ( ang = 0; ang < nang; ang++ )
            {
                PSI_1D(ang)
                    = PC_1D(ang)*DINV_4D(ang,(i-1),(j-1),(k-1));

                PSII_3D(ang,(j-1),(k-1))
                    = 2*PSI_1D(ang) - PSII_3D(ang,(j-1),(k-1));

                PSIJ_3D(ang,(ic-1),(k-1))
                    = 2*PSI_1D(ang) - PSIJ_3D(ang,(ic-1),(k-1));

                if ( NDIMEN == 3 )
                {
                    PSIK_3D(ang,(ic-1),(j-1))
                        = 2*PSI_1D(ang) - PSIK_3D(ang,(ic-1),(j-1));
                }

                if ( VDELT_CONST != 0 )
                {
                    PTR_OUT_6D(ang,(i-1),(j-1),(k-1),(i1-1),(i2-1))
                        = 2*PSI_1D(ang)
                        - PTR_IN_6D(ang,(i-1),(j-1),(k-1),(i1-1),(i2-1));
                }
            }
        }
        else
        {
            // Multi-pass set to zero + rebalance fixup. Determine angles
            // that will need fixup first.
            sum_hv = 0;
            for (ang = 0; ang < nang; ang++)
            {
                for (indx1 = 0; indx1 < 4; indx1++)
                {
                    HV_2D(ang, indx1) = 1;
                    sum_hv += HV_2D(ang,indx1);
                }

                PC_1D(ang) = PC_1D(ang) * DINV_4D(ang,(i-1),(j-1),(k-1));
            }

            // fixup_loop
            while (true)
            {
                sum_hv_tmp = 0;

                for ( ang = 0; ang < nang; ang++ )
                {
                    FXHV_2D(ang,0) =  2*PC_1D(ang) - PSII_3D(ang,(j-1),(k-1));

                    FXHV_2D(ang,1) =  2*PC_1D(ang) - PSIJ_3D(ang,(ic-1),(k-1));

                    if ( NDIMEN == 3 )
                    {
                        FXHV_2D(ang,2) = 2*PC_1D(ang) - PSIK_3D(ang,(ic-1),(j-1));
                    }

                    if ( VDELT_CONST != 0 )
                    {
                        FXHV_2D(ang,3) = 2*PC_1D(ang) - PTR_IN_6D(ang,(i-1),(j-1),(k-1),(i1-1),(i2-1));
                    }

                    for ( indx1 = 0; indx1 < 4; indx1++ )
                    {
                        if ( FXHV_2D(ang,indx1) < 0 )
                        {
                            HV_2D(ang,indx1) = 0;
                        }
                        sum_hv_tmp += HV_2D(ang,indx1);
                    }
                }

                // Exit loop when all angles are fixed up
                if (sum_hv == sum_hv_tmp) break;

                sum_hv = sum_hv_tmp;

                // Recompute balance equation numerator and denominator 
                // and get new cell average flux
                for ( ang = 0; ang < nang; ang++ )
                {
                    PC_1D(ang) = PSII_3D(ang,(j-1),(k-1))
                        * MU_1D(ang) * HI * (1+HV_2D(ang,0))
                        + PSIJ_3D(ang,(ic-1),(k-1))
                        * HJ_1D(ang) * (1+HV_2D(ang,1))
                        + PSIK_3D(ang,(ic-1),(j-1))
                        * HK_1D(ang) * (1+HV_2D(ang,2));

                    if ( VDELT_CONST != 0 )
                    {
                        PC_1D(ang) += VDELT_CONST
                            * PTR_IN_6D(ang,(i-1),(j-1),(k-1),(i1-1),(i2-1))
                            * (1+HV_2D(ang,3));
                    }

                    PC_1D(ang) = PSI_1D(ang) + 0.5*PC_1D(ang);

                    DEN_1D(ang) = T_XS_3D((i-1),(j-1),(k-1))
                        + MU_1D(ang)  * HI * HV_2D(ang,0)
                        + HJ_1D(ang)  * HV_2D(ang,1)
                        + HK_1D(ang)  * HV_2D(ang,2)
                        + VDELT_CONST * HV_2D(ang,3);

                    if ( DEN_1D(ang) > TOLR )
                    {
                        PC_1D(ang) /= DEN_1D(ang);
                    }
                    else
                    {
                        PC_1D(ang) = 0;
                    }
                }
            } // end fixup_loop

            // Fixup done, compute edges
            for (ang = 0; ang < nang; ang++)
            {
                PSI_1D(ang) = PC_1D(ang);

                PSII_3D(ang,(j-1),(k-1)) = FXHV_2D(ang,0) * HV_2D(ang,0);

                PSIJ_3D(ang,(ic-1),(k-1)) = FXHV_2D(ang,1) * HV_2D(ang,1);

                if (NDIMEN == 3)
                {
                    PSIK_3D(ang,(ic-1),(j-1)) = FXHV_2D(ang,2) * HV_2D(ang,2);
                }

                if (VDELT_CONST != 0)
                {
                    PTR_OUT_6D(ang,(i-1),(j-1),(k-1),(i1-1),(i2-1)) = FXHV_2D(ang,3) * HV_2D(ang,3);
                }
            }
        }

        // Clear the flux arrays
        if ( oct == 1 )
        {
            FLUX_4D((i-1),(j-1),(k-1),(g-1)) = 0;

            for ( indx1 = 0; indx1 < (CMOM-1); indx1++ )
            {
                FLUXM_5D(indx1,(i-1),(j-1),(k-1),(g-1)) = 0;
            }
        }

        // Compute the flux moments
        sum_wpsi = 0;

        for (ang = 0; ang < nang; ang++)
        {
            sum_wpsi += W_1D(ang) * PSI_1D(ang);
        }

        FLUX_4D((i-1),(j-1),(k-1),(g-1)) += sum_wpsi;

        for (l = 1; l <= (CMOM-1); l++)
        {
            sum_ecwpsi = 0;

            for ( ang = 0; ang < nang; ang++ )
            {
                sum_ecwpsi += EC_2D(ang,(l))*W_1D(ang)*PSI_1D(ang);
            }

            FLUXM_5D((l-1),(i-1),(j-1),(k-1),(g-1)) += sum_ecwpsi;
        }

        // Calculate min and max scalar fluxes (not used elsewhere currently)
        if (oct == NOCT)
        {
            FMIN = MIN( FMIN, FLUX_3D((i-1),(j-1),(k-1)) );
            FMAX = MAX( FMAX, FLUX_3D((i-1),(j-1),(k-1)) );
        }

        // Save edge fluxes (dummy if checks for unused non-vacuum BCs)
        if (j == j_high)
        {
            if ((j_dir==2 && lasty) ||
                ((j_dir == 1 && firsty) && ibb == 1))
            {
                // CONTINUE
            }
            else
            {
                for (ang = 0; ang < nang; ang++)
                {
                    JB_OUT_3D(ang,(ic-1),(k-1)) = PSIJ_3D(ang,(ic-1),(k-1));
                }
            }
        }
        
        if (k == k_high)
        {
            if ((k_dir == 2 && lastz) ||
                ((k_dir==1 && firstz) && ibf == 1))
            {
                // CONTINUE
            }
            else
            {
                for ( ang = 0; ang < nang; ang++ )
                {
                    KB_OUT_3D(ang,(ic-1),(j-1)) = PSIK_3D(ang,(ic-1),(j-1));
                }
            }
        }

        // Compute leakages (not used elsewhere currently)
        if (((i+i_dir-1) == 1) || ((i+i_dir-1) == (NX+1)))
        {
            sum_wmupsii = 0;

            for (ang = 0; ang < nang; ang++)
            {
                sum_wmupsii += WMU_1D(ang) * PSII_3D(ang,(j-1),(k-1));
            }

            FLKX_3D((i+i_dir-1-1),(j-1),(k-1)) += i_step*sum_wmupsii;
        }
        if ((j_dir == 1 && firsty) || (j_dir == 2 && lasty))
        {
            sum_wetapsij = 0;

            for (ang = 0; ang < nang; ang++)
            {
                sum_wetapsij
                    += WETA_1D(ang) * PSIJ_3D(ang,(ic-1),(k-1));
            }

            FLKY_3D((i-1),(j+j_dir-1-1),(k-1)) += j_step*sum_wetapsij;
        }

        if (((k_dir == 1 && firstz) || (k_dir == 2 && lastz)) && NDIMEN == 3)
        {
            sum_wxipsik = 0;

            for (ang = 0; ang < nang; ang++)
            {
                sum_wxipsik += WXI_1D(ang) * PSIK_3D(ang,(ic-1),(j-1));
            }

            FLKZ_3D((i-1),(j-1),(k+k_dir-1-1))
                += k_step*sum_wxipsik;
        }
    }
}

void dim3_sweep_cuda ( input_data *input_vars, 
    bool firsty, bool lasty, bool firstz, bool lastz,
    geom_data *geom_vars, sn_data *sn_vars,
    data_data *data_vars, control_data *control_vars,
    solvar_data *solvar_vars, dim_sweep_data *dim_sweep_vars,
    int ich, int i_dir, int d1, int d2, int d3, int d4, int j_dir,
    int k_dir, int j_low, int k_low, int j_high, int k_high, int j_step, int k_step,
    int i1, int i2, int oct, int g, int *ierr )
{
    // Local variables
    int i;
    int ang, y_ind, ic_ind, z_ind = 4;

    double psi[NANG], pc[NANG], den[NANG];
    double hv[NANG*4], fxhv[NANG*4];

    double *c_psi[NANG], *c_pc[NANG], *c_den[NANG];
    double *c_hv[NANG*4], *c_fxhv[NANG*4];
    
    // Create GPU-copies of data
    input_data *c_input_vars;
    geom_data *c_geom_vars;
    sn_data *c_sn_vars;
    data_data *c_data_vars;
    control_data *c_control_vars;
    solvar_data *c_solvar_vars;
    dim_sweep_data *c_dim_sweep_vars;

    cudaMalloc(&c_input_vars, sizeof(input_data));
    cudaMalloc(&c_geom_vars, sizeof(geom_data));
    cudaMalloc(&c_sn_vars, sizeof(sn_data));
    cudaMalloc(&c_data_vars, sizeof(data_data));
    cudaMalloc(&c_control_vars, sizeof(control_data));
    cudaMalloc(&c_solvar_vars, sizeof(solvar_data));
    cudaMalloc(&c_dim_sweep_vars, sizeof(dim_sweep_data));

    cudaMalloc(c_psi, NANG * sizeof(double));
    cudaMalloc(c_pc, NANG * sizeof(double));
    cudaMalloc(c_den, NANG * sizeof(double));
    cudaMalloc(c_hv, NANG * 4 * sizeof(double));
    cudaMalloc(c_fxhv, NANG * 4 * sizeof(double));


    // Zero out the outgoing boundary arrays and fixup array
    for ( z_ind = 0; z_ind < NZ; z_ind++ )
    {
        for ( ic_ind = 0; ic_ind < ICHUNK; ic_ind++ )
        {
            for ( ang = 0; ang < NANG; ang++ )
            {
                JB_OUT_3D(ang,ic_ind,z_ind) = 0;
            }
        }
    }

    for ( y_ind = 0; y_ind < NY; y_ind++ )
    {
        for ( ic_ind = 0; ic_ind < ICHUNK; ic_ind++ )
        {
            for ( ang = 0; ang < NANG; ang++ )
            {
                KB_OUT_3D(ang,ic_ind,y_ind) = 0;
            }
        }
    }

    for ( i = 0; i < 4; i++)
    {
        for ( ang = 0; ang < NANG; ang++ )
        {
            FXHV_2D(ang, i) = 0;
        }
    }

    
    cudaMemcpy(c_psi, psi, NANG * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(c_pc, pc, NANG * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(c_den, den, NANG * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(c_hv, hv, NANG * 4 * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(c_fxhv, fxhv, NANG * 4 * sizeof(double), cudaMemcpyHostToDevice);

    cudaMemcpy(c_input_vars, input_vars, sizeof(input_data), cudaMemcpyHostToDevice);
    cudaMemcpy(c_geom_vars, geom_vars, sizeof(geom_data), cudaMemcpyHostToDevice);
    cudaMemcpy(c_sn_vars, sn_vars, sizeof(sn_data), cudaMemcpyHostToDevice);
    cudaMemcpy(c_data_vars, data_vars, sizeof(data_data), cudaMemcpyHostToDevice);
    cudaMemcpy(c_control_vars, control_vars, sizeof(control_data), cudaMemcpyHostToDevice);
    cudaMemcpy(c_solvar_vars, solvar_vars, sizeof(solvar_data), cudaMemcpyHostToDevice);
    cudaMemcpy(c_dim_sweep_vars, dim_sweep_vars, sizeof(dim_sweep_data), cudaMemcpyHostToDevice);

    // Loop over cells along the diagonals. When only 1 diagonal, it's
    // normal sweep order. Otherwise, nested threading performs mini-KBA.
    // diagonal loop
    int d;
    for (d = 1; d <= NDIAG; d++)
    {
        printf("lenc: %d\n", DIAG_1D(d-1).lenc);
        diagonal_loop<<<1, (DIAG_1D(d-1).lenc)>>>(
            c_input_vars,
            firsty, lasty, firstz, lastz,
            c_geom_vars, c_sn_vars,
            c_data_vars, c_control_vars,
            c_solvar_vars, c_dim_sweep_vars,

            ich, i_dir, d1, d2, d3, d4, j_dir,
            k_dir, j_low, k_low, j_high, k_high, j_step, k_step,
            i1, i2, oct, g, ierr,

            input_vars->nang,
            *c_psi, *c_pc, *c_den, *c_hv, *c_fxhv,
            d
        );
    }

    // Copy from device back to hsot
    cudaMemcpy(c_input_vars, input_vars, sizeof(input_data), cudaMemcpyDeviceToHost);
    cudaMemcpy(c_geom_vars, geom_vars, sizeof(geom_data), cudaMemcpyDeviceToHost);
    cudaMemcpy(c_sn_vars, sn_vars, sizeof(sn_data), cudaMemcpyDeviceToHost);
    cudaMemcpy(c_data_vars, data_vars, sizeof(data_data), cudaMemcpyDeviceToHost);
    cudaMemcpy(c_control_vars, control_vars, sizeof(control_data), cudaMemcpyDeviceToHost);
    cudaMemcpy(c_solvar_vars, solvar_vars, sizeof(solvar_data), cudaMemcpyDeviceToHost);
    cudaMemcpy(c_dim_sweep_vars, dim_sweep_vars, sizeof(dim_sweep_data), cudaMemcpyDeviceToHost);

    // Clean up
    cudaFree(c_input_vars);
    cudaFree(c_geom_vars);
    cudaFree(c_sn_vars);
    cudaFree(c_data_vars);
    cudaFree(c_control_vars);
    cudaFree(c_solvar_vars);
    cudaFree(c_dim_sweep_vars);

    cudaFree(c_psi);
    cudaFree(c_pc);
    cudaFree(c_den);
    cudaFree(c_hv);
    cudaFree(c_fxhv);
}