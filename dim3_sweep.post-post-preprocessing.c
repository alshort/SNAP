void dim3_sweep_data_init ( dim_sweep_data *dim_sweep_vars )
{
    dim_sweep_vars->fmin = 0;
    dim_sweep_vars->fmax = 0;
}

// 3-D slab mesh sweeper.
void dim3_sweep ( input_data *input_vars, para_data *para_vars,
                  geom_data *geom_vars, sn_data *sn_vars,
                  data_data *data_vars, control_data *control_vars,
                  solvar_data *solvar_vars, dim_sweep_data *dim_sweep_vars,
                  int ich, int i_dir, int d1, int d2, int d3, int d4, int j_dir,
                  int k_dir, int j_low, int k_low, int j_high, int k_high, int j_step, int k_step,
                  int i1, int i2, int oct, int g, int *ierr )
{
    // Local variables 
    int i_step, d, n, ic, i, j, k, l, ibl, ibr, ibb, ibt, ibf, ibk;

    int z_ind, y_ind, ic_ind, ang, indx1 = 4;

    double sum_hv = 0, sum_hv_tmp = 0, sum_wpsi = 0, sum_ecwpsi = 0,
        sum_wmupsii = 0, sum_wetapsij = 0, sum_wxipsik = 0;

    double psi[input_vars->nang], pc[input_vars->nang], den[input_vars->nang];
    double hv[input_vars->nang*4], fxhv[input_vars->nang*4];

    double vec1_vec2_tmp[input_vars->nang], PSI_2X[input_vars->nang], hv_p1[input_vars->nang*4],
        mu_hv[input_vars->nang], hj_hv[input_vars->nang], hk_hv[input_vars->nang], w_psi[input_vars->nang];

    double unit_vec[indx1];
    for ( i = 0; i < indx1; i++ )
    {
        unit_vec[i] = 1;
    }
 
    // Set up the sweep order in the i-direction.
    i_step = -1;
    if ( i_dir == 2 ) i_step = 1;

    // Zero out the outgoing boundary arrays and fixup array
    for ( z_ind = 0; z_ind < input_vars->nz; z_ind++ )
    {
        for ( ic_ind = 0; ic_ind < input_vars->ichunk; ic_ind++ )
        {
            for ( ang = 0; ang < input_vars->nang; ang++ )
            {
                solvar_vars->jb_out[ (g-1) * input_vars->nz * input_vars->ichunk * input_vars->nang + z_ind * input_vars->ichunk * input_vars->nang + ic_ind * input_vars->nang + ang ] = 0;
            }
        }
    }

    for ( y_ind = 0; y_ind < input_vars->ny; y_ind++ )
    {
        for ( ic_ind = 0; ic_ind < input_vars->ichunk; ic_ind++ )
        {
            for ( ang = 0; ang < input_vars->nang; ang++ )
            {
                solvar_vars->kb_out[ (g-1) * input_vars->ny * input_vars->ichunk * input_vars->nang + y_ind * input_vars->ichunk * input_vars->nang + ic_ind * input_vars->nang + ang ] = 0;
            }
        }
    }

    for ( i = 0; i < 4; i++)
    {
        for ( ang = 0; ang < input_vars->nang; ang++ )
        {
            fxhv[ i*input_vars->nang + ang ] = 0;
        }
    }

    // Loop over cells along the diagonals. When only 1 diagonal, it's
    // normal sweep order. Otherwise, nested threading performs mini-KBA.
    for ( d = 1; d <= geom_vars->ndiag; d++ ) // diagonal loop
    {
        #pragma omp for schedule(static, 1) private(n,ic,i,j,k,l,psi,pc,sum_hv,hv,den)        
        for ( n = 1; n <= (geom_vars->diag_vars[d-1].lenc); n++ ) // line_loop
        {
            ic = geom_vars->diag_vars[d-1].cell_id_vars[n-1].ic;

            if ( i_step < 0 )
            {
                i = ich*input_vars->ichunk - ic + 1;
            }
            else
            {
                i = (ich-1)*input_vars->ichunk + ic;
            }

            if ( i <= input_vars->nx )
            {
                j = geom_vars->diag_vars[d-1].cell_id_vars[n-1].jc;

                if ( j_step < 0 )
                {
                    j = input_vars->ny - j + 1;
                }

                k = geom_vars->diag_vars[d-1].cell_id_vars[n-1].kc;

                if ( k_step < 0 )
                {
                    k = input_vars->nz - k + 1;
                }

                // Left/right boundary conditions, always vacuum.
                ibl = 0;
                ibr = 0;

                if ( (i == input_vars->nx) && (i_step == -1) )
                {
                    for ( ang = 0; ang < input_vars->nang; ang++ )
                    {
                        solvar_vars->psii[ (g-1) * input_vars->nz * input_vars->ny * input_vars->nang + (k-1) * input_vars->ny * input_vars->nang + (j-1) * input_vars->nang + ang ] = 0;
                    }
                }
                else if ( i == 1 && i_step == 1 )
                {
                    switch ( ibl )
                    {
                    case 0:
                    case 1:
                        for ( ang = 0; ang < input_vars->nang; ang++ )
                        {
                            solvar_vars->psii[ (g-1) * input_vars->nz * input_vars->ny * input_vars->nang + (k-1) * input_vars->ny * input_vars->nang + (j-1) * input_vars->nang + ang ] = 0;
                        }
                    }
                }

                // Top/bottom boundary condtions. Vacuum at global boundaries, but
                // set to some incoming flux from neighboring proc. 
                ibb = 0;
                ibt = 0;
                if ( j == j_low )
                {
                    if ( j_dir == 1 && para_vars->lasty )
                    {
                        for ( ang = 0; ang < input_vars->nang; ang++ )
                        {
                            solvar_vars->psij[ (g-1) * input_vars->nz * input_vars->ichunk * input_vars->nang + (k-1) * input_vars->ichunk * input_vars->nang + (ic-1) * input_vars->nang + ang ] = 0;
                        }
                    }
                    else if ( j_dir == 2 && para_vars->firsty )
                    {
                        switch ( ibb )
                        {
                        case 0:
                        case 1:
                            for ( ang = 0; ang < input_vars->nang; ang++ )
                            {
                                solvar_vars->psij[ (g-1) * input_vars->nz * input_vars->ichunk * input_vars->nang + (k-1) * input_vars->ichunk * input_vars->nang + (ic-1) * input_vars->nang + ang ] = 0;
                            }
                        }
                    }
                    else
                    {
                        for ( ang = 0; ang < input_vars->nang; ang++ )
                        {
                            solvar_vars->psij[ (g-1) * input_vars->nz * input_vars->ichunk * input_vars->nang + (k-1) * input_vars->ichunk * input_vars->nang + (ic-1) * input_vars->nang + ang ]
                                = solvar_vars->jb_in[ (g-1) * input_vars->nz * input_vars->ichunk * input_vars->nang + (k-1) * input_vars->ichunk * input_vars->nang + (ic-1) * input_vars->nang + ang ];
                        }
                    }
                }
 
                // Front/back boundary condtions. Vacuum at global boundaries, but
                // set to some incoming flux from neighboring proc.
                ibf = 0;
                ibk = 0;
                if ( k == k_low )
                {
                    if ( (k_dir == 1 && para_vars->lastz) || input_vars->ndimen < 3 )
                    {
                        for ( ang = 0; ang < input_vars->nang; ang++ )
                        {
                            solvar_vars->psik[ (g-1) * input_vars->ny * input_vars->ichunk * input_vars->nang + (j-1) * input_vars->ichunk * input_vars->nang + (ic-1) * input_vars->nang + ang ] = 0;
                        }
                    }
                    else if ( k_dir == 2 && para_vars->firstz )
                    {
                        switch ( ibf )
                        {
                        case 0:
                        case 1:
                            for ( ang = 0; ang < input_vars->nang; ang++ )
                            {
                                solvar_vars->psik[ (g-1) * input_vars->ny * input_vars->ichunk * input_vars->nang + (j-1) * input_vars->ichunk * input_vars->nang + (ic-1) * input_vars->nang + ang ] = 0;
                            }
                        }
                    }

                    else
                    {
                        for ( ang = 0; ang < input_vars->nang; ang++ )
                        {
                            solvar_vars->psik[ (g-1) * input_vars->ny * input_vars->ichunk * input_vars->nang + (j-1) * input_vars->ichunk * input_vars->nang + (ic-1) * input_vars->nang + ang ]
                                = solvar_vars->kb_in[ (g-1) * input_vars->ny * input_vars->ichunk * input_vars->nang + (j-1) * input_vars->ichunk * input_vars->nang + (ic-1) * input_vars->nang + ang ];
                        }
                    }
                }

                // Compute the angular source
                for ( ang = 0; ang < input_vars->nang; ang++ )
                {
                    psi[ang] = solvar_vars->qtot[ (g-1) * input_vars->nz * input_vars->ny * input_vars->nx * sn_vars->cmom + (k-1) * input_vars->ny * input_vars->nx * sn_vars->cmom + (j-1) * input_vars->nx * sn_vars->cmom + (i-1) * sn_vars->cmom + 0 ];

                    if ( input_vars->src_opt == 3 )
                    {
                        psi[ang] +=
                            data_vars->qim[ (g-1) * sn_vars->noct * input_vars->nz * input_vars->ny * input_vars->nx * input_vars->nang + (oct-1) * input_vars->nz * input_vars->ny * input_vars->nx * input_vars->nang + (k-1) * input_vars->ny * input_vars->nx * input_vars->nang + (j-1) * input_vars->nx * input_vars->nang + (i-1) * input_vars->nang + ang ];
                    }
                }

                for ( l = 2; l <= sn_vars->cmom; l++ )
                {
                    for ( ang = 0; ang < input_vars->nang; ang++ )
                    {
                        psi[ang] +=
                            sn_vars->ec[ (oct-1) * sn_vars->cmom * input_vars->nang + (l-1) * input_vars->nang + ang ]
                            *solvar_vars->qtot[ (g-1) * input_vars->nz * input_vars->ny * input_vars->nx * sn_vars->cmom + (k-1) * input_vars->ny * input_vars->nx * sn_vars->cmom + (j-1) * input_vars->nx * sn_vars->cmom + (i-1) * sn_vars->cmom + (l-1) ];
                    }
                }

                // Compute the numerator for the update formula
                for ( ang = 0; ang < input_vars->nang; ang++ )
                {
                    pc[ang] = psi[ang]
                        + solvar_vars->psii[ (g-1) * input_vars->nz * input_vars->ny * input_vars->nang + (k-1) * input_vars->ny * input_vars->nang + (j-1) * input_vars->nang + ang ] *sn_vars->mu[ang]*geom_vars->hi
                        + solvar_vars->psij[ (g-1) * input_vars->nz * input_vars->ichunk * input_vars->nang + (k-1) * input_vars->ichunk * input_vars->nang + (ic-1) * input_vars->nang + ang ]*geom_vars->hj[ang]
                        + solvar_vars->psik[ (g-1) * input_vars->ny * input_vars->ichunk * input_vars->nang + (j-1) * input_vars->ichunk * input_vars->nang + (ic-1) * input_vars->nang + ang ]*geom_vars->hk[ang];

                    if ( data_vars->vdelt[g-1] != 0 )
                    {
                        pc[ang] += data_vars->vdelt[g-1]
                            *solvar_vars->ptr_in[ (i2-1) * sn_vars->noct * input_vars->nz * input_vars->ny * input_vars->nx * input_vars->nang + (i1-1) * input_vars->nz * input_vars->ny * input_vars->nx * input_vars->nang + (k-1) * input_vars->ny * input_vars->nx * input_vars->nang + (j-1) * input_vars->nx * input_vars->nang + (i-1) * input_vars->nang + ang ];
                    }
                }

                // Compute the solution of the center. Use DD for edges. Use fixup
                // if requested.
                if ( input_vars->fixup == 0 )
                {
                    for ( ang = 0; ang < input_vars->nang; ang++ )
                    {
                        psi[ang]
                            = pc[ang]*geom_vars->dinv[ (g-1) * input_vars->nz * input_vars->ny * input_vars->nx * input_vars->nang + (k-1) * input_vars->ny * input_vars->nx * input_vars->nang + (j-1) * input_vars->nx * input_vars->nang + (i-1) * input_vars->nang + ang ];

                        solvar_vars->psii[ (g-1) * input_vars->nz * input_vars->ny * input_vars->nang + (k-1) * input_vars->ny * input_vars->nang + (j-1) * input_vars->nang + ang ]
                            = 2*psi[ang] - solvar_vars->psii[ (g-1) * input_vars->nz * input_vars->ny * input_vars->nang + (k-1) * input_vars->ny * input_vars->nang + (j-1) * input_vars->nang + ang ];

                        solvar_vars->psij[ (g-1) * input_vars->nz * input_vars->ichunk * input_vars->nang + (k-1) * input_vars->ichunk * input_vars->nang + (ic-1) * input_vars->nang + ang ]
                            = 2*psi[ang] - solvar_vars->psij[ (g-1) * input_vars->nz * input_vars->ichunk * input_vars->nang + (k-1) * input_vars->ichunk * input_vars->nang + (ic-1) * input_vars->nang + ang ];

                        if ( input_vars->ndimen == 3 )
                        {
                            solvar_vars->psik[ (g-1) * input_vars->ny * input_vars->ichunk * input_vars->nang + (j-1) * input_vars->ichunk * input_vars->nang + (ic-1) * input_vars->nang + ang ]
                                = 2*psi[ang] - solvar_vars->psik[ (g-1) * input_vars->ny * input_vars->ichunk * input_vars->nang + (j-1) * input_vars->ichunk * input_vars->nang + (ic-1) * input_vars->nang + ang ];
                        }

                        if ( data_vars->vdelt[g-1] != 0 )
                        {
                            solvar_vars->ptr_out[ (i2-1) * sn_vars->noct * input_vars->nz * input_vars->ny * input_vars->nx * input_vars->nang + (i1-1) * input_vars->nz * input_vars->ny * input_vars->nx * input_vars->nang + (k-1) * input_vars->ny * input_vars->nx * input_vars->nang + (j-1) * input_vars->nx * input_vars->nang + (i-1) * input_vars->nang + ang ]
                                = 2*psi[ang]
                                - solvar_vars->ptr_in[ (i2-1) * sn_vars->noct * input_vars->nz * input_vars->ny * input_vars->nx * input_vars->nang + (i1-1) * input_vars->nz * input_vars->ny * input_vars->nx * input_vars->nang + (k-1) * input_vars->ny * input_vars->nx * input_vars->nang + (j-1) * input_vars->nx * input_vars->nang + (i-1) * input_vars->nang + ang ];
                        }
                    }
                }
                else
                {
                    // Multi-pass set to zero + rebalance fixup. Determine angles
                    // that will need fixup first.
                    sum_hv = 0;
                    for (ang = 0; ang < input_vars->nang; ang++)
                    {
                        for (indx1 = 0; indx1 < 4; indx1++)
                        {
                            hv[ indx1*input_vars->nang + ang ] = 1;
                            sum_hv += hv[ indx1*input_vars->nang + ang ];
                        }

                        pc[ang] = pc[ang] * geom_vars->dinv[ (g-1) * input_vars->nz * input_vars->ny * input_vars->nx * input_vars->nang + (k-1) * input_vars->ny * input_vars->nx * input_vars->nang + (j-1) * input_vars->nx * input_vars->nang + (i-1) * input_vars->nang + ang ];
                    }

                    // fixup_loop
                    while (1)
                    {
                        sum_hv_tmp = 0;

                        for ( ang = 0; ang < input_vars->nang; ang++ )
                        {
                            fxhv[ 0*input_vars->nang + ang ] =  2*pc[ang] - solvar_vars->psii[ (g-1) * input_vars->nz * input_vars->ny * input_vars->nang + (k-1) * input_vars->ny * input_vars->nang + (j-1) * input_vars->nang + ang ];

                            fxhv[ 1*input_vars->nang + ang ] =  2*pc[ang] - solvar_vars->psij[ (g-1) * input_vars->nz * input_vars->ichunk * input_vars->nang + (k-1) * input_vars->ichunk * input_vars->nang + (ic-1) * input_vars->nang + ang ];

                            if ( input_vars->ndimen == 3 )
                            {
                                fxhv[ 2*input_vars->nang + ang ] = 2*pc[ang] - solvar_vars->psik[ (g-1) * input_vars->ny * input_vars->ichunk * input_vars->nang + (j-1) * input_vars->ichunk * input_vars->nang + (ic-1) * input_vars->nang + ang ];
                            }

                            if ( data_vars->vdelt[g-1] != 0 )
                            {
                                fxhv[ 3*input_vars->nang + ang ] = 2*pc[ang] - solvar_vars->ptr_in[ (i2-1) * sn_vars->noct * input_vars->nz * input_vars->ny * input_vars->nx * input_vars->nang + (i1-1) * input_vars->nz * input_vars->ny * input_vars->nx * input_vars->nang + (k-1) * input_vars->ny * input_vars->nx * input_vars->nang + (j-1) * input_vars->nx * input_vars->nang + (i-1) * input_vars->nang + ang ];
                            }

                            for ( indx1 = 0; indx1 < 4; indx1++ )
                            {
                                if ( fxhv[ indx1*input_vars->nang + ang ] < 0 )
                                {
                                    hv[ indx1*input_vars->nang + ang ] = 0;
                                }
                                sum_hv_tmp += hv[ indx1*input_vars->nang + ang ];
                            }
                        }
                        
                        // Exit loop when all angles are fixed up
                        if (sum_hv == sum_hv_tmp) break;

                        sum_hv = sum_hv_tmp;

                        // Recompute balance equation numerator and denominator and get
                        // new cell average flux
                        for ( ang = 0; ang < input_vars->nang; ang++ )
                        {
                            pc[ang] = solvar_vars->psii[ (g-1) * input_vars->nz * input_vars->ny * input_vars->nang + (k-1) * input_vars->ny * input_vars->nang + (j-1) * input_vars->nang + ang ]
                                * sn_vars->mu[ang] * geom_vars->hi * (1+hv[ 0*input_vars->nang + ang ])
                                + solvar_vars->psij[ (g-1) * input_vars->nz * input_vars->ichunk * input_vars->nang + (k-1) * input_vars->ichunk * input_vars->nang + (ic-1) * input_vars->nang + ang ]
                                * geom_vars->hj[ang] * (1+hv[ 1*input_vars->nang + ang ])
                                + solvar_vars->psik[ (g-1) * input_vars->ny * input_vars->ichunk * input_vars->nang + (j-1) * input_vars->ichunk * input_vars->nang + (ic-1) * input_vars->nang + ang ]
                                * geom_vars->hk[ang] * (1+hv[ 2*input_vars->nang + ang ]);

                            if ( data_vars->vdelt[g-1] != 0 )
                            {
                                pc[ang] += data_vars->vdelt[g-1]
                                    * solvar_vars->ptr_in[ (i2-1) * sn_vars->noct * input_vars->nz * input_vars->ny * input_vars->nx * input_vars->nang + (i1-1) * input_vars->nz * input_vars->ny * input_vars->nx * input_vars->nang + (k-1) * input_vars->ny * input_vars->nx * input_vars->nang + (j-1) * input_vars->nx * input_vars->nang + (i-1) * input_vars->nang + ang ]
                                    * (1+hv[ 3*input_vars->nang + ang ]);
                            }

                            pc[ang] = psi[ang] + 0.5*pc[ang];

                            den[ang] = solvar_vars->t_xs[ (g-1) * input_vars->nz * input_vars->ny * input_vars->nx + (k-1) * input_vars->ny * input_vars->nx + (j-1) * input_vars->nx + (i-1) ]
                                + sn_vars->mu[ang]  * geom_vars->hi * hv[ 0*input_vars->nang + ang ]
                                + geom_vars->hj[ang]  * hv[ 1*input_vars->nang + ang ]
                                + geom_vars->hk[ang]  * hv[ 2*input_vars->nang + ang ]
                                + data_vars->vdelt[g-1] * hv[ 3*input_vars->nang + ang ];

                            if ( den[ang] > control_vars->tolr )
                            {
                                pc[ang] /= den[ang];
                            }
                            else
                            {
                                pc[ang] = 0;
                            }
                        }
                    } // end fixup_loop
 
                    // Fixup done, compute edges
                    for (ang = 0; ang < input_vars->nang; ang++)
                    {
                        psi[ang] = pc[ang];

                        solvar_vars->psii[ (g-1) * input_vars->nz * input_vars->ny * input_vars->nang + (k-1) * input_vars->ny * input_vars->nang + (j-1) * input_vars->nang + ang ] = fxhv[ 0*input_vars->nang + ang ] * hv[ 0*input_vars->nang + ang ];

                        solvar_vars->psij[ (g-1) * input_vars->nz * input_vars->ichunk * input_vars->nang + (k-1) * input_vars->ichunk * input_vars->nang + (ic-1) * input_vars->nang + ang ] = fxhv[ 1*input_vars->nang + ang ] * hv[ 1*input_vars->nang + ang ];

                        if (input_vars->ndimen == 3)
                        {
                            solvar_vars->psik[ (g-1) * input_vars->ny * input_vars->ichunk * input_vars->nang + (j-1) * input_vars->ichunk * input_vars->nang + (ic-1) * input_vars->nang + ang ] = fxhv[ 2*input_vars->nang + ang ] * hv[ 2*input_vars->nang + ang ];
                        }

                        if (data_vars->vdelt[g-1] != 0)
                        {
                            solvar_vars->ptr_out[ (i2-1) * sn_vars->noct * input_vars->nz * input_vars->ny * input_vars->nx * input_vars->nang + (i1-1) * input_vars->nz * input_vars->ny * input_vars->nx * input_vars->nang + (k-1) * input_vars->ny * input_vars->nx * input_vars->nang + (j-1) * input_vars->nx * input_vars->nang + (i-1) * input_vars->nang + ang ] = fxhv[ 3*input_vars->nang + ang ] * hv[ 3*input_vars->nang + ang ];
                        }
                    }
                }

                // Clear the flux arrays
                if ( oct == 1 )
                {
                    solvar_vars->flux[ (g-1) * input_vars->nz * input_vars->ny * input_vars->nx + (k-1) * input_vars->ny * input_vars->nx + (j-1) * input_vars->nx + (i-1) ] = 0;

                    for ( indx1 = 0; indx1 < (sn_vars->cmom-1); indx1++ )
                    {
                        solvar_vars->fluxm[ (g-1) * input_vars->nz * input_vars->ny * input_vars->nx * (sn_vars->cmom-1) + (k-1) * input_vars->ny * input_vars->nx * (sn_vars->cmom-1) + (j-1) * input_vars->nx * (sn_vars->cmom-1) + (i-1) * (sn_vars->cmom-1) + indx1 ] = 0;
                    }
                }

                // Compute the flux moments
                sum_wpsi = 0;

                for (ang = 0; ang < input_vars->nang; ang++)
                {
                    sum_wpsi += sn_vars->w[ang] * psi[ang];
                }

                solvar_vars->flux[ (g-1) * input_vars->nz * input_vars->ny * input_vars->nx + (k-1) * input_vars->ny * input_vars->nx + (j-1) * input_vars->nx + (i-1) ] += sum_wpsi;

                for (l = 1; l <= (sn_vars->cmom-1); l++)
                {
                    sum_ecwpsi = 0;

                    for ( ang = 0; ang < input_vars->nang; ang++ )
                    {
                        sum_ecwpsi += sn_vars->ec[ (oct-1) * sn_vars->cmom * input_vars->nang + (l) * input_vars->nang + ang ]*sn_vars->w[ang]*psi[ang];
                    }

                    solvar_vars->fluxm[ (g-1) * input_vars->nz * input_vars->ny * input_vars->nx * (sn_vars->cmom-1) + (k-1) * input_vars->ny * input_vars->nx * (sn_vars->cmom-1) + (j-1) * input_vars->nx * (sn_vars->cmom-1) + (i-1) * (sn_vars->cmom-1) + (l-1) ] += sum_ecwpsi;
                }

                // Calculate min and max scalar fluxes (not used elsewhere
                // currently)
                if (oct == sn_vars->noct)
                {
                    dim_sweep_vars->fmin = (((dim_sweep_vars->fmin) < (solvar_vars->flux[ (g-1) * input_vars->nz * input_vars->ny * input_vars->nx + (k-1) * input_vars->ny * input_vars->nx + (j-1) * input_vars->nx + (i-1) ])) ? (dim_sweep_vars->fmin) : (solvar_vars->flux[ (g-1) * input_vars->nz * input_vars->ny * input_vars->nx + (k-1) * input_vars->ny * input_vars->nx + (j-1) * input_vars->nx + (i-1) ]));
                    dim_sweep_vars->fmax = (((dim_sweep_vars->fmax) > (solvar_vars->flux[ (g-1) * input_vars->nz * input_vars->ny * input_vars->nx + (k-1) * input_vars->ny * input_vars->nx + (j-1) * input_vars->nx + (i-1) ])) ? (dim_sweep_vars->fmax) : (solvar_vars->flux[ (g-1) * input_vars->nz * input_vars->ny * input_vars->nx + (k-1) * input_vars->ny * input_vars->nx + (j-1) * input_vars->nx + (i-1) ]));
                }
                
                // Save edge fluxes (dummy if checks for unused non-vacuum BCs)
                if (j == j_high)
                {
                    if ((j_dir==2 && para_vars->lasty) ||
                        ((j_dir == 1 && para_vars->firsty) && ibb == 1))
                    {
                        
                    }
                    else
                    {
                        for (ang = 0; ang < input_vars->nang; ang++)
                        {
                            solvar_vars->jb_out[ (g-1) * input_vars->nz * input_vars->ichunk * input_vars->nang + (k-1) * input_vars->ichunk * input_vars->nang + (ic-1) * input_vars->nang + ang ] = solvar_vars->psij[ (g-1) * input_vars->nz * input_vars->ichunk * input_vars->nang + (k-1) * input_vars->ichunk * input_vars->nang + (ic-1) * input_vars->nang + ang ];
                        }
                    }
                }
                
                if (k == k_high)
                {
                    if ((k_dir == 2 && para_vars->lastz) ||
                        ((k_dir==1 && para_vars->firstz) && ibf == 1))
                    {
                        
                    }
                    else
                    {
                        for ( ang = 0; ang < input_vars->nang; ang++ )
                        {
                            solvar_vars->kb_out[ (g-1) * input_vars->ny * input_vars->ichunk * input_vars->nang + (j-1) * input_vars->ichunk * input_vars->nang + (ic-1) * input_vars->nang + ang ] = solvar_vars->psik[ (g-1) * input_vars->ny * input_vars->ichunk * input_vars->nang + (j-1) * input_vars->ichunk * input_vars->nang + (ic-1) * input_vars->nang + ang ];
                        }
                    }
                }

                // Compute leakages (not used elsewhere currently)
                if (((i+i_dir-1) == 1) || ((i+i_dir-1) == (input_vars->nx+1)))
                {
                    sum_wmupsii = 0;

                    for (ang = 0; ang < input_vars->nang; ang++)
                    {
                        sum_wmupsii += sn_vars->wmu[ang] * solvar_vars->psii[ (g-1) * input_vars->nz * input_vars->ny * input_vars->nang + (k-1) * input_vars->ny * input_vars->nang + (j-1) * input_vars->nang + ang ];
                    }

                    solvar_vars->flkx[ (g-1) * input_vars->nz * input_vars->ny * (input_vars->nx+1) + (k-1) * input_vars->ny * (input_vars->nx+1) + (j-1) * (input_vars->nx+1) + (i+i_dir-1-1) ] += i_step*sum_wmupsii;
                }
                if ((j_dir == 1 && para_vars->firsty) || (j_dir == 2 && para_vars->lasty))
                {
                    sum_wetapsij = 0;

                    for (ang = 0; ang < input_vars->nang; ang++)
                    {
                        sum_wetapsij
                            += sn_vars->weta[ang] * solvar_vars->psij[ (g-1) * input_vars->nz * input_vars->ichunk * input_vars->nang + (k-1) * input_vars->ichunk * input_vars->nang + (ic-1) * input_vars->nang + ang ];
                    }

                    solvar_vars->flky[ (g-1) * input_vars->nz * (input_vars->ny+1) * input_vars->nx + (k-1) * (input_vars->ny+1) * input_vars->nx + (j+j_dir-1-1) * input_vars->nx + (i-1) ] += j_step*sum_wetapsij;
                }

                if (((k_dir == 1 && para_vars->firstz) || (k_dir == 2 && para_vars->lastz)) && input_vars->ndimen == 3)
                {
                    sum_wxipsik = 0;

                    for (ang = 0; ang < input_vars->nang; ang++)
                    {
                        sum_wxipsik += sn_vars->wxi[ang] * solvar_vars->psik[ (g-1) * input_vars->ny * input_vars->ichunk * input_vars->nang + (j-1) * input_vars->ichunk * input_vars->nang + (ic-1) * input_vars->nang + ang ];
                    }

                    solvar_vars->flkz[ (g-1) * (input_vars->nz+1) * input_vars->ny * input_vars->nx + (k+k_dir-1-1) * input_vars->ny * input_vars->nx + (j-1) * input_vars->nx + (i-1) ]
                        += k_step*sum_wxipsik;
                }
            }
        } 
    } 
}