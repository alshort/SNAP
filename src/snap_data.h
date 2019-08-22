// input.c
/***********************************************************************
 * Input Variables
 *
 * Parallel processing inputs:
 * npey       - number of PEs in y-x plain
 * npez       - number of PEs in z-x plain
 * ichunk     -
 * nthreads   - number of OpenMP threads
 * nnested    - number of nested threads
 *
 * Geometry inputs:
 * ndimen     - number of spatial dimensions 1/2/3
 * nx         - number of x-dir spatial cells (global)
 * ny         - number of y-dir spatial cells
 *              (global on input, reset to per PE in setup)
 * nz         - number of z-dir spatial cells
 *              (global on input, reset to per PE in setup)
 * lx         - total length of x domain
 * ly         - total length of y domain
 * lz         - total length of z domain
 *
 * Sn inputs:
 * nmom       - number of discrete ordinates per octant
 * nang       - scattering order
 *
 * Data inputs:
 * ng         - number of groups
 * mat_opt    - material layout, 0/1/2=homogeneous/center/corner, with
 *              two materials, and material 2 nowhere/center/corner
 * src_opt    - source layout, 0/1/2=homogenous/src-center/src-corner=
 *              source everywhere/center of problem/corner, strength=10.0
 * scatp      - 0/1=no/yes print the full scattering matrix to file 'slgg'
 *
 * Control inputs:
 * epsi       -
 * tf         -
 * iitm       -
 * oitm       -
 * timedep    -
 * nsteps     -
 * it_det     -
 * fluxp      -
 * fixup      -
 ***********************************************************************/
typedef struct input_data_
{
    // Parallel processing inputs
    int npey, npez, ichunk, nthreads, nnested;

    // Geometry inputs
    int ndimen, nx, ny, nz;
    double lx, ly, lz;

    // Sn inputs
    int nmom, nang;

    // Data inputs
    int ng, mat_opt, src_opt, scatp;

    // Control inputs
    double epsi, tf;
    int iitm, oitm, timedep, nsteps, it_det, fluxp, fixup;

} input_data;

// time.c
/***********************************************************************
 * Time run-time variables
 *
 * tsnap    - total SNAP run time
 * tparset  - parallel environment setup time
 * tinp     - input run time
 * tset     - setup run time
 * tslv     - total solution run time
 * tparam   - time for setting up solve parameters
 * totrsrc  - time for outer source computations
 * tinners  - total time spent on inner iterations
 * tinrsrc  - time for inner source computations
 * tsweeps  - time for transport sweeps, including angular sourc compute
 * tinrmisc - time for miscellaneous inner ops
 * tslvmisc - time for miscellaneous solution ops
 * tout     - output run time
 * tgrind   - transport grind time
 ***********************************************************************/
typedef struct time_data_
{
    double tsnap, tparset, tinp, tset, tslv, tparam, totrsrc, tinners,
        tinrsrc, tsweeps, tinrmisc, tslvmisc, tout, tgrind;

} time_data;

// geom.c
/***********************************************************************
 * Geometry run-time variables
 *
 * ny_gl    - global number of y-dir spatial cells
 * nz_gl    - global number of z-dir spatial cells
 * jlb      - global index of local lower y bound
 * jub      - global index of local upper y bound
 * klb      - global index of local lower z bound
 * kub      - global index of local upper z bound
 *
 * dx       - x width of spatial cell
 * dy       - y width of spatial cell
 * dz       - z width of spatial cell
 *
 * nc       - number of i-chunks, nx/ichunk
 *
 * hi       - Spatial DD x-coefficient
 * hj(nang) - Spatial DD y-coefficient
 * hk(nang) - Spatial DD z-coefficient
 *
 * dinv(nang,nx,ny,nz,ng) - Sweep denominator, pre-computed/inverted
 *
 * ndiag    - number of diagonals of mini-KBA sweeps in nested threading
 ***********************************************************************/
typedef struct cell_id_type_
{
    int ic, jc, kc;

} cell_id_type;

typedef struct diag_type_
{
    int lenc;

    cell_id_type *cell_id_vars;

} diag_type;

typedef struct geom_data_
{
    int ny_gl, nz_gl, jlb, jub, klb, kub, nc, ndiag;

    double dx, dy, dz, hi;

    double *hj, *hk; // 1D arrays

    double *dinv; // 5D array

    diag_type *diag_vars; // 1D array

} geom_data;


// sn.c
/***********************************************************************
 * SN run-time variables
 *
 * cmom       - computational number of moments according to nmom & ndimen
 * noct       - number of directional octants
 * mu(nang)   - x direction cosines
 * eta(nang)  - y direction cosines
 * xi(nang)   - z direction cosines
 * w(nang)    - angle weights
 *
 * wmu(nang)  - w*mu
 * weta(nang) - w*eta
 * wxi(nang)  - w*xi
 *
 * ec(nang,cmom,noct) - Scattering expansion coefficients
 * lma(nmom)          - number of 'm' moments per order l
 ***********************************************************************/
typedef struct sn_data_
{
    int cmom, noct;

    double *mu, *eta, *xi, *w, *wmu,*weta, *wxi; // 1-D arrays

    double *ec; // 3-D array

    int *lma;

} sn_data;

// data.c
/***********************************************************************
 * Data run-time variables
 *
 * v(ng)         - mock velocity array
 * nmat          - number of materials
 * mat(nx,ny,nz) - material identifier array
 *
 * qi(nx,ny,nz,ng)             - fixed source array for src_opt<3
 * qim(nang,nx,ny,nz,noct,ng)  - fixed source array for src_opt>=3

 * sigt(nmat,ng)          - total interaction
 * siga(nmat,ng)          - absorption
 * sigs(nmat,ng)          - scattering, total
 * slgg(nmat,nmom,ng,ng)  - scattering matrix, all moments/groups
 * vdelt(ng)              - time-absorption coefficient
 ***********************************************************************/
typedef struct data_data_
{
    int nmat;
    int *mat;  // 3-D array

    double *v, *vdelt;          // 1-D arrays
    double *sigt, *siga, *sigs; // 2-D arrays
    double *qi, *slgg;          // 4-D arrays
    double *qim;                // 6-D arrays

} data_data;

// conrol.c
/***********************************************************************
 * control run-time variables
 *
 * dt       - time-step size
 *
 * tolr      - parameter, small number used for determining how to
 *             compute flux error
 * dfmxi(ng) - max error of inner iteration
 * dfmxo     - max error of outer iteration
 *
 * inrdone(ng)  - logical for inners being complete
 * otrdone      - logical for outers being complete
 ***********************************************************************/
typedef struct control_data_
{
    bool otrdone;
    bool *inrdone;     // 1-D array

    double dt, dfmxo, tolr;

    double *dfmxi;     // 1-D array

} control_data;

// mms.c
/***********************************************************************
 * mms variables
 *
 * ref_flux(nx,ny,nz,ng)          - Manufactured solution
 * ref_fluxm(cmom-1,nx,ny,nz,ng)  - Manufactured solution moments
 *
 * a_const       - i function constant
 * b_const       - j function constant
 * c_const       - k function constant
 *
 * ib(nx+1)      - i cell boundaries
 * jb(ny+1)      - j cell boundaries
 * kb(nz+1)      - k cell boundaries
 ***********************************************************************/
typedef struct mms_data_
{
    double *ref_flux;   // 4-D array
    double *ref_fluxm;  // 5-D array

    double a_const, b_const, c_const;

    double *ib, *jb, *kb; // 1-D arrays

} mms_data;

// solvar.c
/***********************************************************************
 * Solvar Module variables
 *
 * ptr_in(nang,nx,ny,nz,noct,ng)   - Incoming time-edge flux pointer
 * ptr_out(nang,nx,ny,nz,noct,ng)  - Outgoing time-edge flux pointer
 *
 * flux(nx,ny,nz,ng)          - Scalar flux moments array
 * fluxpo(nx,ny,nz,ng)        - Previous outer copy of scalar flux array
 * fluxpi(nx,ny,nz,ng)        - Previous inner copy of scalar flux array
 * fluxm(cmom-1,nx,ny,nz,ng)  - Flux moments array
 *
 * q2grp(cmom,nx,ny,nz,ng)  - Out-of-group scattering + fixed sources
 * qtot(cmom,nx,ny,nz,ng)   - Total source: q2grp + within-group source
 *
 * t_xs(nx,ny,nz,ng)       - Total cross section on mesh
 * a_xs(nx,ny,nz,ng)       - Absorption cross section on mesh
 * s_xs(nmom,nx,ny,nz,ng)  - In-group scattering cross section on mesh
 *
 * psii(nang,ny,nz,ng)     - Working psi_x array
 * psij(nang,ichunk,nz,ng) - Working psi_y array
 * psik(nang,ichunk,ny,ng) - Working psi_z array
 *
 * jb_in(nang,ichunk,nz,ng)  - y-dir boundary flux in from comm
 * jb_out(nang,ichunk,nz,ng) - y-dir boundary flux out to comm
 * kb_in(nang,ichunk,ny,ng)  - z-dir boundary flux in from comm
 * kb_out(nang,ichunk,ny,ng) - z-dir boundary flux out to comm
 *
 * flkx(nx+1,ny,nz,ng)     - x-dir leakage array
 * flky(nx,ny+1,nz,ng)     - y-dir leakage array
 * flkz(nx,ny,nz+1,ng)     - z-dir leakage array
 *
 ***********************************************************************/
typedef struct solvar_data_
{
    double *flux, *fluxpo, *fluxpi, *t_xs, *a_xs, // 4-D arrays
        *psii, *psij, *psik, *jb_in, *jb_out,
        *kb_in, *kb_out, *flkx, *flky, *flkz;

    double *qtot, *q2grp, *fluxm, *s_xs;          // 5-D arrays

    double *ptr_in, *ptr_out;                     // 6-D arrays
} solvar_data;


// dim1_sweep.c and dim3_sweep.c
/***********************************************************************
 * dim_sweep module variable
 *
 * fmin        - min scalar flux. Dummy for now, not used elsewhere.
 * fmax        - max scalar flux. Dummy for now, not used elsewhere.
 ***********************************************************************/
typedef struct dim_sweep_data_
{
    double fmin, fmax;
} dim_sweep_data;


// sweep.c
/***********************************************************************
 * Module variables
 *
 * mtag               - message tag
 * <y or z>p_snd      - rank of process sending to
 * <y or z>p_rcv      - rank of process receiving from
 * incoming<y or z>   - logical determining if proc will receive a message
 * outgoing<y or z>   - logical determining if proc will send a message
 ***********************************************************************/
typedef struct sweep_data_
{
    int mtag, yp_snd, yp_rcv, zp_snd, zp_rcv;

    bool incomingy, incomingz, outgoingy, outgoingz;

} sweep_data;