/***********************************************************************
 * Header file for c version of SNAP.
 ***********************************************************************/
#ifndef _SNAP_H
#define _SNAP_H

#include <stdlib.h>
#include <stdio.h>
#include <malloc.h>
#include <mpi.h>
#include <omp.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>

#include "snap_data.h"

#ifdef USEMKL
#include "mkl.h"
#define INTEL_BB 64
#define VML_ACCURACY VML_EP
#define VML_HANDLING VML_FTZDAZ_OFF
#define VML_ERROR VML_ERRMODE_DEFAULT
#define VECLEN_MIN 40
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/***********************************************************************
 * Typedef functions
 ***********************************************************************/
// plib.c
/***********************************************************************
 * Parallel Run-time Varibales
 *
 * Note all ranks are zero based
 *
 * root       - root process for comm_snap, 0
 *
 * nproc      - number of MPI processes
 * iproc      - rank of calling process in base communicator
 *
 * comm_snap  - base communicator, duplicated from MPI_COMM_WORLD
 * comm_space - SDD communicator, ndimen-1 grid for 2-D (x-y) or
 *              3-D (x-y-z) problems. Non-existent for 1-D (x) problems.
 * sproc      - rank of calling process in comm_space
 *
 * ycomm      - y-dimension process communicator
 * zcomm      - z-dimension process communicator
 * yproc      - PE column in SDD 2-D PE mesh (comm_space)
 * zproc      - PE row in SDD 2-D PE mesh (comm_space)
 * firsty     - logical determining if lowest yproc
 * lasty      - logical determining if highest yproc
 * firstz     - logical determining if lowest zproc
 * lastz      - logical determining if highest zproc
 * ylop       - rank of preceding yproc in ycomm
 * yhip       - rank of succeeding yproc in ycomm
 * zlop       - rank of preceding zproc in zcomm
 * zhip       - rank of succeeding zproc in zcomm
 *
 * g_off      - group offset for message tags
 *
 * thread_level       - level of MPI thread support
 * thread_single      - MPI_THREAD_SINGLE
 * thread_funneled    - MPI_THREAD_FUNNELED
 * thread_serialized  - MPI_THREAD_SERIALIZED
 * thread_multiple    - MPI_THREAD_MULTIPLE
 * lock               - OpenMP lock
 *
 * num_grth   - minimum number of nthreads and ng; used to ensure loop
 *              over groups with communications is sized properly
 * do_nested  - true/false use nested threading, i.e., mini-KBA
 ***********************************************************************/
typedef struct para_data_
{
    int root, g_off;

    MPI_Comm comm_snap, comm_space;

    int nproc, iproc, sproc, ycomm,
        zcomm, yproc, zproc, ylop, yhip, zlop, zhip, thread_level,
        thread_single, thread_funneled, thread_serialized, thread_multiple,
        max_threads, num_grth;

    omp_lock_t lock;

    bool firsty, lasty, firstz, lastz, do_nested;

} para_data;



/***********************************************************************
 * Function prototypes
 ***********************************************************************/
/* plib.c: Contains the variables that control parallel decomposition and the
 subroutines for parallel enviornment setup. Only module that requires MPI
 library interaction except for time.c (MPI_WTIME). */
// Constructor for para_data type
void para_data_init ( para_data *para_vars );

// Initialize the MPI process environment
void pinit ( int argc, char *argv[], para_data *para_vars,
             double *time, int *ierr );

// Call to execute an MPI_Barrier
int barrier ( MPI_Comm comm );

// Setup the SDD communicator
void pcomm_set( int npey, int npez,  para_data *para_vars, int *ierr );

// Call to execute MPI_Finalize
void pend ( void );

// All reduce global max value (integer)
int glmax_i ( int *value, MPI_Comm comm );

// All reduce global max value (double)
int glmax_d ( double *value, MPI_Comm comm );

// All reduce global max value (double) for 1-d array
int glmax_d_1d ( double *value, int dlen, MPI_Comm comm );

// All reduce global min value (integer)
int glmin_i ( int *value, MPI_Comm comm );

// All reduce global min value (double)
int glmin_d ( double *value, MPI_Comm comm );

// Broadcast (integer) scalar
int bcast_i_scalar ( int *value, MPI_Comm comm, int bproc, int nproc );

// Broadcast (double) scalar
int bcast_d_scalar ( double *value, MPI_Comm comm,
                     int bproc, int nproc );

// Broadcast (integer) 1-d array
int bcast_i_1d ( int *value, int ilen, MPI_Comm comm,
                 int bproc, int nproc );

// Broadcast (double) 1-d array
int bcast_d_1d ( double *value, int dlen, MPI_Comm comm,
                 int bproc, int nproc );

// Send a rank-2 double presicision array
int psend_d_2d ( double *value, int d1, int d2, MPI_Comm comm,
                 int proc, int myproc, int mtag);

// Send a rank-3 double presicision array
int psend_d_3d ( double *value, int d1, int d2, int d3, MPI_Comm comm,
                 int proc, int myproc, int mtag);

// Receive a rank-2 double presicision array
int precv_d_2d ( double *value, int d1, int d2, MPI_Comm comm,
                 int proc, int myproc, int mtag);

// Receive a rank-3 double presicision array
int precv_d_3d ( double *value, int d1, int d2, int d3, MPI_Comm comm,
                 int proc, int myproc, int mtag);

// Return rank of proc defined by coordinates of Cartesian communicator
int cartrank ( int *coord, int *rank, MPI_Comm comm );

// Setup the number of OpenMP threads. Check if any proc is exceeding
// max threads. Reset and report if so.
void pinit_omp( MPI_Comm comm, int *nthreads, int nnested,
                bool *do_nested, int *ierr, char **error );

// Operate on an OpenMP lock
void plock_omp ( char *dowhat, omp_lock_t *lock );

// Return thread number of caller
int thread_num( void );


// utils.c
int cmdarg ( int argc, char *argv[], char **inputFile, char **outputFile,
             char **error, int iproc, int root );

int open_file ( FILE **fp, char *fileName, char *fileAction,
                char **error, int iproc, int root );

int close_file ( FILE *fp, char *fileName, char **error,
                 int iproc, int root );

void print_error ( FILE *fp, char *error, int iproc, int root );

int string_empty ( char *stringName );

void stop_run ( int inputFlag, int solveFlag, int statusFlag, para_data *para_vars,
                sn_data *sn_vars, data_data *data_vars, mms_data *mms_vars,
                geom_data *geom_vars, solvar_data *solvar_vars, control_data *control_vars );


// dealloc.c
void dealloc_input ( int selectFlag, sn_data *sn_vars,
                     data_data *data_vars, mms_data *mms_vars );

void dealloc_solve ( int selectFlag, geom_data *geom_vars,
                     solvar_data *solvar_vars, control_data *control_vars );


// version.c
void version_print ( FILE *fp_out );


// input.c
void input_data_init ( input_data *input_vars );

int read_input ( FILE *fp_in, FILE *fp_out, input_data *input_vars,
                 para_data *para_vars, time_data *time_vars );

void get_input_value ( char *lineData, char *valueID, char **tmpData );

void input_echo ( input_data *input_vars, FILE *fp_out );

int input_check ( FILE *fp_out, input_data *input_vars, para_data *para_vars );

int var_bcast ( input_data *input_vars, para_data *para_vars );


// time.c
void time_data_init ( time_data *time_vars );

double wtime ( void );

void time_summ ( FILE *fp_out, time_data *time_vars );


//setup.c
void setup ( input_data *input_vars, para_data *para_vars, time_data *time_vars,
             geom_data *geom_vars, sn_data *sn_vars, data_data *data_vars,
             solvar_data *solvar_vars, control_data *control_vars,
             mms_data *mms_vars, FILE *fp_out, int *ierr, char **error );

void setup_alloc( input_data *input_vars, para_data *para_vars, sn_data *sn_vars,
                  data_data *data_vars, int *flg, int *ierr, char **error );

void setup_delta( input_data *input_vars, geom_data *geom_vars,
                  control_data *control_vars );

void setup_vel ( input_data *input_vars, data_data *data_vars );

void setup_angle ( input_data *input_vars, sn_data *sn_vars);

void setup_mat ( input_data *input_vars, geom_data *geom_vars, data_data *data_vars,
                 int *i1, int *i2, int *j1, int *j2, int *k1, int *k2 );

void setup_data ( input_data *input_vars, data_data *data_vars);

void setup_src ( input_data *input_vars, para_data *para_vars, geom_data *geom_vars,
                 sn_data *sn_vars, data_data *data_vars, control_data *control_vars,
                 mms_data *mms_vars, int *i1, int *i2, int *j1, int *j2, int *k1,
                 int *k2, int *ierr, char **error );

void setup_echo ( FILE *fp_out, input_data *input_vars, para_data *para_vars,
                  geom_data *geom_vars, data_data *data_vars, sn_data *sn_vars,
                  control_data *control_vars, int mis, int mie, int mjs, int mje, int mks,
                  int mke, int qis, int qie, int qjs, int qje, int qks, int qke );

void setup_scatp( input_data *input_vars, para_data *para_vars,
                  data_data *data_vars, int *ierr, char **error );


// geom.c
void geom_data_init ( geom_data *geom_vars );

void geom_alloc ( input_data *input_vars, geom_data *geom_vars, int *ierr );

void geom_dealloc ( geom_data *geom_vars );

void param_calc ( input_data *input_vars, sn_data *sn_vars,
                  solvar_data *solvar_vars, data_data *data_vars,
                  geom_data *geom_vars, int ng_indx );

void diag_setup ( input_data *input_vars, para_data *para_vars,
                  geom_data *geom_vars, int *ierr, char **error );


// sn.c
void sn_data_init ( sn_data *sn_vars );

void sn_allocate ( sn_data *sn_vars, input_data *input_vars, int *ierr );

void sn_deallocate ( sn_data *sn_vars );

void expcoeff ( input_data *input_vars, sn_data *sn_vars, int *ndimen );


/* data.c: contains the variables and setup subroutines for the mock
 cross section data. It establishes the number of groups and constructs
 the cross section arrays.*/
// Constructor for data_data type
void data_data_init (data_data *data_vars );

// Allocate data module arrays
void data_allocate ( data_data *data_vars, input_data *input_vars,
                     sn_data *sn_vars, int *ierr );

// Deallocate the data module arrays
void data_deallocate ( data_data *data_vars );


// control.c
void control_data_init ( control_data *control_vars );

void control_alloc ( input_data *input_vars, control_data *control_vars, int *ierr);

void control_dealloc ( control_data *control_vars );


// mms.c
void mms_data_init ( mms_data *mms_vars );

void mms_setup ( input_data *input_vars, para_data *para_vars, geom_data *geom_vars,
                 data_data *data_vars, sn_data *sn_vars, control_data *control_vars,
                 mms_data *mms_vars, int *ierr, char **error );

void mms_allocate ( input_data *input_vars, sn_data *sn_vars, mms_data *mms_vars,
                    int *ierr, char **error );

void mms_deallocate (  mms_data *mms_vars );

void mms_cells ( input_data *input_vars, geom_data *geom_vars, mms_data *mms_vars );

void mms_flux_1 ( input_data *input_vars, geom_data *geom_vars,
                  sn_data *sn_vars, mms_data *mms_vars );

void mms_trigint ( char *trig, int lc, double d, double del,
                   double *cb, double *fn );

void mms_src_1( input_data *input_vars, geom_data *geom_vars, data_data *data_vars,
                sn_data *sn_vars, mms_data *mms_vars );

void mms_flux_1_2 ( input_data *input_vars, control_data *control_vars,
                    mms_data *mms_vars );

void mms_verify_1 ( input_data *input_vars, para_data *para_vars, control_data *control_vars,
                    mms_data *mms_vars, solvar_data *solvar_vars, FILE *fp_out );


// translv.c
void translv ( input_data *input_vars, para_data *para_vars, time_data *time_vars,
               geom_data *geom_vars, sn_data *sn_vars, data_data *data_vars,
               control_data *control_vars, solvar_data *solvar_vars, mms_data *mms_vars,
               sweep_data *sweep_vars, dim_sweep_data *dim_sweep_vars,
               FILE *fp_out, int *ierr, char **error );

// solvar.c
void solvar_data_init ( solvar_data *solvar_vars );

void solvar_alloc ( input_data *input_vars, sn_data* sn_vars, solvar_data *solvar_vars,
                    int *ierr );

void solvar_dealloc ( solvar_data *solvar_vars );


// dim1_sweep.c
void dim1_sweep_data_init ( dim_sweep_data *dim_sweep_vars );

void dim1_sweep ( input_data *input_vars, geom_data *geom_vars, sn_data *sn_vars, data_data *data_vars,
                  control_data *control_vars, solvar_data *solvar_vars, dim_sweep_data *dim_sweep_vars,
                  int id, int oct, int d1, int d2, int d3, int d4, int i1, int i2, int g, int *ierr );

// dim3_sweep.c
// void dim3_sweep_data_init ( dim_sweep_data *dim_sweep_vars );

// void dim3_sweep ( input_data *input_vars, para_data *para_vars,
//                   geom_data *geom_vars, sn_data *sn_vars,
//                   data_data *data_vars, control_data *control_vars,
//                   solvar_data *solvar_vars, dim_sweep_data *dim_sweep_vars,
//                   int ich, int id, int d1, int d2, int d3, int d4, int jd,
//                   int kd, int jlo, int klo, int jhi, int khi, int jst, int kst,
//                   int i1, int i2, int oct, int g, int *ierr );


// sweep.c
void sweep_data_init (sweep_data *sweep_vars );

void sweep ( input_data *input_vars, para_data *para_vars, geom_data *geom_vars,
             sn_data *sn_vars, data_data *data_vars, control_data *control_vars,
             solvar_data *solvar_vars, sweep_data *sweep_vars,
             dim_sweep_data *dim_sweep_vars, int *ierr );
/*
void sweep_recv_bdry ( para_data *para_vars, sweep_data *sweep_vars,solvar_data *solvar_vars,
                       input_data *input_vars,
                       double *value, int d1, int d2, int d3, MPI_Comm comm,
                       int proc, int myproc, int g, int iop );

void sweep_send_bdry ( para_data *para_vars, sweep_data *sweep_vars, solvar_data *solvar_vars,
                       input_data *input_vars,
                       double *value, int d1, int d2, int d3, MPI_Comm comm,
                       int proc, int myproc, int g, int iop );
*/
void sweep_recv_bdry ( input_data *input_vars, para_data *para_vars,
                       solvar_data *solvar_vars, sweep_data *sweep_vars, int g, int iop );

void sweep_send_bdry ( input_data *input_vars, para_data *para_vars,
                       solvar_data *solvar_vars, sweep_data *sweep_vars, int g, int iop );


// octsweep.c
void octsweep ( input_data *input_vars, para_data *para_vars, geom_data *geom_vars, sn_data *sn_vars,
                data_data *data_vars, control_data *control_vars,
                solvar_data *solvar_vars, dim_sweep_data *dim_sweep_vars,
                int g, int iop, int jd, int kd, int jlo, int jhi, int jst,
                int klo, int khi, int kst, int *ierr );


// expxs.c
#ifdef USEMACRO
void expxs_reg ( double *xs, int *map, double *cs, int nmat_size, int nx_size,
                 int ny_size, int nz_size, int ng_indx, int ng_size );

void expxs_slgg ( double *scat, int *map, double *cs, int nmat_size, int nmom_size,
                  int nx_size, int ny_size, int nz_size, int ng_indx, int ng_size );
#else
void expxs_reg ( input_data *input_vars, data_data *data_vars, solvar_data *solvar_vars,
                 double *cs, int ng_indx, int gp_indx, int l_indx );

void expxs_slgg ( input_data *input_vars, data_data *data_vars,
                  solvar_data *solvar_vars, int ng_indx );
#endif


// outer.c
int outer ( input_data *input_vars, para_data *para_vars, geom_data *geom_vars,
            time_data *time_vars, data_data *data_vars, sn_data *sn_vars,
            control_data *control_vars, solvar_data *solvar_vars,
            sweep_data *sweep_vars, dim_sweep_data *dim_sweep_vars,
            FILE *fp_out, int *ierr );

void otr_src ( input_data *input_vars, data_data *data_vars,
               sn_data *sn_vars, solvar_data *solvar_vars, int *ierr );

void otr_src_scat ( input_data *input_vars, data_data *data_vars, sn_data *sn_vars,
                    solvar_data *solvar_vars, int g, int gp, int *ierr );

void otr_conv ( input_data *input_vars, para_data *para_vars, control_data *control_vars,
                solvar_data *solvar_vars, int *ierr );


// inner.c
void inner ( input_data *input_vars, para_data *para_vars, geom_data *geom_vars,
             time_data *time_vars, data_data *data_vars, sn_data *sn_vars,
             control_data *control_vars, solvar_data *solvar_vars,
             sweep_data *sweep_vars, dim_sweep_data *dim_sweep_vars, int inno,
             int *iits, FILE *fp_out, int *ierr );

void inr_src ( input_data *input_vars, sn_data *sn_vars,
               control_data *control_vars, solvar_data *solvar_vars, para_data *para_vars );

void inr_src_scat ( input_data *input_vars, sn_data *sn_vars,
                    solvar_data *solvar_vars, int g, para_data *para_vars );

void inr_conv ( input_data *input_vars, para_data *para_vars,
                control_data *control_vars, solvar_data *solvar_vars,
                int inno, int *iits, FILE *fp_out, int *ierr );


// output.c
void output ( input_data *input_vars, para_data *para_vars, time_data *time_vars,
              geom_data *geom_vars, data_data *data_vars, sn_data *sn_vars,
              control_data *control_vars, mms_data *mms_vars, solvar_data *solvar_vars,
              sweep_data *sweep_vars, FILE *fp_out, int *ierr, char **error );

//void output_send ( input_data *input_vars, para_data *para_vars,
//                   control_data *control_vars, sweep_data *sweep_vars,
//                   double *fprnt, int *ierr );

//void output_recv ( input_data *input_vars, para_data *para_vars,
//                   control_data *control_vars, sweep_data *sweep_vars,
//                   double *fprnt, int *ierr );

void output_send ( int dim1, int dim2, MPI_Comm comm, int root, int sproc,
                   int mtag, double *fprnt, int *ierr );

void output_recv ( int dim1, int dim2, MPI_Comm comm, int proc, int sproc,
                   int mtag, double *fprnt, int *ierr );


void output_flux_file ( input_data *input_vars, para_data *para_vars,
                        geom_data *geom_vars, data_data *data_vars,
                        sn_data *sn_vars, control_data *control_vars,
                        mms_data *mms_vars, solvar_data *solvar_vars,
                        sweep_data *sweep_vars, int klb, int kub,
                        int *ierr, char **error, FILE *fp_out );


/***********************************************************************
 * Memory allocation macros
 ***********************************************************************/
#ifdef USEMKL

#define ALLOC_STR(PNTR, STRLEN, IERR)                           \
    if ( !PNTR )                                                \
    {                                                           \
        PNTR = (char *) mkl_malloc ( STRLEN, INTEL_BB );        \
        if (!PNTR )                                             \
        {                                                       \
            perror("ALLOC_STR");                                \
            fprintf(stderr,                                     \
                    "Allocation failed for " #PNTR              \
                    ".  Terminating...\n");                     \
            *IERR = 1; /* exit(-1); */                         \
        }                                                       \
    }                                                           \
    else                                                        \
    {                                                           \
        PNTR = (char *) mkl_realloc ( PNTR, STRLEN );           \
        if (!PNTR )                                             \
        {                                                       \
            perror("ALLOC_STR");                                \
            fprintf(stderr,                                     \
                    "Re-allocation failed for " #PNTR           \
                    ".  Terminating...\n");                     \
            *IERR = 1; /* exit(-1); */                         \
        }                                                       \
    }

#define REALLOC_STR(PNTR, STRLEN, IERR)                 \
    PNTR = (char *) mkl_realloc ( PNTR, STRLEN );       \
    if (!PNTR)                                          \
    {                                                   \
        perror("REALLOC_STR");                          \
        fprintf(stderr,                                 \
                "Re-allocation failed for " #PNTR       \
                ".  Terminating...\n");                 \
        *IERR = 1; /* exit(-1); */                      \
    }

#define FREE(PNTR)                                     \
    if ( PNTR )                                        \
    {                                                  \
        mkl_free ( PNTR );                             \
    }

#define REALLOC_2D(PNTR, NUMX, NUMY, TYPE, IERR)                        \
    PNTR = (TYPE *)mkl_realloc(PNTR, NUMX*NUMY*sizeof(TYPE));                 \
    if (!PNTR)                                                          \
    {                                                                   \
        perror("REALLOC_2D");                                             \
        fprintf(stderr,                                                 \
                "Reallocation failed for " #PNTR ".  Terminating...\n");  \
        *IERR = 1; /* exit(-1); */                                     \
    }


#define ALLOC_1D(PNTR, NUM, TYPE, IERR)                                 \
    PNTR = (TYPE *)mkl_calloc(NUM, sizeof(TYPE), INTEL_BB);             \
    if (!PNTR)                                                          \
    {                                                                   \
        perror("ALLOC_1D");                                             \
        fprintf(stderr,                                                 \
                "Allocation failed for " #PNTR ".  Terminating...\n");  \
        *IERR = 1; /* exit(-1); */                                     \
    }

#define ALLOC_2D(PNTR, NUMX, NUMY, TYPE, IERR)                          \
    PNTR = (TYPE *)mkl_calloc((NUMX) * (NUMY),                          \
                              sizeof(TYPE), INTEL_BB);                  \
    if (!PNTR)                                                          \
    {                                                                   \
        perror("ALLOC_2D");                                             \
        fprintf(stderr,                                                 \
                "Allocation failed for " #PNTR ".  Terminating...\n");  \
        *IERR = 1; /* exit(-1); */                                     \
    }

#define ALLOC_3D(PNTR, NUMX, NUMY, NUMZ, TYPE, IERR)                    \
    PNTR = (TYPE *)mkl_calloc((NUMX) * (NUMY) * (NUMZ),                 \
                              sizeof(TYPE), INTEL_BB);                  \
    if (!PNTR)                                                          \
    {                                                                   \
        perror("ALLOC_3D");                                             \
        fprintf(stderr,                                                 \
                "Allocation failed for " #PNTR ".  Terminating...\n");  \
        *IERR = 1; /* exit(-1); */                                     \
    }

#define ALLOC_4D(PNTR, NUMW, NUMX, NUMY, NUMZ, TYPE, IERR)              \
    PNTR = (TYPE *)mkl_calloc((NUMW) * (NUMX) * (NUMY) * (NUMZ),        \
                              sizeof(TYPE), INTEL_BB);                  \
    if (!PNTR)                                                          \
    {                                                                   \
        perror("ALLOC_4D");                                             \
        fprintf(stderr,                                                 \
                "Allocation failed for " #PNTR ".  Terminating...\n");  \
        *IERR = 1; /* exit(-1); */                                     \
    }

#define ALLOC_5D(PNTR, NUMV, NUMW, NUMX, NUMY, NUMZ, TYPE, IERR)        \
    PNTR = (TYPE *)mkl_calloc((NUMV) * (NUMW)                           \
                              * (NUMX) * (NUMY) * (NUMZ),               \
                              sizeof(TYPE), INTEL_BB);                  \
    if (!PNTR)                                                          \
    {                                                                   \
        perror("ALLOC_5D");                                             \
        fprintf(stderr,                                                 \
                "Allocation failed for " #PNTR ".  Terminating...\n");  \
        *IERR = 1; /* exit(-1); */                                     \
    }

#define ALLOC_6D(PNTR, NUMU, NUMV, NUMW, NUMX, NUMY, NUMZ, TYPE, IERR)  \
    PNTR = (TYPE *)mkl_calloc((NUMU) * (NUMV) *(NUMW)                   \
                              * (NUMX) * (NUMY) * (NUMZ),               \
                              sizeof(TYPE), INTEL_BB);                  \
    if (!PNTR)                                                          \
    {                                                                   \
        perror("ALLOC_6D");                                             \
        fprintf(stderr,                                                 \
                "Allocation failed for " #PNTR ".  Terminating...\n");  \
        *IERR = 1; /* exit(-1); */                                     \
    }

#else

#define ALLOC_STR(PNTR, STRLEN, IERR)                   \
    if ( !PNTR )                                        \
    {                                                   \
        PNTR = (char *) malloc ( STRLEN );              \
        if (!PNTR )                                     \
        {                                               \
            perror("ALLOC_STR");                        \
            fprintf(stderr,                             \
                    "Allocation failed for " #PNTR      \
                    ".  Terminating...\n");             \
            *IERR = 1; /* exit(-1); */                 \
        }                                               \
    }                                                   \
    else                                                \
    {                                                   \
        PNTR = (char *) realloc ( PNTR, STRLEN );       \
        if (!PNTR )                                     \
        {                                               \
            perror("ALLOC_STR");                        \
            fprintf(stderr,                             \
                    "Re-allocation failed for " #PNTR   \
                    ".  Terminating...\n");             \
            *IERR = 1; /* exit(-1); */                 \
        }                                               \
    }

#define REALLOC_STR(PNTR, STRLEN, IERR)                 \
    PNTR = (char *) realloc ( PNTR, STRLEN );           \
    if (!PNTR)                                          \
    {                                                   \
        perror("REALLOC_STR");                          \
        fprintf(stderr,                                 \
                "Re-allocation failed for " #PNTR       \
                ".  Terminating...\n");                 \
        *IERR = 1; /* exit(-1); */                     \
    }

#define FREE(PNTR)                              \
    if ( PNTR )                                 \
    {                                           \
        free ( PNTR );                          \
    }

#define REALLOC_2D(PNTR, NUMX, NUMY, TYPE, IERR)                        \
    PNTR = (TYPE *)realloc(PNTR, NUMX*NUMY*sizeof(TYPE));              \
    if (!PNTR)                                                          \
    {                                                                   \
     perror("REALLOC_2D");                                                \
     fprintf(stderr,                                                    \
                 "Reallocation failed for " #PNTR ". Terminating...\n");  \
     *IERR = 1; /* exit(-1); */                                         \
     }

#define ALLOC_1D(PNTR, NUM, TYPE, IERR)                                 \
                PNTR = (TYPE *)calloc(NUM, sizeof(TYPE));               \
                if (!PNTR)                                              \
                {                                                       \
                    perror("ALLOC_1D");                                 \
                    fprintf(stderr,                                     \
                            "Allocation failed for " #PNTR ". Terminating...\n"); \
                    *IERR = 1; /* exit(-1); */                         \
                }

#define ALLOC_2D(PNTR, NUMX, NUMY, TYPE, IERR)                          \
    PNTR = (TYPE *)calloc((NUMX) * (NUMY), sizeof(TYPE));               \
    if (!PNTR)                                                          \
    {                                                                   \
     perror("ALLOC_2D");                                                \
     fprintf(stderr,                                                    \
                 "Allocation failed for " #PNTR ". Terminating...\n");  \
     *IERR = 1; /* exit(-1); */                                        \
     }

#define ALLOC_3D(PNTR, NUMX, NUMY, NUMZ, TYPE, IERR)                    \
    PNTR = (TYPE *)calloc((NUMX) * (NUMY) * (NUMZ), sizeof(TYPE));      \
    if (!PNTR)                                                          \
    {                                                                   \
        perror("ALLOC_3D");                                             \
        fprintf(stderr,                                                 \
                "Allocation failed for " #PNTR ".  Terminating...\n");  \
        *IERR = 1; /* exit(-1); */                                     \
    }

#define ALLOC_4D(PNTR, NUMW, NUMX, NUMY, NUMZ, TYPE, IERR)              \
    PNTR = (TYPE *)calloc((NUMW) * (NUMX) * (NUMY) * (NUMZ),            \
                          sizeof(TYPE));                                \
    if (!PNTR)                                                          \
    {                                                                   \
        perror("ALLOC_4D");                                             \
        fprintf(stderr,                                                 \
                "Allocation failed for " #PNTR ".  Terminating...\n");  \
        *IERR = 1; /* exit(-1); */                                     \
    }

#define ALLOC_5D(PNTR, NUMV, NUMW, NUMX, NUMY, NUMZ, TYPE, IERR)        \
    PNTR = (TYPE *)calloc((NUMV) * (NUMW) * (NUMX) * (NUMY) * (NUMZ),   \
                              sizeof(TYPE));                            \
    if (!PNTR)                                                          \
    {                                                                   \
        perror("ALLOC_5D");                                             \
        fprintf(stderr,                                                 \
                "Allocation failed for " #PNTR ".  Terminating...\n");  \
        *IERR = 1; /* exit(-1); */                                     \
    }

#define ALLOC_6D(PNTR, NUMU, NUMV, NUMW, NUMX, NUMY, NUMZ, TYPE, IERR)  \
    PNTR = (TYPE *)calloc((NUMU) * (NUMV) *(NUMW)                       \
                          * (NUMX) * (NUMY) * (NUMZ),                   \
                              sizeof(TYPE));                            \
    if (!PNTR)                                                          \
    {                                                                   \
     perror("ALLOC_6D");                                                \
     fprintf(stderr,                                                    \
                 "Allocation failed for " #PNTR ".  Terminating...\n"); \
     *IERR = 1; /* exit(-1); */                                         \
     }

#endif

// plib.c
/* root */
#define ROOT_PARA(PARA)              PARA->root
// Assuming 'PARA=para_vars'
#define ROOT                         ROOT_PARA(para_vars)

/* g_off */
#define G_OFF_PARA(PARA)             PARA->g_off
// Assuming 'PARA=para_vars'
#define G_OFF                        G_OFF_PARA(para_vars)

/* comm_snap */
#define COMM_SNAP_PARA(PARA)         PARA->comm_snap
// Assuming 'PARA=para_vars'
#define COMM_SNAP                    COMM_SNAP_PARA(para_vars)

/* comm_space */
#define COMM_SPACE_PARA(PARA)        PARA->comm_space
// Assuming 'PARA=para_vars'
#define COMM_SPACE                   COMM_SPACE_PARA(para_vars)

/* nproc */
#define NPROC_PARA(PARA)             PARA->nproc
// Assuming 'PARA=para_vars'
#define NPROC                        NPROC_PARA(para_vars)

/* iproc */
#define IPROC_PARA(PARA)             PARA->iproc
// Assuming 'PARA=para_vars'
#define IPROC                        IPROC_PARA(para_vars)

/* sproc */
#define SPROC_PARA(PARA)             PARA->sproc
// Assuming 'PARA=para_vars'
#define SPROC                        SPROC_PARA(para_vars)

/* ycomm */
#define YCOMM_PARA(PARA)             PARA->ycomm
// Assuming 'PARA=para_vars'
#define YCOMM                        YCOMM_PARA(para_vars)

/* zcomm */
#define ZCOMM_PARA(PARA)             PARA->zcomm
// Assuming 'PARA=para_vars'
#define ZCOMM                        ZCOMM_PARA(para_vars)

/* yproc */
#define YPROC_PARA(PARA)             PARA->yproc
// Assuming 'PARA=para_vars'
#define YPROC                        YPROC_PARA(para_vars)

/* zproc */
#define ZPROC_PARA(PARA)             PARA->zproc
// Assuming 'PARA=para_vars'
#define ZPROC                        ZPROC_PARA(para_vars)

/* ylop */
#define YLOP_PARA(PARA)              PARA->ylop
// Assuming 'PARA=para_vars'
#define YLOP                         YLOP_PARA(para_vars)

/* yhip */
#define YHIP_PARA(PARA)              PARA->yhip
// Assuming 'PARA=para_vars'
#define YHIP                         YHIP_PARA(para_vars)

/* zlop */
#define ZLOP_PARA(PARA)              PARA->zlop
// Assuming 'PARA=para_vars'
#define ZLOP                         ZLOP_PARA(para_vars)

/* zhip */
#define ZHIP_PARA(PARA)              PARA->zhip
// Assuming 'PARA=para_vars'
#define ZHIP                         ZHIP_PARA(para_vars)

/* thread_level */
#define THREAD_LEVEL_PARA(PARA)      PARA->thread_level
// Assuming 'PARA=para_vars'
#define THREAD_LEVEL                 THREAD_LEVEL_PARA(para_vars)

/* thread_single */
#define THREAD_SINGLE_PARA(PARA)     PARA->thread_single
// Assuming 'PARA=para_vars'
#define THREAD_SINGLE                THREAD_SINGLE_PARA(para_vars)

/* thread_funneled */
#define THREAD_FUNNELED_PARA(PARA)   PARA->thread_funneled
// Assuming 'PARA=para_vars'
#define THREAD_FUNNELED              THREAD_FUNNELED_PARA(para_vars)

/* thread_serialized */
#define THREAD_SERIALIZED_PARA(PARA) PARA->thread_serialized
// Assuming 'PARA=para_vars'
#define THREAD_SERIALIZED            THREAD_SERIALIZED_PARA(para_vars)

/* thread_multiple */
#define THREAD_MULTIPLE_PARA(PARA)   PARA->thread_multiple
// Assuming 'PARA=para_vars'
#define THREAD_MULTIPLE              THREAD_MULTIPLE_PARA(para_vars)

/* max_threads */
#define MAX_THREADS_PARA(PARA)       PARA->max_threads
// Assuming 'PARA=para_vars'
#define MAX_THREADS                  MAX_THREADS_PARA(para_vars)

/* lock */
#define LOCK_PARA(PARA)              PARA->lock
// Assuming 'PARA=para_vars'
#define LOCK                         LOCK_PARA(para_vars)

/* num_grth */
#define NUM_GRTH_PARA(PARA)          PARA->num_grth
// Assuming 'PARA=para_vars'
#define NUM_GRTH                     NUM_GRTH_PARA(para_vars)

/* firsty */
#define FIRSTY_PARA(PARA)            PARA->firsty
// Assuming 'PARA=para_vars'
#define FIRSTY                       FIRSTY_PARA(para_vars)

/* lasty */
#define LASTY_PARA(PARA)             PARA->lasty
// Assuming 'PARA=para_vars'
#define LASTY                        LASTY_PARA(para_vars)

/* firstz */
#define FIRSTZ_PARA(PARA)            PARA->firstz
// Assuming 'PARA=para_vars'
#define FIRSTZ                       FIRSTZ_PARA(para_vars)

/* lastz */
#define LASTZ_PARA(PARA)             PARA->lastz
// Assuming 'PARA=para_vars'
#define LASTZ                        LASTZ_PARA(para_vars)

/* do_nested */
#define DO_NESTED_PARA(PARA)         PARA->do_nested
// Assuming 'PARA=para_vars'
#define DO_NESTED                    DO_NESTED_PARA(para_vars)


#endif