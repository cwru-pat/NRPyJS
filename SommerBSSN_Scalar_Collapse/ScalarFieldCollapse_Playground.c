//This is being written into a C-file 
//Emscripten compiles the C code
// Step P0: Define REAL and NGHOSTS; and declare CFL_FACTOR. This header is generated in NRPy+.
#include "ScalarFieldCollapse_Playground_REAL__NGHOSTS__CFL_FACTOR.h"

#include "rfm_files/rfm_struct__declare.h"

#include "declare_Cparameters_struct.h"

#include "emscripten.h"

// All SIMD intrinsics used in SIMD-enabled C code loops are defined here:
#include "SIMD/SIMD_intrinsics.h"

// Step P1: Import needed header files
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "time.h"
#include "stdint.h" // Needed for Windows GCC 6.x compatibility
#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884L
#endif
#ifndef M_SQRT1_2
#define M_SQRT1_2 0.707106781186547524400844362104849039L
#endif
#define wavespeed 1.0 // Set CFL-based "wavespeed" to 1.0.
#define alpha_threshold (2e-3) // Value below which we rule gravitational collapse has happened

// Step P2: Declare the IDX4S(gf,i,j,k) macro, which enables us to store 4-dimensions of
//           data in a 1D array. In this case, consecutive values of "i"
//           (all other indices held to a fixed value) are consecutive in memory, where
//           consecutive values of "j" (fixing all other indices) are separated by
//           Nxx_plus_2NGHOSTS0 elements in memory. Similarly, consecutive values of
//           "k" are separated by Nxx_plus_2NGHOSTS0*Nxx_plus_2NGHOSTS1 in memory, etc.
#define IDX4S(g,i,j,k) \
( (i) + Nxx_plus_2NGHOSTS0 * ( (j) + Nxx_plus_2NGHOSTS1 * ( (k) + Nxx_plus_2NGHOSTS2 * (g) ) ) )
#define IDX4ptS(g,idx) ( (idx) + (Nxx_plus_2NGHOSTS0*Nxx_plus_2NGHOSTS1*Nxx_plus_2NGHOSTS2) * (g) )
#define IDX3S(i,j,k) ( (i) + Nxx_plus_2NGHOSTS0 * ( (j) + Nxx_plus_2NGHOSTS1 * ( (k) ) ) )
#define LOOP_REGION(i0min,i0max, i1min,i1max, i2min,i2max) \
  for(int i2=i2min;i2<i2max;i2++) for(int i1=i1min;i1<i1max;i1++) for(int i0=i0min;i0<i0max;i0++)
#define LOOP_ALL_GFS_GPS(ii) _Pragma("omp parallel for") \
  for(int (ii)=0;(ii)<Nxx_plus_2NGHOSTS_tot*NUM_EVOL_GFS;(ii)++)

// Step P3: Set UUGF and VVGF macros, as well as xx_to_Cart()
#include "boundary_conditions/gridfunction_defines.h"

// Step P4: Set xx_to_Cart(const paramstruct *restrict params,
//                     REAL *restrict xx[3],
//                     const int i0,const int i1,const int i2,
//                     REAL xCart[3]),
//           which maps xx->Cartesian via
//    {xx[0][i0],xx[1][i1],xx[2][i2]}->{xCart[0],xCart[1],xCart[2]}
#include "xx_to_Cart.h"

// Step P5: Defines set_Nxx_dxx_invdx_params__and__xx(const int EigenCoord, const int Nxx[3],
//                                       paramstruct *restrict params, REAL *restrict xx[3]),
//          which sets params Nxx,Nxx_plus_2NGHOSTS,dxx,invdx, and xx[] for
//          the chosen Eigen-CoordSystem if EigenCoord==1, or
//          CoordSystem if EigenCoord==0.
#include "set_Nxx_dxx_invdx_params__and__xx.h"

// Step P6: Include basic functions needed to impose curvilinear
//          parity and boundary conditions.
#include "boundary_conditions/CurviBC_include_Cfunctions.h"


// Step P7: Implement the algorithm for upwinding.
//          *NOTE*: This upwinding is backwards from
//          usual upwinding algorithms, because the
//          upwinding control vector in BSSN (the shift)
//          acts like a *negative* velocity.
//#define UPWIND_ALG(UpwindVecU) UpwindVecU > 0.0 ? 1.0 : 0.0

// Step P8: Include function for enforcing detgammahat constraint.
#include "enforce_detgammahat_constraint.h"

// Step P9: Find the CFL-constrained timestep
#include "find_timestep.h"

// Step P10: Declare initial data input struct:
//           stores data from initial data solver,
//           so they can be put on the numerical grid.
typedef struct __ID_inputs {
    int interp_stencil_size;
    int numlines_in_file;
    REAL *r_arr,*sf_arr,*psi4_arr,*alpha_arr;
} ID_inputs;

// Part P11: Declare all functions for setting up ScalarField initial data.
/* Routines to interpolate the ScalarField solution and convert to ADM & T^{munu}: */
#include "../nrpytutorial/ScalarField/ScalarField_interp.h"
#include "ID_scalarfield_ADM_quantities.h"
#include "ID_scalarfield_spherical.h"
#include "ID_scalarfield_xx0xx1xx2_to_BSSN_xx0xx1xx2.h"
#include "ID_scalarfield.h"

/* Next perform the basis conversion and compute all needed BSSN quantities */
#include "ID_ADM_xx0xx1xx2_to_BSSN_xx0xx1xx2__ALL_BUT_LAMBDAs.h"
#include "ID_BSSN__ALL_BUT_LAMBDAs.h"
#include "ID_BSSN_lambdas.h"



float fr[30000];
float fsf[30000];
float fpsi4[30000];
float falpha[30000];


// Step P12: Set the generic driver function for setting up BSSN initial data
void initial_data(const paramstruct *restrict params,const bc_struct *restrict bcstruct,
                  const rfm_struct *restrict rfmstruct,
                  REAL *restrict xx[3], REAL *restrict auxevol_gfs, REAL *restrict in_gfs) {
#include "set_Cparameters.h"
//We will change this to getter functions to grab the initial data from Leos initial conditions
    // Step 1: Set up ScalarField initial data
    // Step 1.a: Read ScalarField initial data from data file
    // Open the data file:
    // [Deleted to hard code in initial conditions]
    // Count the number of lines in the data file:
    
    
    // Allocate space for all data arrays:
    int numlines_in_file = 30000;
    REAL *r_arr     = (REAL *)malloc(sizeof(REAL)*numlines_in_file);
    REAL *sf_arr    = (REAL *)malloc(sizeof(REAL)*numlines_in_file);
    REAL *psi4_arr  = (REAL *)malloc(sizeof(REAL)*numlines_in_file);
    REAL *alpha_arr = (REAL *)malloc(sizeof(REAL)*numlines_in_file);

    // Read from the data file, filling in arrays
    // read_datafile__set_arrays() may be found in ScalarField/ScalarField_interp.h
    #include "../SFID.h"
    //REAL *r_arr = r_arrjs;
    //REAL *sf_arr = sf_arrjs;
    //REAL *psi4_arr = psi4_arrjs;
    //REAL *alpha_arr = alpha_arrjs;
    
    //fprintf(stderr,"R_arr %e");
    
    for(int ln =0; ln<numlines_in_file;ln++){
        r_arr[ln] = fr[ln];
        sf_arr[ln] = fsf[ln];
        psi4_arr[ln] = fpsi4[ln];
        alpha_arr[ln] = falpha[ln];
    }
    //fprintf(stderr,"INITIAL_DATA %e %e %e %e TEST\n", r_arr[29999], sf_arr[29999], psi4_arr[29999],alpha_arr[29999]);
    
    const int interp_stencil_size = 12;
    ID_inputs SF_in;
    SF_in.interp_stencil_size = interp_stencil_size;
    SF_in.numlines_in_file    = numlines_in_file;
    SF_in.r_arr               = r_arr;
    SF_in.sf_arr              = sf_arr;
    SF_in.psi4_arr            = psi4_arr;
    SF_in.alpha_arr           = alpha_arr;

    // Step 1.b: Interpolate data from data file to set BSSN gridfunctions
    ID_scalarfield(params,xx,SF_in, in_gfs);
    ID_BSSN__ALL_BUT_LAMBDAs(params,xx,SF_in, in_gfs);
    //apply_bcs_curvilinear(params, bcstruct, NUM_EVOL_GFS, evol_gf_parity, in_gfs);
    apply_bcs_sommerfeld(params, xx[3], bcstruct, NUM_EVOL_GFS, evol_gf_parity, in_gfs, rhs_gfs); 
    enforce_detgammahat_constraint(rfmstruct, params,                   in_gfs);
    ID_BSSN_lambdas(params, xx, in_gfs);
    //apply_bcs_curvilinear(params, bcstruct, NUM_EVOL_GFS, evol_gf_parity, in_gfs);
    apply_bcs_sommerfeld(params, xx[3], bcstruct, NUM_EVOL_GFS, evol_gf_parity, in_gfs, rhs_gfs); 
    enforce_detgammahat_constraint(rfmstruct, params,                   in_gfs);

    free(r_arr);
    free(sf_arr);
    free(psi4_arr);
    free(alpha_arr);
}

// Step P11: Declare function for evaluating Hamiltonian constraint (diagnostic)
#include "Hamiltonian_constraint.h"

// Step P12: Declare rhs_eval function, which evaluates BSSN RHSs
#include "rhs_eval.h"

// Step P13: Declare Ricci_eval function, which evaluates Ricci tensor
#include "Ricci_eval.h"

//#include "NRPyCritCol_regridding.h"

REAL rho_max = 0.0;
//This is where the code diverges from the tutorial.
// (These used to be defined in main function, but to make them global they were moved outside main, and getter/setter 
// functions are used to modify them)
// Step P15: Declare global pointers and variables to be referenced in getter function.
//           This must be done in order to access variables created in the initialize
//           function in the stepfoward function.
int dimensions[4] = {28, 8, 2, 18}; //#Filler values that are updated later on

    
int arrNGHOSTS[4];
paramstruct params_p;
rfm_struct rfmstruct_p;
bc_struct bcstruct_p;
REAL N_final_p, output_every_N_p, dt_p;
//Contains scalar & metric fields, but also psi_4, so dont need the psi4 ones.
// Make resolution + cfl factor global variables. Generally make most variables global. 
REAL *y_n_gfs_p, *auxevol_gfs_p, *k_odd_gfs_p, *k_even_gfs_p, *y_nplus1_running_total_gfs_p, *xx_p[3];

// initialize() function:
// Step 0: Set up grid structure, allocate memory for gridfunctions, set up coordinates
// Step 1: Set up initial data to an exact solutiondef
// Step 2: Initialize global pointers and variables

// Main function replaced by initialize and getter functions. Emscripten connects to Javascript.
int EMSCRIPTEN_KEEPALIVE initialize(REAL CFL_FACTOR){
    
    //fprintf(stderr,"Running initialize!!\n");
    paramstruct params;
#include "set_Cparameters_default.h"

    // Step 0a: Set up numerical grid structure, first in space...
    int n1 = dimensions[0];
    int n2 = dimensions[1];
    int n3 = dimensions[2];
    const int Nxx[3] = {n1,n2,n3};

    // Step 0b: Set free parameters, overwriting Cparameters defaults
    //          by hand or with command-line input, as desired.
#include "free_parameters.h"

    // Step 0c: Uniform coordinate grids are stored to *xx[3]
    REAL *xx[3];
    // Step 0c.i: Set bcstruct
    bc_struct bcstruct;
    {
        int EigenCoord = 1;
        // Step 0c.ii: Call set_Nxx_dxx_invdx_params__and__xx(), which sets
        //             params Nxx,Nxx_plus_2NGHOSTS,dxx,invdx, and xx[] for the
        //             chosen Eigen-CoordSystem.
        
        set_Nxx_dxx_invdx_params__and__xx(EigenCoord, Nxx, &params, xx);
        
        // Step 0c.iii: Set Nxx_plus_2NGHOSTS_tot
#include "set_Cparameters-nopointer.h"
        const int Nxx_plus_2NGHOSTS_tot = Nxx_plus_2NGHOSTS0*Nxx_plus_2NGHOSTS1*Nxx_plus_2NGHOSTS2;
        // Step 0d: Find ghostzone mappings; set up bcstruct
#include "boundary_conditions/driver_bcstruct.h"
        // Step 0d.i: Free allocated space for xx[][] array
        for(int i=0;i<3;i++) free(xx[i]);
    }
    
    // Step 0e: Call set_Nxx_dxx_invdx_params__and__xx(), which sets
    //          params Nxx,Nxx_plus_2NGHOSTS,dxx,invdx, and xx[] for the
    //          chosen (non-Eigen) CoordSystem.
    int EigenCoord = 0;
    set_Nxx_dxx_invdx_params__and__xx(EigenCoord, Nxx, &params, xx);

    // Step 0f: Set all C parameters "blah" for params.blah, including
    //          Nxx_plus_2NGHOSTS0 = params.Nxx_plus_2NGHOSTS0, etc.
    #include "set_Cparameters-nopointer.h"
    const int Nxx_plus_2NGHOSTS_tot = Nxx_plus_2NGHOSTS0*Nxx_plus_2NGHOSTS1*Nxx_plus_2NGHOSTS2;

    // Step 0g: Time coordinate parameters
    const REAL t_final = 16; // Run until the user says stop, but let them be able to pause and continue, with reset button.

    // Step 0h: Set timestep based on smallest proper distance between gridpoints and CFL factor
    REAL dt = find_timestep(&params, xx);
    
    N_final_p = (int)(t_final / dt + 0.5); // The number of points in time.
                                             // Add 0.5 to account for C rounding down
                                             // typecasts to integers.
    //The Collision Code:
    //REAL out_approx_every_t = 0.2;
    //output_every_N_p = (int)(out_approx_every_t*((REAL)N_final_p)/t_final);
    //The Collapse code is below:
    //Every Nth step writes it to a file. For our version, instead of writing to a file just 
    //            1. Stores data 2. animates
    // We probably wont need to store data.
    REAL output_every_N = 20;
    
    

    // Step 0i: Error out if the number of auxiliary gridfunctions outnumber evolved gridfunctions.
    //              This is a limitation of the RK method. You are always welcome to declare & allocate
    //              additional gridfunctions by hand.
    
    if(NUM_AUX_GFS > NUM_EVOL_GFS) {
        fprintf(stderr,"Error: NUM_AUX_GFS > NUM_EVOL_GFS. Either reduce the number of auxiliary gridfunctions,\n");
        fprintf(stderr,"       or allocate (malloc) by hand storage for *diagnostic_output_gfs. \n");
        exit(1);
    }

    // Step 0j: Allocate memory for gridfunctions
#include "MoLtimestepping/RK_Allocate_Memory.h"
    REAL *restrict auxevol_gfs = (REAL *)malloc(sizeof(REAL) * NUM_AUXEVOL_GFS * Nxx_plus_2NGHOSTS_tot);

    // Step 0k: Set up precomputed reference metric arrays
    // Step 0k.i: Allocate space for precomputed reference metric arrays.
#include "rfm_files/rfm_struct__malloc.h"

    // Step 0k.ii: Define precomputed reference metric arrays.
    {
    #include "set_Cparameters-nopointer.h"
    #include "rfm_files/rfm_struct__define.h"
    }

    // Step 1a: Set up initial data to an exact solution
    initial_data(&params,&bcstruct, &rfmstruct, xx, auxevol_gfs, y_n_gfs);

    // Step 1b: Apply boundary conditions, as initial data
    //          are sometimes ill-defined in ghost zones.
    //          E.g., spherical initial data might not be
    //          properly defined at points where r=-1.
    //apply_bcs_curvilinear(&params, &bcstruct, NUM_EVOL_GFS,evol_gf_parity, y_n_gfs);
    apply_bcs_sommerfeld(&params, &xx[3], &bcstruct, NUM_EVOL_GFS, evol_gf_parity, y_n_gfs, &rhs_gfs); 
    enforce_detgammahat_constraint(&rfmstruct, &params, y_n_gfs);

    //Step 2: Assign pointers/Initialize global variables
    arrNGHOSTS[0] = NGHOSTS;
    arrNGHOSTS[1] = Nxx_plus_2NGHOSTS0;
    arrNGHOSTS[2] = Nxx_plus_2NGHOSTS1;
    arrNGHOSTS[3] = Nxx_plus_2NGHOSTS2;
    dt_p = dt; 
    //These appear to be all the pointers needed 
    params_p = params;
    rfmstruct_p = rfmstruct;
    bcstruct_p = bcstruct;
    y_n_gfs_p = y_n_gfs;
    auxevol_gfs_p = auxevol_gfs;
    k_odd_gfs_p = k_odd_gfs;
    k_even_gfs_p = k_even_gfs;
    y_nplus1_running_total_gfs_p = y_nplus1_running_total_gfs;
    xx_p[0]=xx[0];
    xx_p[1]=xx[1];
    xx_p[2]=xx[2];
    //printf("%f",xx_p[0]);

    return 0;
}

// stepForward() function:
// Step 1: Define and initialize variables from initialize() function so they can be used in the RK-like Method
//         of Lines timestepping algorithm
// Step 2: Step forward one timestep (t -> t+dt) in time using chosen RK-like MoL timestepping algorithm
// Second half of main, does one single timescript that Javascript loops
void EMSCRIPTEN_KEEPALIVE stepForward(j){

    // Step 1: Redefine and initialize variables. In order to call each time-step one by one, we need to redefine
    //         some variables used in the MoL timestepping algorithm with the saved values from the initialization
    //         step
    int Nxx_plus_2NGHOSTS0 = arrNGHOSTS[1];
    int Nxx_plus_2NGHOSTS1 = arrNGHOSTS[2];
    int Nxx_plus_2NGHOSTS2 = arrNGHOSTS[3];
    int Nxx_plus_2NGHOSTS_tot = Nxx_plus_2NGHOSTS0*Nxx_plus_2NGHOSTS1*Nxx_plus_2NGHOSTS2;
    int N_final = N_final_p;
    int output_every_N = output_every_N_p;
    //REAL *restrict diagnostic_output_gfs = diagnostic_output_gfs_p;
    REAL dt = dt_p;
    paramstruct params = params_p;
    rfm_struct rfmstruct = rfmstruct_p;
    bc_struct bcstruct = bcstruct_p;
    REAL *y_n_gfs = y_n_gfs_p;
    REAL *restrict auxevol_gfs = auxevol_gfs_p;
    REAL *k_odd_gfs = k_odd_gfs_p;
    REAL *k_even_gfs = k_even_gfs_p;
    REAL *restrict y_nplus1_running_total_gfs = y_nplus1_running_total_gfs_p;
    REAL *xx[3];
    xx[0]=xx_p[0];
    xx[1]=xx_p[1];
    xx[2]=xx_p[2];
    //if(j){
    //    size_t len = sizeof(y_n_gfs_p)/sizeof(y_n_gfs_p[0]);
    //    fprintf(stderr, "The length of the array is: %zu\n", len);
    //}
#include "boundary_conditions/driver_bcstruct.h"

    // Step 2: Step forward one timestep (t -> t+dt) in time using
    //           chosen RK-like MoL timestepping algorithm
    // Step 3.a: Output 2D data file periodically, for visualization
  //  if(n%output_every_N == 0) {
        // Evaluate Hamiltonian constraint violation
   //     Hamiltonian_constraint(&rfmstruct, &params, y_n_gfs,auxevol_gfs, diagnostic_output_gfs);

        //char filename[100];
        //sprintf(filename,"out%d-%08d.txt",Nxx[0],n);
        //const int i1mid=Nxx_plus_2NGHOSTS1/2;
        //const int i2mid=Nxx_plus_2NGHOSTS2/2;
        //FILE *fp = fopen(filename, "w");
        //for( int i0=NGHOSTS;i0<Nxx_plus_2NGHOSTS0-NGHOSTS;i0++) {
        //    const int idx  = IDX3S(i0,i1mid,i2mid);
        //    const REAL xx0 = xx[0][i0];
        //    REAL xCart[3];
        //    xx_to_Cart(&params,xx,i0,i1mid,i2mid,xCart);
        //    const REAL rr = sqrt( xCart[0]*xCart[0] + xCart[1]*xCart[1] + xCart[2]*xCart[2] );
        //    fprintf(fp,"%e %e %e %e %e %e %e\n",xx0,rr,
        //            y_n_gfs[IDX4ptS(SFGF,idx)],y_n_gfs[IDX4ptS(SFMGF,idx)],
        //            y_n_gfs[IDX4ptS(ALPHAGF,idx)],y_n_gfs[IDX4ptS(CFGF,idx)],
        //            log10(fabs(diagnostic_output_gfs[IDX4ptS(HGF,idx)])));
        //}
        //fclose(fp);
    //}

#include "MoLtimestepping/RK_MoL.h"
}



// Getter functions used to access the data in javascript after compliling.

// getNFinal(): returns final time-step value
REAL EMSCRIPTEN_KEEPALIVE getNFinal(){
    return N_final_p;
}

// setNFinal(): sets the final time-step value
void EMSCRIPTEN_KEEPALIVE setNFinal(int i){
    N_final_p = i;
}

// getdt(): gets the value for dt
REAL EMSCRIPTEN_KEEPALIVE getdt(){
    return dt_p;
}

// getNGHOSTS(): returns desired dimensions including ghosts shells
REAL EMSCRIPTEN_KEEPALIVE getNGHOSTS(int i){
    return arrNGHOSTS[i];
}

// getFxVal(): returns desired function value at specfic point in space, index i is determined from the
//             IDX3S and IDX4ptS arrays
REAL EMSCRIPTEN_KEEPALIVE getFxVal(int i){
    return y_n_gfs_p[i];
}

// getIDX3S(): returns IDX3S index for a given point
REAL EMSCRIPTEN_KEEPALIVE getIDX3S(int i0, int i1, int i2){
    int Nxx_plus_2NGHOSTS0 = arrNGHOSTS[1];
    int Nxx_plus_2NGHOSTS1 = arrNGHOSTS[2];
    int Nxx_plus_2NGHOSTS2 = arrNGHOSTS[3];
    return IDX3S(i0,i1,i2);
}

// getIDX4ptS(): returns IDX4ptS index to be used in getFxVal() using information from getIDX3S()
REAL EMSCRIPTEN_KEEPALIVE getIDX4ptS(int fx, REAL idx){
    int Nxx_plus_2NGHOSTS0 = arrNGHOSTS[1];
    int Nxx_plus_2NGHOSTS1 = arrNGHOSTS[2];
    int Nxx_plus_2NGHOSTS2 = arrNGHOSTS[3];
    return IDX4ptS(fx,idx);
}

// getCart(): returns the cartesian version of a given spherical coordinate
REAL EMSCRIPTEN_KEEPALIVE getCart(int i0, int i1, int i2, int index){
    REAL xCart[3];
    xx_to_Cart(&params_p,xx_p,i0,i1,i2,xCart);
    return xCart[index];
}

// setDim(): setter function allowing the user to change the resolution of the simulation
void EMSCRIPTEN_KEEPALIVE setDim(int dim1, int dim2, int dim3){
    dimensions[0] = dim1;
    dimensions[1] = dim2;
    dimensions[2] = dim3;
    dimensions[3]= dim1/4;
}

// getDim(): getter function for the dimensions/resolution of the simulation
int EMSCRIPTEN_KEEPALIVE getDim(int i){
    return dimensions[i];
}

// setInitials(): setter functions for the initial conditions
void EMSCRIPTEN_KEEPALIVE setInitR(float i, int j){
    fr[j] = i;
}
void EMSCRIPTEN_KEEPALIVE setInitSf(float i, int j){
    fsf[j] = i;
}
void EMSCRIPTEN_KEEPALIVE setInitPsi4(float i, int j){
    fpsi4[j] = i;
}
void EMSCRIPTEN_KEEPALIVE setInitAlpha(float i, int j){
    falpha[j] = i;
}
// A very crude attempt that hopefully returns the scalar field values
//void EMSCRIPTEN_KEEPALIVE getSf(int i){
//    return sf[i];
//} (It did not return the scalar field)

// runsim(): function that runs the initialize function with a specific resolution and CFL factor
void EMSCRIPTEN_KEEPALIVE runsim(){
    CFL_FACTOR = 1.0;
    
    fprintf(stderr,"Running Sim3!\n");
    
    initialize(CFL_FACTOR);
    
    fprintf(stderr,"Finished Sim!\n");
}
