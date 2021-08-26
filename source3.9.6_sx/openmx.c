/*****************************************************************************

  Ver. 3.9 (26/July/2019)

  OpenMX (Open source package for Material eXplorer) is a program package
  for linear scaling density functional calculations of large-scale materials.
  Almost time-consuming parts can be performed in O(N) operations where N is
  the number of atoms (or basis orbitals). Thus, the program might be useful
  for studies of large-scale materials.
  The distribution of this program package follows the practice of
  the GNU General Public Licence (GPL).

  OpenMX is based on 

   *  Local density and generalized gradient approximation (LDA, LSDA, GGA)
      to the exchange-corellation term
   *  Norm-conserving pseudo potentials
   *  Variationally optimized pseudo atomic basis orbitals
   *  Solution of Poisson's equation using FFT
   *  Evaluation of two-center integrals using Fourier transformation
   *  Evaluation of three-center integrals using fixed real space grids
   *  Simple mixing, direct inversion in the interative subspace (DIIS),
      and Guaranteed-reduction Pulay's methods for SCF calculations.
   *  Solution of the eigenvalue problem using O(N) methods
   *  ...

  See also our website (http://www.openmx-square.org/)
  for recent developments.


    **************************************************************
     Copyright

     Taisuke Ozaki

     Present (23/Sep./2019) official address

       Institute for Solid State Physics, University of Tokyo,
       Kashiwanoha 5-1-5, Kashiwa, Chiba 277-8581, Japan

       e-mail: t-ozaki@issp.u-tokyo.ac.jp
    **************************************************************
 
*****************************************************************************/

/**********************************************************************
  openmx.c:

     openmx.c is the main routine of OpenMX.

  Log of openmx.c:

     5/Oct/2003  Released by T.Ozaki

***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
/*  stat section */
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
/*  end stat section */
#include "openmx_common.h"
#include "mpi.h"
#include <omp.h>

#include "tran_prototypes.h"
#include "tran_variables.h"


int main(int argc, char *argv[]) 
{ 
  static int numprocs,myid;
  static int MD_iter,i,j,po,ip;
  static char fileMemory[YOUSO10]; 
  double TStime,TEtime;
  MPI_Comm mpi_comm_parent;

  /* MPI initialize */

  mpi_comm_level1 = MPI_COMM_WORLD; 
  MPI_COMM_WORLD1 = MPI_COMM_WORLD; 

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD1,&numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD1,&myid);
  NUMPROCS_MPI_COMM_WORLD = numprocs;
  MYID_MPI_COMM_WORLD = myid;
  Num_Procs = numprocs;

  /* check if OpenMX was called by MPI_spawn. */

  MPI_Comm_get_parent(&mpi_comm_parent);
  if (mpi_comm_parent!=MPI_COMM_NULL) MPI_spawn_flag = 1;
  else                                MPI_spawn_flag = 0;

  /* for measuring elapsed time */

  dtime(&TStime);

  /* check argv */

  if (argc==1){

    if (myid==Host_ID) printf("\nCould not find an input file.\n\n");
    MPI_Finalize(); 
    exit(0);
  } 

  /* initialize Runtest_flag */

  Runtest_flag = 0;

  /****************************************************
    ./openmx -nt # 

    specifies the number of threads in parallelization
    by OpenMP
  ****************************************************/
  
  openmp_threads_num = 1; /* default */

  po = 0;
  if (myid==Host_ID){
    for (i=0; i<argc; i++){
      if ( strcmp(argv[i],"-nt")==0 ){
        po = 1;
        ip = i;
      }
    }
  }

  MPI_Bcast(&po, 1, MPI_INT, Host_ID, MPI_COMM_WORLD1);
  MPI_Bcast(&ip, 1, MPI_INT, Host_ID, MPI_COMM_WORLD1);

  if ( (argc-1)<(ip+1) ){
    if (myid==Host_ID){
      printf("cannot find the number of threads\n");
    }
    MPI_Finalize();
    exit(0);
  }

  if ( po==1 ){
    openmp_threads_num = atoi(argv[ip+1]);

    if (openmp_threads_num<=0){ 
      if (myid==Host_ID){
        printf("check the number of threads\n");
      }
      MPI_Finalize();
      exit(0);
    }
  }

  omp_set_num_threads(openmp_threads_num);  

  if (myid==Host_ID){
    printf("\nThe number of threads in each node for OpenMP parallelization is %d.\n\n",openmp_threads_num);
  }

  /****************************************************
    ./openmx -show directory 

    showing PAOs and VPS used in input files stored in 
    "directory".
  ****************************************************/

  if (strcmp(argv[1],"-show")==0){
    Show_DFT_DATA(argv);
    exit(0);
  }

  /****************************************************
    ./openmx -maketest

    making of *.out files in order to check whether 
    OpenMX normally runs on many platforms or not.
  ****************************************************/

  if ( (argc==2 || argc==3) && strcmp(argv[1],"-maketest")==0){
    Maketest("S",argc,argv);
    exit(0);
  }

  /****************************************************
   ./openmx -runtest

   check whether OpenMX normally runs on many platforms
   or not by comparing the stored *.out and generated
   *.out on your machine.
  ****************************************************/

  if ( strcmp(argv[1],"-runtest")==0){
    Runtest("S",argc,argv);
  }

  /****************************************************
   ./openmx -maketestL

    making of *.out files in order to check whether 
    OpenMX normally runs for relatively large systems
    on many platforms or not
  ****************************************************/

  if ( (argc==2 || argc==3) && strcmp(argv[1],"-maketestL")==0){
    Maketest("L",argc,argv);
    exit(0);
  }

  /****************************************************
   ./openmx -maketestL2

    making of *.out files in order to check whether 
    OpenMX normally runs for large systems on many 
    platforms or not
  ****************************************************/

  if ( (argc==2 || argc==3) && strcmp(argv[1],"-maketestL2")==0){
    Maketest("L2",argc,argv);
    exit(0);
  }

  /****************************************************
   ./openmx -maketestL3

    making of *.out files in order to check whether 
    OpenMX normally runs for large systems on many 
    platforms or not
  ****************************************************/

  if ( (argc==2 || argc==3) && strcmp(argv[1],"-maketestL3")==0){
    Maketest("L3",argc,argv);
    exit(0);
  }

  /****************************************************
   ./openmx -runtestL

   check whether OpenMX normally runs for relatively
   large systems on many platforms or not by comparing
   the stored *.out and generated *.out on your machine.
  ****************************************************/

  if (strcmp(argv[1],"-runtestL")==0){
    Runtest("L",argc,argv);
  }

  /****************************************************
   ./openmx -runtestL2

   check whether OpenMX normally runs for large systems 
   on many platforms or not by comparing the stored *.out 
   and generated *.out on your machine.
  ****************************************************/

  if (strcmp(argv[1],"-runtestL2")==0){
    Runtest("L2",argc,argv);
  }

  /****************************************************
   ./openmx -runtestL3

   check whether OpenMX normally runs for small, medium size, 
   and large systems on many platforms or not by comparing 
   the stored *.out and generated *.out on your machine.
  ****************************************************/

  if (strcmp(argv[1],"-runtestL3")==0){
    Runtest("L3",argc,argv);
  }

  /*******************************************************
   check memory leak by monitoring actual used memory size
  *******************************************************/

  if ( (argc==2 || argc==3) && strcmp(argv[1],"-mltest")==0){
    Memory_Leak_test(argc,argv);
    exit(0);
  }

  /****************************************************
   ./openmx -maketestG

    making of *.out files in order to check whether 
    OpenMX normally runs for geometry optimization
    on many platforms or not
  ****************************************************/

  if ( (argc==2 || argc==3) && strcmp(argv[1],"-maketestG")==0){
    Maketest("G",argc,argv);
    exit(0);
  }

  /****************************************************
   ./openmx -maketestC

    making of *.out files in order to check whether 
    OpenMX normally runs for geometry optimization
    on many platforms or not
  ****************************************************/

  if ( (argc==2 || argc==3) && strcmp(argv[1],"-maketestC")==0){
    Maketest("C",argc,argv);
    exit(0);
  }

  /****************************************************
   ./openmx -runtestG

   check whether OpenMX normally runs for geometry 
   optimization on many platforms or not by comparing
   the stored *.out and generated *.out on your machine.
  ****************************************************/

  if (strcmp(argv[1],"-runtestG")==0){
    Runtest("G",argc,argv);
  }

  /****************************************************
   ./openmx -runtestC

   check whether OpenMX normally runs for simultaneous 
   optimization for cell and geometry on many platforms 
   or not by comparing the stored *.out and generated
   *.out on your machine.
  ****************************************************/

  if (strcmp(argv[1],"-runtestC")==0){
    Runtest("C",argc,argv);
  }

  /****************************************************
   ./openmx -maketestWF

    making of *.out files in order to check whether 
    OpenMX normally runs for generation of Wannier 
    functions on many platforms or not
  ****************************************************/

  if ( (argc==2 || argc==3) && strcmp(argv[1],"-maketestWF")==0){
    Maketest("WF",argc,argv);
    exit(0);
  }

  /****************************************************
   ./openmx -runtestWF

   check whether OpenMX normally runs for generating 
   Wannier functions on many platforms or not by comparing
   the stored *.out and generated *.out on your machine.
  ****************************************************/

  if (strcmp(argv[1],"-runtestWF")==0){
    Runtest("WF",argc,argv);
  }

  /****************************************************
   ./openmx -maketestNEGF

    making of *.out files in order to check whether 
    OpenMX normally runs for NEGF calculations 
    on many platforms or not
  ****************************************************/

  if ( (argc==2 || argc==3) && strcmp(argv[1],"-maketestNEGF")==0){
    Maketest("NEGF",argc,argv);
    exit(0);
  }

  /****************************************************
   ./openmx -runtestNEGF

    check whether OpenMX normally runs for NEGF calculations 
    on many platforms or not
  ****************************************************/

  if (strcmp(argv[1],"-runtestNEGF")==0){
    Runtest("NEGF",argc,argv);
    MPI_Finalize();
    exit(0);
  }

  /********************************************************************
   ./openmx -maketestCDDF

    making of *.df_re and *.df_im files in order to check whether
    OpenMX normally runs for CDDF calculations on many platforms or not
  *********************************************************************/

  if ( (argc==2 || argc==3) && strcmp(argv[1],"-maketestCDDF")==0){
    Maketest("CDDF",argc,argv);
    exit(0);
  }

  /****************************************************
   ./openmx -runtestCDDF

    check whether OpenMX normally runs for CDDF calculations 
    on many platforms or not
  ****************************************************/

  if (strcmp(argv[1],"-runtestCDDF")==0){
    Runtest("CDDF",argc,argv);
    MPI_Finalize();
    exit(0);
  }

  /********************************************************************
   ./openmx -maketestDCLNO

    making of *.out files in order to check whether 
    OpenMX normally runs for DC-LNO calculations
    on many platforms or not
  *********************************************************************/

  if ( (argc==2 || argc==3) && strcmp(argv[1],"-maketestDCLNO")==0){
    Maketest("DCLNO",argc,argv);
    exit(0);
  }

  /****************************************************
   ./openmx -runtestDCLNO

   check whether OpenMX normally runs for DCLNO calculations 
   on many platforms or not
  ****************************************************/

  if (strcmp(argv[1],"-runtestDCLNO")==0){
    Runtest("DCLNO",argc,argv);
    MPI_Finalize();
    exit(0);
  }

  /*******************************************************
   check consistency between analytic and numerical forces
  *******************************************************/

  if ( (argc==3 || argc==4) && strcmp(argv[1],"-forcetest")==0){

    if      (strcmp(argv[2],"0")==0) force_flag = 0; 
    else if (strcmp(argv[2],"1")==0) force_flag = 1; 
    else if (strcmp(argv[2],"2")==0) force_flag = 2; 
    else if (strcmp(argv[2],"3")==0) force_flag = 3; 
    else if (strcmp(argv[2],"4")==0) force_flag = 4; 
    else if (strcmp(argv[2],"5")==0) force_flag = 5; 
    else if (strcmp(argv[2],"6")==0) force_flag = 6; 
    else if (strcmp(argv[2],"7")==0) force_flag = 7;
    else if (strcmp(argv[2],"8")==0) force_flag = 8;
    else {
      printf("unsupported flag for -forcetest\n");
      exit(0);
    }

    Force_test(argc,argv);
    exit(0);
  }

  /*********************************************************
   check consistency between analytic and numerical stress
  *********************************************************/

  if ( (argc==3 || argc==4) && strcmp(argv[1],"-stresstest")==0){

    if      (strcmp(argv[2],"0")==0) stress_flag = 0; 
    else if (strcmp(argv[2],"1")==0) stress_flag = 1; 
    else if (strcmp(argv[2],"2")==0) stress_flag = 2; 
    else if (strcmp(argv[2],"3")==0) stress_flag = 3; 
    else if (strcmp(argv[2],"4")==0) stress_flag = 4; 
    else if (strcmp(argv[2],"5")==0) stress_flag = 5; 
    else if (strcmp(argv[2],"6")==0) stress_flag = 6; 
    else if (strcmp(argv[2],"7")==0) stress_flag = 7;
    else if (strcmp(argv[2],"8")==0) stress_flag = 8;
    else {
      printf("unsupported flag for -stresstest\n");
      exit(0);
    }

    Stress_test(argc,argv);
    MPI_Finalize();
    exit(0);
  }

  /*******************************************************
    check the NEB calculation or not, and if yes, go to 
    the NEB calculation.
  *******************************************************/

  if (neb_check(argv)) neb(argc,argv);

  /*******************************************************
   allocation of CompTime and show the greeting message 
  *******************************************************/

  CompTime = (double**)malloc(sizeof(double*)*numprocs); 
  for (i=0; i<numprocs; i++){
    CompTime[i] = (double*)malloc(sizeof(double)*30); 
    for (j=0; j<30; j++) CompTime[i][j] = 0.0;
  }

  if (myid==Host_ID){  
    printf("\n*******************************************************\n"); 
    printf("*******************************************************\n"); 
    printf(" Welcome to OpenMX   Ver. %s                           \n",Version_OpenMX); 
    printf(" Copyright (C), 2002-2019, T. Ozaki                    \n"); 
    printf(" OpenMX comes with ABSOLUTELY NO WARRANTY.             \n"); 
    printf(" This is free software, and you are welcome to         \n"); 
    printf(" redistribute it under the constitution of the GNU-GPL.\n");
    printf("*******************************************************\n"); 
    printf("*******************************************************\n\n"); 
  } 

  Init_List_YOUSO();
  remake_headfile = 0;
  ScaleSize = 1.2; 

  /****************************************************
                   Read the input file
  ****************************************************/

  init_alloc_first();

  CompTime[myid][1] = readfile(argv);
  MPI_Barrier(MPI_COMM_WORLD1);

  /* initialize PrintMemory routine */

  sprintf(fileMemory,"%s%s.memory%i",filepath,filename,myid);
  PrintMemory(fileMemory,0,"init"); 
  PrintMemory_Fix();

  /* initialize */

  init();

  /* for DFTD-vdW by okuno */
  /* for version_dftD by Ellner*/
  if(dftD_switch==1 && version_dftD==2) DFTDvdW_init();
  if(dftD_switch==1 && version_dftD==3) DFTD3vdW_init();

  /* check "-mltest2" mode */

  po = 0;
  if (myid==Host_ID){
    for (i=0; i<argc; i++){
      if ( strcmp(argv[i],"-mltest2")==0 ){
        po = 1;
        ip = i;
      }
    }
  }

  MPI_Bcast(&po, 1, MPI_INT, Host_ID, MPI_COMM_WORLD1);
  MPI_Bcast(&ip, 1, MPI_INT, Host_ID, MPI_COMM_WORLD1);

  if ( po==1 ) ML_flag = 1;
  else         ML_flag = 0;

  /* check "-forcetest2" mode */

  po = 0;
  if (myid==Host_ID){
    for (i=0; i<argc; i++){
      if ( strcmp(argv[i],"-forcetest2")==0 ){
        po = 1;
        ip = i;
      }
    }
  }

  MPI_Bcast(&po, 1, MPI_INT, Host_ID, MPI_COMM_WORLD1);
  MPI_Bcast(&ip, 1, MPI_INT, Host_ID, MPI_COMM_WORLD1);

  if ( po==1 ){
    force_flag = atoi(argv[ip+1]);
    ForceConsistency_flag = 1;
  }

  /* check force consistency 
     the number of processes 
     should be less than 2.
  */

  if (ForceConsistency_flag==1){

    Check_Force(argv);
    CompTime[myid][20] = OutData(argv[1]);
    Merge_LogFile(argv[1]);
    Free_Arrays(0);
    MPI_Finalize();
    exit(0); 
    return 0;
  }

  /* check "-stresstest2" mode */

  po = 0;
  if (myid==Host_ID){


    for (i=0; i<argc; i++){
      if ( strcmp(argv[i],"-stresstest2")==0 ){

        po = 1;
        ip = i;
      }
    }
  }

  MPI_Bcast(&po, 1, MPI_INT, Host_ID, MPI_COMM_WORLD1);
  MPI_Bcast(&ip, 1, MPI_INT, Host_ID, MPI_COMM_WORLD1);

  if ( po==1 ){
    stress_flag = atoi(argv[ip+1]);
    StressConsistency_flag = 1;
  }

  /* check stress consistency 
     the number of processes 
     should be less than 2.
  */

  if (StressConsistency_flag==1){
    Check_Stress(argv);
    CompTime[myid][20] = OutData(argv[1]);
    Merge_LogFile(argv[1]);
    Free_Arrays(0);
    MPI_Finalize();
    exit(0); 
    return 0;
  }

  /****************************************************
      SCF-DFT calculations, MD and geometrical
      optimization.
  ****************************************************/

  MD_iter = 1;
  Temp_MD_iter = 1;

  do {

    if (MD_switch==12)
      CompTime[myid][2] += truncation(1,1);  /* EvsLC */
    else if (MD_cellopt_flag==1)
      CompTime[myid][2] += truncation(1,1);  /* cell optimization */
    else 
      CompTime[myid][2] += truncation(MD_iter,1);

#ifdef _NEC
  MPI_Barrier(mpi_comm_level1);
#endif

    if (ML_flag==1 && myid==Host_ID) Get_VSZ(MD_iter);

    if (Solver==4) {
      TRAN_Calc_GridBound( mpi_comm_level1, atomnum, WhatSpecies, Spe_Atom_Cut1,
                           Ngrid1, Grid_Origin, Gxyz, tv, gtv, rgtv, Left_tv, Right_tv );

      /* output: TRAN_region[], TRAN_grid_bound */
    }

    if (Solver!=4 || TRAN_SCF_skip==0){

      CompTime[myid][3] += DFT(MD_iter,(MD_iter-1)%orbitalOpt_per_MDIter+1);
      iterout(MD_iter+MD_Current_Iter,MD_TimeStep*(MD_iter+MD_Current_Iter-1),filepath,filename);

      /* MD or geometry optimization */
      if (ML_flag==0) CompTime[myid][4] += MD_pac(MD_iter,argv[1]);
    }
    else{
      MD_Opt_OK = 1;
    }

    MD_iter++;
    Temp_MD_iter++;

  } while(MD_Opt_OK==0 && (MD_iter+MD_Current_Iter)<=MD_IterNumber);

  if ( TRAN_output_hks ) {
     /* left is dummy */
     TRAN_RestartFile(mpi_comm_level1, "write","left",filepath,TRAN_hksoutfilename);
  }

  /****************************************************
               calculate Voronoi charge
  ****************************************************/
 
  if (Voronoi_Charge_flag==1) Voronoi_Charge();

  /****************************************************
        calculate Voronoi orbital magnetic moment
  ****************************************************/
 
  if (Voronoi_OrbM_flag==1) Voronoi_Orbital_Moment();

  /****************************************************
        output analysis of decomposed energies
  ****************************************************/

  if (Energy_Decomposition_flag==1) Output_Energy_Decomposition();

  /****************************************************
  making of a file *.frac for the fractional coordinates
  ****************************************************/

  Make_FracCoord(argv[1]);

  /****************************************************
   generate Wannier functions added by Hongming Weng
  ****************************************************/

  /* hmweng */
  if(Wannier_Func_Calc){
    if (myid==Host_ID) printf("Calling Generate_Wannier...\n");fflush(0);

    Generate_Wannier(argv[1]);
  }

  /****************************************************
      population analysis based on atomic orbitals 
              resembling Wannier functions
  ****************************************************/

  if(pop_anal_aow_flag){
    if (myid==Host_ID) printf("Population analysis based on atomic orbitals resembling Wannier functions\n");fflush(0);

    Population_Analysis_Wannier2(argv);
  }

  /*********************************************************
     Electronic transport calculations based on NEGF:
     transmission, current, eigen channel analysis, and
     real space analysis of current
  *********************************************************/

  if (Solver==4 && TRAN_analysis==1) {

    /* if SCF is skipped, calculate values of basis functions on each grid */ 
    if (1<=TRAN_SCF_skip) i = Set_Orbitals_Grid(0); 

    if (SpinP_switch==3) {
      TRAN_Main_Analysis_NC( mpi_comm_level1, argc, argv, Matomnum, M2G,  
                             GridN_Atom, GridListAtom, CellListAtom, 
                             Orbs_Grid, TNumGrid );
    }
    else {
      TRAN_Main_Analysis( mpi_comm_level1, argc, argv, Matomnum, M2G,  
                          GridN_Atom, GridListAtom, CellListAtom, 
                          Orbs_Grid, TNumGrid );
    }
  }

  /*********************************************************
   calculations of core level spectra:

   CLE_Type =-1; NONE 
   CLE_Type = 0; XANES0: single particle calculation
   CLE_Type = 1; XANES1: core excitation by delta-SCF
   CLE_Type = 2; XANES2: core excitation with a valence 
                         excitation by delta-SCF
  **********************************************************/

  /* calc. of overlap matrix with the position operator */

  if (0<=CLE_Type){
    Set_OLP_p(OLP_p);

    /* XANES:  single particle calculation */
    if (CLE_Type==0){
      /* XANES0(); */
    }

  }

  /****************************************************
                  Making of output files
  ****************************************************/

  if (OutData_bin_flag) 
    CompTime[myid][20] = OutData_Binary(argv[1]);
  else 
    CompTime[myid][20] = OutData(argv[1]);

  /****************************************************
    write connectivity, Hamiltonian, overlap, density
    matrices, and etc. to a file, filename.scfout 
  ****************************************************/

  if (HS_fileout==1) SCF2File("write",argv[1]);

  /* elapsed time */

  dtime(&TEtime);
  CompTime[myid][0] = TEtime - TStime;
  Output_CompTime();
  for (i=0; i<numprocs; i++){
    free(CompTime[i]);
  }
  free(CompTime);

  /* merge log files */
  Merge_LogFile(argv[1]);

  /* free arrays for NEGF */

  if (Solver==4){

    TRAN_Deallocate_Atoms();
    TRAN_Deallocate_RestartFile("left");
    TRAN_Deallocate_RestartFile("right");
  }

  /* free arrays */

  Free_Arrays(0);

  /* print memory */

  PrintMemory("total",0,"sum");

  MPI_Barrier(MPI_COMM_WORLD1);
  if (myid==Host_ID){
    printf("\nThe calculation was normally finished.\n");fflush(stdout);
  }

  /* if OpenMX is called by MPI_spawn. */

  if (MPI_spawn_flag==1){

    MPI_Comm_get_parent(&mpi_comm_parent);

    if(mpi_comm_parent!=MPI_COMM_NULL){
      //MPI_Comm_disconnect(&mpi_comm_parent);
    }

    fclose(MPI_spawn_stream);

    MPI_Finalize();
  }   
  else{
    MPI_Finalize();
    exit(0);
  }

  return 0;
}
