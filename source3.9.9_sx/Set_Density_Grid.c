/**********************************************************************
  Set_Density_Grid.c:

     Set_Density_Grid.c is a subroutine to calculate a charge density 
     on grid by one-particle wave functions.

  Log of Set_Density_Grid.c:

     22/Nov/2001  Released by T. Ozaki
     19/Apr/2013  Modified by A.M. Ito     

***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "openmx_common.h"
#include "mpi.h"
#include <omp.h>

#define  measure_time   0

#define _LONG_VEC2
#define _LONG_VEC3
#define AITUNE_T2DEN

double Set_Density_Grid(int Cnt_kind, int Calc_CntOrbital_ON, double *****CDM, double **Density_Grid_B0)
{
  static int firsttime=1;
  int al,L0,Mul0,M0,p,size1,size2;
  int Gc_AN,Mc_AN,Mh_AN,LN,AN,BN,CN;
  int n1,n2,n3,k1,k2,k3,N3[4];
  int Cwan,NO0,NO1,Rn,N,Hwan,i,j,k,n;
  int NN_S,NN_R;
  unsigned long long int N2D,n2D,GN; 
  int Max_Size,My_Max;
  int size_Tmp_Den_Grid;
  int size_Den_Snd_Grid_A2B;
  int size_Den_Rcv_Grid_A2B;
  int h_AN,Gh_AN,Rnh,spin,Nc,GRc,Nh,Nog;
  int Nc_0,Nc_1,Nc_2,Nc_3,Nh_0,Nh_1,Nh_2,Nh_3;

  double threshold;
  double tmp0,tmp1,sk1,sk2,sk3,tot_den,sum;
  double tmp0_0,tmp0_1,tmp0_2,tmp0_3;
  double sum_0,sum_1,sum_2,sum_3;
  double d1,d2,d3,cop,sip,sit,cot;
  double x,y,z,Cxyz[4];
  double TStime,TEtime;
  double ***Tmp_Den_Grid;
  double **Den_Snd_Grid_A2B;
  double **Den_Rcv_Grid_A2B;
  double *tmp_array;
  double *tmp_array2;
  
  int *Snd_Size,*Rcv_Size;
  int numprocs,myid,tag=999,ID,IDS,IDR;
  double Stime_atom, Etime_atom;
  double time0,time1,time2;

  MPI_Status stat;
  MPI_Request request;
  MPI_Status *stat_send;
  MPI_Status *stat_recv;
  MPI_Request *request_send;
  MPI_Request *request_recv;

  /* for OpenMP */
  int OMPID,Nthrds;

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);
  
  dtime(&TStime);

  /* allocation of arrays */

  size_Tmp_Den_Grid = 0;
  Tmp_Den_Grid = (double***)malloc(sizeof(double**)*(SpinP_switch+1)); 
  for (i=0; i<(SpinP_switch+1); i++){
    Tmp_Den_Grid[i] = (double**)malloc(sizeof(double*)*(Matomnum+1)); 
    Tmp_Den_Grid[i][0] = (double*)malloc(sizeof(double)*1); 
    for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
      Gc_AN = F_M2G[Mc_AN];
      Tmp_Den_Grid[i][Mc_AN] = (double*)malloc(sizeof(double)*GridN_Atom[Gc_AN]);
	  
      /* AITUNE */
      for (Nc=0; Nc<GridN_Atom[Gc_AN]; Nc++){
	Tmp_Den_Grid[i][Mc_AN][Nc] = 0.0;
      }
      
      size_Tmp_Den_Grid += GridN_Atom[Gc_AN];
    }
  }

  size_Den_Snd_Grid_A2B = 0; 
  Den_Snd_Grid_A2B = (double**)malloc(sizeof(double*)*numprocs);
  for (ID=0; ID<numprocs; ID++){
    Den_Snd_Grid_A2B[ID] = (double*)malloc(sizeof(double)*Num_Snd_Grid_A2B[ID]*(SpinP_switch+1));
    size_Den_Snd_Grid_A2B += Num_Snd_Grid_A2B[ID]*(SpinP_switch+1);
  }  

  size_Den_Rcv_Grid_A2B = 0;   
  Den_Rcv_Grid_A2B = (double**)malloc(sizeof(double*)*numprocs);
  for (ID=0; ID<numprocs; ID++){
    Den_Rcv_Grid_A2B[ID] = (double*)malloc(sizeof(double)*Num_Rcv_Grid_A2B[ID]*(SpinP_switch+1));
    size_Den_Rcv_Grid_A2B += Num_Rcv_Grid_A2B[ID]*(SpinP_switch+1);   
  }

  /* PrintMemory */

  if (firsttime==1){
    PrintMemory("Set_Density_Grid: AtomDen_Grid",    sizeof(double)*size_Tmp_Den_Grid, NULL);
    PrintMemory("Set_Density_Grid: Den_Snd_Grid_A2B",sizeof(double)*size_Den_Snd_Grid_A2B, NULL);
    PrintMemory("Set_Density_Grid: Den_Rcv_Grid_A2B",sizeof(double)*size_Den_Rcv_Grid_A2B, NULL);
    firsttime = 0;
  }

  /****************************************************
                when orbital optimization
  ****************************************************/

  if (Calc_CntOrbital_ON==1 && Cnt_kind==0 && Cnt_switch==1){
      
    for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
       
      dtime(&Stime_atom);
      
      /* COrbs_Grid */
 
      Gc_AN = M2G[Mc_AN];
      Cwan = WhatSpecies[Gc_AN];
      NO0 = Spe_Total_CNO[Cwan]; 
      for (Nc=0; Nc<GridN_Atom[Gc_AN]; Nc++){

        al = -1;
	for (L0=0; L0<=Spe_MaxL_Basis[Cwan]; L0++){
	  for (Mul0=0; Mul0<Spe_Num_CBasis[Cwan][L0]; Mul0++){
	    for (M0=0; M0<=2*L0; M0++){

	      al++;
	      tmp0 = 0.0;

	      for (p=0; p<Spe_Specified_Num[Cwan][al]; p++){
	        j = Spe_Trans_Orbital[Cwan][al][p];  
	        tmp0 += CntCoes[Mc_AN][al][p]*Orbs_Grid[Mc_AN][Nc][j];/* AITUNE */
	      }

	      COrbs_Grid[Mc_AN][al][Nc] = (Type_Orbs_Grid)tmp0;
	    }
	  }
        }
      }

      dtime(&Etime_atom);
      time_per_atom[Gc_AN] += Etime_atom - Stime_atom;
    }

    /**********************************************
     MPI:

     COrbs_Grid    
    ***********************************************/

    /* allocation of arrays  */
    Snd_Size = (int*)malloc(sizeof(int)*numprocs); 
    Rcv_Size = (int*)malloc(sizeof(int)*numprocs); 

    /* find data size for sending and receiving */

    My_Max = -10000;
    for (ID=0; ID<numprocs; ID++){

      IDS = (myid + ID) % numprocs;
      IDR = (myid - ID + numprocs) % numprocs;

      if (ID!=0){
        /*  sending size */
        if (F_Snd_Num[IDS]!=0){
          /* find data size */ 
          size1 = 0; 
          for (n=0; n<F_Snd_Num[IDS]; n++){
            Gc_AN = Snd_GAN[IDS][n];
            Cwan = WhatSpecies[Gc_AN];
            size1 += GridN_Atom[Gc_AN]*Spe_Total_CNO[Cwan];
          }

          Snd_Size[IDS] = size1;
          MPI_Isend(&size1, 1, MPI_INT, IDS, tag, mpi_comm_level1, &request);
        }
        else{
          Snd_Size[IDS] = 0;
        }

        /* receiving size */
        if (F_Rcv_Num[IDR]!=0){
          MPI_Recv(&size2, 1, MPI_INT, IDR, tag, mpi_comm_level1, &stat);
          Rcv_Size[IDR] = size2;
        }
        else{
          Rcv_Size[IDR] = 0;
        }
        if (F_Snd_Num[IDS]!=0) MPI_Wait(&request,&stat);
      } 
      else{
        Snd_Size[IDS] = 0;
        Rcv_Size[IDR] = 0;
      }

      if (My_Max<Snd_Size[IDS]) My_Max = Snd_Size[IDS];
      if (My_Max<Rcv_Size[IDR]) My_Max = Rcv_Size[IDR];

    }  

    MPI_Allreduce(&My_Max, &Max_Size, 1, MPI_INT, MPI_MAX, mpi_comm_level1);
    /* allocation of arrays */ 
    tmp_array  = (double*)malloc(sizeof(double)*Max_Size);
    tmp_array2 = (double*)malloc(sizeof(double)*Max_Size);

    /* send and receive COrbs_Grid */

    for (ID=0; ID<numprocs; ID++){

      IDS = (myid + ID) % numprocs;
      IDR = (myid - ID + numprocs) % numprocs;

      if (ID!=0){

        /* sending of data */ 

        if (F_Snd_Num[IDS]!=0){

          /* find data size */
          size1 = Snd_Size[IDS];

          /* multidimentional array to vector array */
          k = 0; 
          for (n=0; n<F_Snd_Num[IDS]; n++){
            Mc_AN = Snd_MAN[IDS][n];
            Gc_AN = Snd_GAN[IDS][n];
            Cwan = WhatSpecies[Gc_AN];
            NO0 = Spe_Total_CNO[Cwan]; 
            for (i=0; i<NO0; i++){
              for (Nc=0; Nc<GridN_Atom[Gc_AN]; Nc++){
                tmp_array[k] = COrbs_Grid[Mc_AN][i][Nc];
                k++;
              }          
            }
          } 

          /* MPI_Isend */
          MPI_Isend(&tmp_array[0], size1, MPI_DOUBLE, IDS,
                    tag, mpi_comm_level1, &request);
        }

        /* receiving of block data */

        if (F_Rcv_Num[IDR]!=0){

          /* find data size */
          size2 = Rcv_Size[IDR]; 

          /* MPI_Recv */
          MPI_Recv(&tmp_array2[0], size2, MPI_DOUBLE, IDR, tag, mpi_comm_level1, &stat);

          k = 0;
          Mc_AN = F_TopMAN[IDR] - 1;
          for (n=0; n<F_Rcv_Num[IDR]; n++){
            Mc_AN++;
            Gc_AN = Rcv_GAN[IDR][n];
            Cwan = WhatSpecies[Gc_AN];
            NO0 = Spe_Total_CNO[Cwan]; 

            for (i=0; i<NO0; i++){
              for (Nc=0; Nc<GridN_Atom[Gc_AN]; Nc++){
                COrbs_Grid[Mc_AN][i][Nc] = tmp_array2[k];
                k++;
              }          
            }
          }
        }
        if (F_Snd_Num[IDS]!=0) MPI_Wait(&request,&stat);
      } 
    }  

    /* freeing of arrays  */
    free(tmp_array);
    free(tmp_array2);
    free(Snd_Size);
    free(Rcv_Size);
  }

  /**********************************************
              calculate Tmp_Den_Grid
  ***********************************************/
#if measure_time
  dtime(&time1);
#endif


  /* AITUNE ========================== */ 
  int OneD_Nloop = 0;
  int ai_MaxNc = 0;
  for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
    int Gc_AN = M2G[Mc_AN];    
    for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
      OneD_Nloop++;
      if(ai_MaxNc < GridN_Atom[Gc_AN]) {ai_MaxNc = GridN_Atom[Gc_AN];}
    }
  }  
  /* ai_MaxNc is maximum of GridN_Atom[] */

  int gNthrds;
#pragma omp parallel
  {
    gNthrds = omp_get_num_threads();
  }

  
  /* ========================== AITUNE */ 

#pragma omp parallel shared(myid,G2ID,Orbs_Grid_FNAN,List_YOUSO,time_per_atom,Tmp_Den_Grid,Orbs_Grid,COrbs_Grid,Cnt_switch,Cnt_kind,GListTAtoms2,GListTAtoms1,NumOLG,CDM,SpinP_switch,WhatSpecies,ncn,F_G2M,natn,Spe_Total_CNO,M2G) private(OMPID,Nthrds,Mc_AN,h_AN,Stime_atom,Etime_atom,Gc_AN,Cwan,NO0,Gh_AN,Mh_AN,Rnh,Hwan,NO1,spin,i,j,Nog,Nc_0,Nc_1,Nc_2,Nc_3,Nh_0,Nh_1,Nh_2,Nh_3,sum_0,sum_1,sum_2,sum_3,tmp0_0,tmp0_1,tmp0_2,tmp0_3,Nc,Nh,sum,tmp0)
  {

    
    /* get info. on OpenMP */ 

    OMPID = omp_get_thread_num();
    Nthrds = omp_get_num_threads();



    #pragma omp for
    for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){

      dtime(&Stime_atom);


      /* set data on Mc_AN */

      const int Gc_AN = M2G[Mc_AN];
      const int Cwan = WhatSpecies[Gc_AN];
      const int NO0 = Spe_Total_CNO[Cwan]; 
	  
      
      for (int spin=0; spin<=SpinP_switch; spin++){	
        for (int Nc=0; Nc<GridN_Atom[Gc_AN]; Nc++){
          Tmp_Den_Grid[spin][Mc_AN][Nc] = 0.0;
        }
      }





#ifdef AITUNE_T2DEN
      int max_NumOLG = 0;
      for (int h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
        if(max_NumOLG < NumOLG[Mc_AN][h_AN]) max_NumOLG = NumOLG[Mc_AN][h_AN];
      }

      double** T2_Den_Grid = (double**)malloc(sizeof(double*)*(SpinP_switch+1));
      {
        double* pnt = (double*)malloc(sizeof(double)*(SpinP_switch+1)*max_NumOLG);
        for (int spin=0; spin<=SpinP_switch; spin++){ 
          T2_Den_Grid[spin] = pnt + max_NumOLG * spin;
        }
        
      }
#endif



      for (int h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){

	/* set data on h_AN */
    
	const int Gh_AN = natn[Gc_AN][h_AN];
	const int Mh_AN = F_G2M[Gh_AN];
	const int Rnh = ncn[Gc_AN][h_AN];
	const int Hwan = WhatSpecies[Gh_AN];
	const int NO1 = Spe_Total_CNO[Hwan];


#ifdef AITUNE_T2DEN
  {
    double* pnt = T2_Den_Grid[0];
    for (int Nog=0; Nog<(SpinP_switch+1)*max_NumOLG; Nog++){
      pnt[Nog] = 0.0;
    }
  }
#endif

#ifdef _LONG_VEC3
#define Type_Orbs_Grid_T   double   //Type_Orbs_Grid
  
  Type_Orbs_Grid_T** Orbs_Grid_FNAN_t3 = (Type_Orbs_Grid_T**)malloc(sizeof(Type_Orbs_Grid_T*)*NO1); 
  Type_Orbs_Grid_T* pnt = (Type_Orbs_Grid_T*)malloc(sizeof(Type_Orbs_Grid_T)*(NumOLG[Mc_AN][h_AN])*NO1); 
  for (int j=0; j<NO1; j++){
    Orbs_Grid_FNAN_t3[j] = pnt + NumOLG[Mc_AN][h_AN] * j;
	  if (G2ID[Gh_AN]!=myid){
	    for (int Nog=0; Nog<NumOLG[Mc_AN][h_AN]; Nog++){
	      Orbs_Grid_FNAN_t3[j][Nog] = Orbs_Grid_FNAN[Mc_AN][h_AN][Nog][j];
	    }
    }else{
      for (int Nog=0; Nog<NumOLG[Mc_AN][h_AN]; Nog++){
        const int Nh = GListTAtoms2[Mc_AN][h_AN][Nog]; 
	      Orbs_Grid_FNAN_t3[j][Nog] = Orbs_Grid[Mh_AN][Nh][j];
	    }
    }
  }

  
#elif defined _LONG_VEC2
#define Type_Orbs_Grid_T   double   //Type_Orbs_Grid

  Type_Orbs_Grid_T** Orbs_Grid_FNAN_t2 = (Type_Orbs_Grid_T**)malloc(sizeof(Type_Orbs_Grid_T*)*NO1); 
  Type_Orbs_Grid_T* pnt = (Type_Orbs_Grid_T*)malloc(sizeof(Type_Orbs_Grid_T)*(NumOLG[Mc_AN][h_AN])*NO1); 
  for (int j=0; j<NO1; j++){
    Orbs_Grid_FNAN_t2[j] = pnt + NumOLG[Mc_AN][h_AN] * j;
	  if (G2ID[Gh_AN]!=myid){
	    for (int Nog=0; Nog<NumOLG[Mc_AN][h_AN]; Nog++){
	      Orbs_Grid_FNAN_t2[j][Nog] = Orbs_Grid_FNAN[Mc_AN][h_AN][Nog][j];
	    }
    }else{
      for (int Nog=0; Nog<NumOLG[Mc_AN][h_AN]; Nog++){
        const int Nh = GListTAtoms2[Mc_AN][h_AN][Nog]; 
	      Orbs_Grid_FNAN_t2[j][Nog] = Orbs_Grid[Mh_AN][Nh][j];
	    }
    }
  }


  Type_Orbs_Grid_T** Orbs_Grid_SELF_t2 = (Type_Orbs_Grid_T**)malloc(sizeof(Type_Orbs_Grid_T*)*NO0); 
  
  Type_Orbs_Grid_T* pnt0 = (Type_Orbs_Grid_T*)malloc(sizeof(Type_Orbs_Grid_T)*NumOLG[Mc_AN][h_AN]*NO0);
  for (int i=0; i<NO0; i++){
    Orbs_Grid_SELF_t2[i] = pnt0 + NumOLG[Mc_AN][h_AN] * i;
    /*
    for (int Nog=0; Nog<NumOLG[Mc_AN][h_AN]; Nog++){
      const int Nc = GListTAtoms1[Mc_AN][h_AN][Nog];
      Orbs_Grid_SELF_t2[i][Nog] = Orbs_Grid[Mc_AN][Nc][i];
    }
    */
  }
  
#endif

  /*******************************************************
     Note:
     Important things for turning, 
     Nc, the elements in GListTAtoms1 are not duplicated.
     Nh, the elements in GListTAtoms2 are not duplicated.  
     Therefore, memory access by Nc and Nh is independent.
     We can use pragma ivdep
   
   * *******************************************************/
	
	if (Cnt_kind==0 && Cnt_switch==1){
    
#if 1
  if(NO1 <= 30){
    for (int spin=0; spin<=SpinP_switch; spin++){ 

      for (int i=0; i<NO0; i++){

#pragma ivdep
#pragma _NEC vector
#pragma _NEC loop_count(1000)
#pragma _NEC ivdep
	        for (Nog=0; Nog<NumOLG[Mc_AN][h_AN]; Nog++){//average vector length is about 1000 //
            const int Nc = GListTAtoms1[Mc_AN][h_AN][Nog];
            const int Nh = GListTAtoms2[Mc_AN][h_AN][Nog]; 
      
            const double orbs_i = COrbs_Grid[Mc_AN][i][Nc];
            double sum = 0.0;

#pragma _NEC nointerchange
#pragma _NEC unroll(30)
#pragma _NEC ivdep
	        for (int j=0; j<30; j++){
            if(j<NO1){
              const double orbs_j = COrbs_Grid[Mh_AN][j][Nh];
              sum += orbs_j*CDM[spin][Mc_AN][h_AN][i][j];
            }
          }
          
#ifdef AITUNE_T2DEN
            T2_Den_Grid[spin][Nog] += orbs_i*sum;
#else
            Tmp_Den_Grid[spin][Mc_AN][Nc] += orbs_i*sum;
#endif
          
          
	      }
      }
    }	

  }else{
    for (int spin=0; spin<=SpinP_switch; spin++){ 

      for (int i=0; i<NO0; i++){
	      for (int j=0; j<NO1; j++){

#pragma ivdep
#pragma _NEC vector
#pragma _NEC loop_count(1000)
#pragma _NEC ivdep
	        for (Nog=0; Nog<NumOLG[Mc_AN][h_AN]; Nog++){//average vector length is about 1000 //
            const int Nc = GListTAtoms1[Mc_AN][h_AN][Nog];
            const int Nh = GListTAtoms2[Mc_AN][h_AN][Nog]; 
      
            const double orbs_i = COrbs_Grid[Mc_AN][i][Nc];
            const double orbs_j = COrbs_Grid[Mh_AN][j][Nh];
#ifdef AITUNE_T2DEN
            T2_Den_Grid[spin][Nog] += orbs_i*orbs_j*CDM[spin][Mc_AN][h_AN][i][j];
#else
            Tmp_Den_Grid[spin][Mc_AN][Nc] += orbs_i*orbs_j*CDM[spin][Mc_AN][h_AN][i][j];
#endif
            
          }
	      }
      }
    }	
  }
#else
    for (int spin=0; spin<=SpinP_switch; spin++){ 

      for (int i=0; i<NO0; i++){
	      for (int j=0; j<NO1; j++){

#pragma ivdep
#pragma _NEC vector
#pragma _NEC loop_count(1000)
#pragma _NEC ivdep
	        for (Nog=0; Nog<NumOLG[Mc_AN][h_AN]; Nog++){//average vector length is about 1000 //
            const int Nc = GListTAtoms1[Mc_AN][h_AN][Nog];
            const int Nh = GListTAtoms2[Mc_AN][h_AN][Nog]; 
      
            const double orbs_i = COrbs_Grid[Mc_AN][i][Nc];
            const double orbs_j = COrbs_Grid[Mh_AN][j][Nh];
#ifdef AITUNE_T2DEN
            T2_Den_Grid[spin][Nog] += orbs_i*orbs_j*CDM[spin][Mc_AN][h_AN][i][j];
#else
            Tmp_Den_Grid[spin][Mc_AN][Nc] += orbs_i*orbs_j*CDM[spin][Mc_AN][h_AN][i][j];
#endif
            
          }
	      }
      }
    }	
#endif      
  }
	

#ifdef _LONG_VEC3
  else {

  if(NO1 <= 30){

    for (int spin=0; spin<=SpinP_switch; spin++){ 

#pragma _NEC loop_count(30)
#pragma _NEC ivdep
      for (int i=0; i<NO0; i++){

#pragma ivdep
#pragma _NEC vector
#pragma _NEC loop_count(1000)
#pragma _NEC ivdep
	        for (int Nog=0; Nog<NumOLG[Mc_AN][h_AN]; Nog++){//average vector length is about 1000 //
	          const int Nc = GListTAtoms1[Mc_AN][h_AN][Nog];
            const double orbs_i = Orbs_Grid[Mc_AN][Nc][i];
	          
            //const int Nh = GListTAtoms2[Mc_AN][h_AN][Nog]; 
            //const double orbs_i = Orbs_Grid_SELF_t2[i][Nog];

            double sum = 0.0;

#pragma _NEC nointerchange            
#pragma _NEC unroll(30)
#pragma _NEC ivdep
	        for (int j=0; j<30; j++){            
            if(j < NO1){
              const double orbs_j = Orbs_Grid_FNAN_t3[j][Nog];    //order of [Nog][j] is slow, probably cash load as one line//
	            //const double orbs_j = Orbs_Grid_FNAN_t3[Nog][j];
              //const double orbs_j = Orbs_Grid[Mh_AN][Nh][j];            
              sum += orbs_j*CDM[spin][Mc_AN][h_AN][i][j];
            }
          }

#ifdef AITUNE_T2DEN
          T2_Den_Grid[spin][Nog] += orbs_i*sum;
#else
          //const int Nc = GListTAtoms1[Mc_AN][h_AN][Nog];
          Tmp_Den_Grid[spin][Mc_AN][Nc] += orbs_i*sum;
#endif
	      }
	    }
    }
  }else{
    for (int spin=0; spin<=SpinP_switch; spin++){ 

#pragma _NEC loop_count(30)
#pragma _NEC ivdep
      for (int i=0; i<NO0; i++){
#pragma _NEC loop_count(30)
#pragma _NEC ivdep
	      for (int j=0; j<NO1; j++){            
#pragma omp simd
#pragma ivdep
#pragma _NEC vector
#pragma _NEC loop_count(1000)
#pragma _NEC ivdep
	        for (int Nog=0; Nog<NumOLG[Mc_AN][h_AN]; Nog++){//average vector length is about 1000 //
	          //const double orbs_i = Orbs_Grid_t2[i][Nc];
            const int Nc = GListTAtoms1[Mc_AN][h_AN][Nog];
            const double orbs_i = Orbs_Grid[Mc_AN][Nc][i];
            const double orbs_j = Orbs_Grid_FNAN_t3[j][Nog];  
	          //const double orbs_j = Orbs_Grid_FNAN_t3[Nog][j];
	          

            //const double orbs_i = Orbs_Grid_SELF_t2[i][Nog];
	          //const double orbs_j = Orbs_Grid_FNAN_t2[j][Nog];
#ifdef AITUNE_T2DEN
            T2_Den_Grid[spin][Nog] += orbs_i*orbs_j*CDM[spin][Mc_AN][h_AN][i][j];
#else
            //const int Nc = GListTAtoms1[Mc_AN][h_AN][Nog];
            Tmp_Den_Grid[spin][Mc_AN][Nc] += orbs_i*orbs_j*CDM[spin][Mc_AN][h_AN][i][j];
#endif
          }
	      }
	    }
    }
  }
  }

  
  free(Orbs_Grid_FNAN_t3[0]);
  free(Orbs_Grid_FNAN_t3);
  

#elif defined _LONG_VEC2  //order before ver 3.6//
  else{

    for (int spin=0; spin<=SpinP_switch; spin++){ 

#pragma _NEC loop_count(30)
#pragma _NEC ivdep
      for (int i=0; i<NO0; i++){
#pragma _NEC loop_count(30)
#pragma _NEC ivdep
	      for (int j=0; j<NO1; j++){            
#pragma omp simd
#pragma ivdep
#pragma _NEC vector
#pragma _NEC loop_count(1000)
#pragma _NEC ivdep
	        for (int Nog=0; Nog<NumOLG[Mc_AN][h_AN]; Nog++){//average vector length is about 1000 //
	          //const double orbs_i = Orbs_Grid_t2[i][Nc];

            const double orbs_i = Orbs_Grid_SELF_t2[i][Nog];
	          const double orbs_j = Orbs_Grid_FNAN_t2[j][Nog];
#ifdef AITUNE_T2DEN
            T2_Den_Grid[spin][Nog] += orbs_i*orbs_j*CDM[spin][Mc_AN][h_AN][i][j];
#else
            const int Nc = GListTAtoms1[Mc_AN][h_AN][Nog];
            Tmp_Den_Grid[spin][Mc_AN][Nc] += orbs_i*orbs_j*CDM[spin][Mc_AN][h_AN][i][j];
#endif
          }
	      }
	    }
    }
  }


  free(Orbs_Grid_FNAN_t2[0]);
  free(Orbs_Grid_FNAN_t2);
  free(Orbs_Grid_SELF_t2[0]);
  free(Orbs_Grid_SELF_t2);
  
#else
  else if (G2ID[Gh_AN]==myid){
    for (int spin=0; spin<=SpinP_switch; spin++){ 

#pragma _NEC loop_count(30)
#pragma _NEC ivdep
      for (int i=0; i<NO0; i++){
#pragma omp simd
#pragma ivdep
#pragma _NEC vector
#pragma _NEC loop_count(1000)
#pragma _NEC ivdep
	        for (int Nog=0; Nog<NumOLG[Mc_AN][h_AN]; Nog++){//average vector length is about 1000 //
	          const int Nc = GListTAtoms1[Mc_AN][h_AN][Nog];
            const int Nh = GListTAtoms2[Mc_AN][h_AN][Nog]; 
#pragma _NEC loop_count(30)
#pragma _NEC ivdep
	      for (int j=0; j<NO1; j++){            
            const double orbs_i = Orbs_Grid[Mc_AN][Nc][i];
	          const double orbs_j = Orbs_Grid[Mh_AN][Nh][j];

            Tmp_Den_Grid[spin][Mc_AN][Nc] += orbs_i*orbs_j*CDM[spin][Mc_AN][h_AN][i][j];
          }
	      }
	    }
    }
  }
  else{

#pragma _NEC loop_count(2)
    for (int spin=0; spin<=SpinP_switch; spin++){
#pragma _NEC loop_count(30)
#pragma _NEC ivdep
      for (int i=0; i<NO0; i++){
#pragma omp simd
#pragma ivdep
#pragma _NEC vector
#pragma _NEC loop_count(1000)
#pragma _NEC ivdep
	        for (int Nog=0; Nog<NumOLG[Mc_AN][h_AN]; Nog++){  //average vector length is about 1000 //
#pragma _NEC loop_count(30)
#pragma _NEC ivdep
	      for (int j=0; j<NO1; j++){
            const int Nc = GListTAtoms1[Mc_AN][h_AN][Nog];
            const int Nh = GListTAtoms2[Mc_AN][h_AN][Nog]; 
	          
            const double orbs_i = Orbs_Grid[Mc_AN][Nc][i];
	          const double orbs_j = Orbs_Grid_FNAN[Mc_AN][h_AN][Nog][j];

            Tmp_Den_Grid[spin][Mc_AN][Nc] += orbs_i*orbs_j*CDM[spin][Mc_AN][h_AN][i][j];
          }
	      }
	    }
    }
  }

#endif    
	

    
	

#ifdef AITUNE_T2DEN
#pragma _NEC loop_count(2)
  for (int spin=0; spin<=SpinP_switch; spin++){
#pragma _NEC vector
#pragma _NEC loop_count(1000)
#pragma _NEC ivdep
    for (int Nog=0; Nog<NumOLG[Mc_AN][h_AN]; Nog++){  //average vector length is about 1000 //
      const int Nc = GListTAtoms1[Mc_AN][h_AN][Nog];
      Tmp_Den_Grid[spin][Mc_AN][Nc] += T2_Den_Grid[spin][Nog];      
    }
  }
#endif

      } /* h_AN */


#ifdef AITUNE_T2DEN
  free(T2_Den_Grid[0]);
  free(T2_Den_Grid);
#endif
    

      dtime(&Etime_atom);
      time_per_atom[Gc_AN] += Etime_atom - Stime_atom;

    } /* Mc_AN */

    /* freeing of arrays */ 



#pragma omp flush(Tmp_Den_Grid)

  } /* #pragma omp parallel */


#if measure_time
  dtime(&time2);
  if(myid==0 && measure_time){
    printf("Time for Part1=%18.5f\n",(time2-time1));fflush(stdout);
  }
#endif

  /******************************************************
      MPI communication from the partitions A to B 
  ******************************************************/

  /* copy Tmp_Den_Grid to Den_Snd_Grid_A2B */

  for (ID=0; ID<numprocs; ID++) Num_Snd_Grid_A2B[ID] = 0;
  
  N2D = Ngrid1*Ngrid2;

  for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){

    Gc_AN = M2G[Mc_AN];
#ifndef _NEC
    for (AN=0; AN<GridN_Atom[Gc_AN]; AN++){

      GN = GridListAtom[Mc_AN][AN];
      GN2N(GN,N3);
      n2D = N3[1]*Ngrid2 + N3[2];
      ID = (int)(n2D*(unsigned long long int)numprocs/N2D);

      if (SpinP_switch==0){
        Den_Snd_Grid_A2B[ID][Num_Snd_Grid_A2B[ID]] = Tmp_Den_Grid[0][Mc_AN][AN];
      }
      else if (SpinP_switch==1){
        Den_Snd_Grid_A2B[ID][Num_Snd_Grid_A2B[ID]*2+0] = Tmp_Den_Grid[0][Mc_AN][AN];
        Den_Snd_Grid_A2B[ID][Num_Snd_Grid_A2B[ID]*2+1] = Tmp_Den_Grid[1][Mc_AN][AN];
      }
      else if (SpinP_switch==3){
        Den_Snd_Grid_A2B[ID][Num_Snd_Grid_A2B[ID]*4+0] = Tmp_Den_Grid[0][Mc_AN][AN];
        Den_Snd_Grid_A2B[ID][Num_Snd_Grid_A2B[ID]*4+1] = Tmp_Den_Grid[1][Mc_AN][AN];
        Den_Snd_Grid_A2B[ID][Num_Snd_Grid_A2B[ID]*4+2] = Tmp_Den_Grid[2][Mc_AN][AN];
        Den_Snd_Grid_A2B[ID][Num_Snd_Grid_A2B[ID]*4+3] = Tmp_Den_Grid[3][Mc_AN][AN];
      }

      Num_Snd_Grid_A2B[ID]++;
    }
#else  //_NEC
      if ( 0 == SpinP_switch){
    for (AN=0; AN<GridN_Atom[Gc_AN]; AN++){

      GN = GridListAtom[Mc_AN][AN];
      GN2N(GN,N3);
      n2D = N3[1]*Ngrid2 + N3[2];
      ID = (int)(n2D*(unsigned long long int)numprocs/N2D);

        Den_Snd_Grid_A2B[ID][Num_Snd_Grid_A2B[ID]] = Tmp_Den_Grid[0][Mc_AN][AN];

      Num_Snd_Grid_A2B[ID]++;
    }
	}
      else if (1 == SpinP_switch){

    for (AN=0; AN<GridN_Atom[Gc_AN]; AN++){

      GN = GridListAtom[Mc_AN][AN];
      GN2N(GN,N3);
      n2D = N3[1]*Ngrid2 + N3[2];
      ID = (int)(n2D*(unsigned long long int)numprocs/N2D);

        Den_Snd_Grid_A2B[ID][Num_Snd_Grid_A2B[ID]*2+0] = Tmp_Den_Grid[0][Mc_AN][AN];
        Den_Snd_Grid_A2B[ID][Num_Snd_Grid_A2B[ID]*2+1] = Tmp_Den_Grid[1][Mc_AN][AN];

      Num_Snd_Grid_A2B[ID]++;
    }
      }
      else if (3 == SpinP_switch){
    for (AN=0; AN<GridN_Atom[Gc_AN]; AN++){

      GN = GridListAtom[Mc_AN][AN];
      GN2N(GN,N3);
      n2D = N3[1]*Ngrid2 + N3[2];
      ID = (int)(n2D*(unsigned long long int)numprocs/N2D);


        Den_Snd_Grid_A2B[ID][Num_Snd_Grid_A2B[ID]*4+0] = Tmp_Den_Grid[0][Mc_AN][AN];
        Den_Snd_Grid_A2B[ID][Num_Snd_Grid_A2B[ID]*4+1] = Tmp_Den_Grid[1][Mc_AN][AN];
        Den_Snd_Grid_A2B[ID][Num_Snd_Grid_A2B[ID]*4+2] = Tmp_Den_Grid[2][Mc_AN][AN];
        Den_Snd_Grid_A2B[ID][Num_Snd_Grid_A2B[ID]*4+3] = Tmp_Den_Grid[3][Mc_AN][AN];


      Num_Snd_Grid_A2B[ID]++;
    }
      }
#endif //_NEC
  }    

  /* MPI: A to B */  

  request_send = malloc(sizeof(MPI_Request)*NN_A2B_S);
  request_recv = malloc(sizeof(MPI_Request)*NN_A2B_R);
  stat_send = malloc(sizeof(MPI_Status)*NN_A2B_S);
  stat_recv = malloc(sizeof(MPI_Status)*NN_A2B_R);

  NN_S = 0;
  NN_R = 0;

  tag = 999;
  for (ID=1; ID<numprocs; ID++){

    IDS = (myid + ID) % numprocs;
    IDR = (myid - ID + numprocs) % numprocs;

    if (Num_Snd_Grid_A2B[IDS]!=0){
      MPI_Isend( &Den_Snd_Grid_A2B[IDS][0], Num_Snd_Grid_A2B[IDS]*(SpinP_switch+1), 
	         MPI_DOUBLE, IDS, tag, mpi_comm_level1, &request_send[NN_S]);
      NN_S++;
    }

    if (Num_Rcv_Grid_A2B[IDR]!=0){
      MPI_Irecv( &Den_Rcv_Grid_A2B[IDR][0], Num_Rcv_Grid_A2B[IDR]*(SpinP_switch+1), 
  	         MPI_DOUBLE, IDR, tag, mpi_comm_level1, &request_recv[NN_R]);
      NN_R++;
    }
  }

  if (NN_S!=0) MPI_Waitall(NN_S,request_send,stat_send);
  if (NN_R!=0) MPI_Waitall(NN_R,request_recv,stat_recv);

  free(request_send);
  free(request_recv);
  free(stat_send);
  free(stat_recv);

  /* for myid */
  for (i=0; i<Num_Rcv_Grid_A2B[myid]*(SpinP_switch+1); i++){
    Den_Rcv_Grid_A2B[myid][i] = Den_Snd_Grid_A2B[myid][i];
  }

  /******************************************************
   superposition of rho_i to calculate charge density 
   in the partition B.
  ******************************************************/

  /* initialize arrays */
  
  for (spin=0; spin<(SpinP_switch+1); spin++){
    for (BN=0; BN<My_NumGridB_AB; BN++){
      Density_Grid_B0[spin][BN] = 0.0;
    }
  }
  
  /* superposition of densities rho_i */

  for (ID=0; ID<numprocs; ID++){

    for (LN=0; LN<Num_Rcv_Grid_A2B[ID]; LN++){

      BN    = Index_Rcv_Grid_A2B[ID][3*LN+0];      
      Gc_AN = Index_Rcv_Grid_A2B[ID][3*LN+1];        
      GRc   = Index_Rcv_Grid_A2B[ID][3*LN+2]; 

      if (Solver!=4 || (Solver==4 && atv_ijk[GRc][1]==0 )){

	/* spin collinear non-polarization */
	if ( SpinP_switch==0 ){
	  Density_Grid_B0[0][BN] += Den_Rcv_Grid_A2B[ID][LN];
	}

	/* spin collinear polarization */
	else if ( SpinP_switch==1 ){
	  Density_Grid_B0[0][BN] += Den_Rcv_Grid_A2B[ID][LN*2  ];
	  Density_Grid_B0[1][BN] += Den_Rcv_Grid_A2B[ID][LN*2+1];
	} 

	/* spin non-collinear */
	else if ( SpinP_switch==3 ){
	  Density_Grid_B0[0][BN] += Den_Rcv_Grid_A2B[ID][LN*4  ];
	  Density_Grid_B0[1][BN] += Den_Rcv_Grid_A2B[ID][LN*4+1];
	  Density_Grid_B0[2][BN] += Den_Rcv_Grid_A2B[ID][LN*4+2];
	  Density_Grid_B0[3][BN] += Den_Rcv_Grid_A2B[ID][LN*4+3];
	} 

      } /* if (Solve!=4.....) */           

    } /* AN */ 
  } /* ID */  

  /****************************************************
   Conjugate complex of Density_Grid[3][MN] due to
   difference in the definition between density matrix
   and charge density
  ****************************************************/

  if (SpinP_switch==3){

    for (BN=0; BN<My_NumGridB_AB; BN++){
      Density_Grid_B0[3][BN] = -Density_Grid_B0[3][BN]; 
    }
  }

  /******************************************************
             MPI: from the partitions B to D
  ******************************************************/

  Density_Grid_Copy_B2D(Density_Grid_B0);

  /* freeing of arrays */

  for (i=0; i<(SpinP_switch+1); i++){
    for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){
      free(Tmp_Den_Grid[i][Mc_AN]);
    }
    free(Tmp_Den_Grid[i]);
  }
  free(Tmp_Den_Grid);

  for (ID=0; ID<numprocs; ID++){
    free(Den_Snd_Grid_A2B[ID]);
  }  
  free(Den_Snd_Grid_A2B);

  for (ID=0; ID<numprocs; ID++){
    free(Den_Rcv_Grid_A2B[ID]);
  }
  free(Den_Rcv_Grid_A2B);

  /* elapsed time */
  dtime(&TEtime);
  time0 = TEtime - TStime;
  if(myid==0 && measure_time) printf("time0=%18.5f\n",time0);

  return time0;
}



void Data_Grid_Copy_B2C_2(double **data_B, double **data_C)
{
  static int firsttime=1;
  int CN,BN,LN,spin,i,gp,NN_S,NN_R;
  double *Work_Array_Snd_Grid_B2C;
  double *Work_Array_Rcv_Grid_B2C;
  int numprocs,myid,tag=999,ID,IDS,IDR;
  MPI_Status stat;
  MPI_Request request;
  MPI_Status *stat_send;
  MPI_Status *stat_recv;
  MPI_Request *request_send;
  MPI_Request *request_recv;

  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  /* allocation of arrays */
  
  Work_Array_Snd_Grid_B2C = (double*)malloc(sizeof(double)*GP_B2C_S[NN_B2C_S]*(SpinP_switch+1)); 
  Work_Array_Rcv_Grid_B2C = (double*)malloc(sizeof(double)*GP_B2C_R[NN_B2C_R]*(SpinP_switch+1)); 

  if (firsttime==1){
    PrintMemory("Data_Grid_Copy_B2C_2: Work_Array_Snd_Grid_B2C",
		sizeof(double)*GP_B2C_S[NN_B2C_S]*(SpinP_switch+1), NULL);
    PrintMemory("Data_Grid_Copy_B2C_2: Work_Array_Rcv_Grid_B2C",
		sizeof(double)*GP_B2C_R[NN_B2C_R]*(SpinP_switch+1), NULL);
    firsttime = 0;
  }

  /******************************************************
             MPI: from the partitions B to C
  ******************************************************/

  request_send = malloc(sizeof(MPI_Request)*NN_B2C_S);
  request_recv = malloc(sizeof(MPI_Request)*NN_B2C_R);
  stat_send = malloc(sizeof(MPI_Status)*NN_B2C_S);
  stat_recv = malloc(sizeof(MPI_Status)*NN_B2C_R);

  NN_S = 0;
  NN_R = 0;

  /* MPI_Irecv */

  for (ID=0; ID<NN_B2C_R; ID++){

    IDR = ID_NN_B2C_R[ID];
    gp = GP_B2C_R[ID];

    if (IDR!=myid){ 
      MPI_Irecv( &Work_Array_Rcv_Grid_B2C[(SpinP_switch+1)*gp], Num_Rcv_Grid_B2C[IDR]*(SpinP_switch+1),
                 MPI_DOUBLE, IDR, tag, mpi_comm_level1, &request_recv[NN_R]);
      NN_R++;
    }

  }

  /* MPI_Isend */

  for (ID=0; ID<NN_B2C_S; ID++){

    IDS = ID_NN_B2C_S[ID];
    gp = GP_B2C_S[ID];

    /* copy Density_Grid_B0 to Work_Array_Snd_Grid_B2C */

    for (LN=0; LN<Num_Snd_Grid_B2C[IDS]; LN++){
      BN = Index_Snd_Grid_B2C[IDS][LN];

      if (SpinP_switch==0){
        Work_Array_Snd_Grid_B2C[gp+LN]       = data_B[0][BN];
      }
      else if (SpinP_switch==1){
        Work_Array_Snd_Grid_B2C[2*gp+2*LN+0] = data_B[0][BN];
        Work_Array_Snd_Grid_B2C[2*gp+2*LN+1] = data_B[1][BN];
      }
      else if (SpinP_switch==3){
        Work_Array_Snd_Grid_B2C[4*gp+4*LN+0] = data_B[0][BN];
        Work_Array_Snd_Grid_B2C[4*gp+4*LN+1] = data_B[1][BN];
        Work_Array_Snd_Grid_B2C[4*gp+4*LN+2] = data_B[2][BN];
        Work_Array_Snd_Grid_B2C[4*gp+4*LN+3] = data_B[3][BN];
      }
    } /* LN */        

    if (IDS!=myid){
      MPI_Isend( &Work_Array_Snd_Grid_B2C[(SpinP_switch+1)*gp], Num_Snd_Grid_B2C[IDS]*(SpinP_switch+1), 
		 MPI_DOUBLE, IDS, tag, mpi_comm_level1, &request_send[NN_S]);
      NN_S++;
    }
  }

  /* MPI_Waitall */

  if (NN_S!=0) MPI_Waitall(NN_S,request_send,stat_send);
  if (NN_R!=0) MPI_Waitall(NN_R,request_recv,stat_recv);

  free(request_send);
  free(request_recv);
  free(stat_send);
  free(stat_recv);

  /* copy Work_Array_Rcv_Grid_B2C to data_C */

  for (ID=0; ID<NN_B2C_R; ID++){

    IDR = ID_NN_B2C_R[ID];

    if (IDR==myid){

      gp = GP_B2C_S[ID];

      for (LN=0; LN<Num_Rcv_Grid_B2C[IDR]; LN++){

	CN = Index_Rcv_Grid_B2C[IDR][LN];

	if (SpinP_switch==0){
	  data_C[0][CN] = Work_Array_Snd_Grid_B2C[gp+LN];
	}     
	else if (SpinP_switch==1){
	  data_C[0][CN] = Work_Array_Snd_Grid_B2C[2*gp+2*LN+0];
	  data_C[1][CN] = Work_Array_Snd_Grid_B2C[2*gp+2*LN+1];
	}     
	else if (SpinP_switch==3){
	  data_C[0][CN] = Work_Array_Snd_Grid_B2C[4*gp+4*LN+0];
	  data_C[1][CN] = Work_Array_Snd_Grid_B2C[4*gp+4*LN+1];
	  data_C[2][CN] = Work_Array_Snd_Grid_B2C[4*gp+4*LN+2];
	  data_C[3][CN] = Work_Array_Snd_Grid_B2C[4*gp+4*LN+3];
	}
      } /* LN */   

    }
    else {

      gp = GP_B2C_R[ID];

      for (LN=0; LN<Num_Rcv_Grid_B2C[IDR]; LN++){
	CN = Index_Rcv_Grid_B2C[IDR][LN];

	if (SpinP_switch==0){
	  data_C[0][CN] = Work_Array_Rcv_Grid_B2C[gp+LN];
	}
	else if (SpinP_switch==1){
	  data_C[0][CN] = Work_Array_Rcv_Grid_B2C[2*gp+2*LN+0];
	  data_C[1][CN] = Work_Array_Rcv_Grid_B2C[2*gp+2*LN+1];
	}     
	else if (SpinP_switch==3){
	  data_C[0][CN] = Work_Array_Rcv_Grid_B2C[4*gp+4*LN+0];
	  data_C[1][CN] = Work_Array_Rcv_Grid_B2C[4*gp+4*LN+1];
	  data_C[2][CN] = Work_Array_Rcv_Grid_B2C[4*gp+4*LN+2];
	  data_C[3][CN] = Work_Array_Rcv_Grid_B2C[4*gp+4*LN+3];
	}
      }
    }
  }

  /* if (SpinP_switch==0), 
     copy data_B[0] to data_B[1]
     copy data_C[0] to data_C[1]
  */

  if (SpinP_switch==0){
    for (BN=0; BN<My_NumGridB_AB; BN++){
      data_B[1][BN] = data_B[0][BN]; 
    }

    for (CN=0; CN<My_NumGridC; CN++){
      data_C[1][CN] = data_C[0][CN]; 
    }
  }

  /* freeing of arrays */
  free(Work_Array_Snd_Grid_B2C);
  free(Work_Array_Rcv_Grid_B2C);
}



void Data_Grid_Copy_B2C_1(double *data_B, double *data_C)
{
  static int firsttime=1;
  int CN,BN,LN,spin,i,gp,NN_S,NN_R;
  double *Work_Array_Snd_Grid_B2C;
  double *Work_Array_Rcv_Grid_B2C;
  int numprocs,myid,tag=999,ID,IDS,IDR;
  MPI_Status stat;
  MPI_Request request;
  MPI_Status *stat_send;
  MPI_Status *stat_recv;
  MPI_Request *request_send;
  MPI_Request *request_recv;

  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  /* allocation of arrays */
  
  Work_Array_Snd_Grid_B2C = (double*)malloc(sizeof(double)*GP_B2C_S[NN_B2C_S]); 
  Work_Array_Rcv_Grid_B2C = (double*)malloc(sizeof(double)*GP_B2C_R[NN_B2C_R]); 

  if (firsttime==1){
    PrintMemory("Data_Grid_Copy_B2C_1: Work_Array_Snd_Grid_B2C",
		sizeof(double)*GP_B2C_S[NN_B2C_S], NULL);
    PrintMemory("Data_Grid_Copy_B2C_1: Work_Array_Rcv_Grid_B2C",
		sizeof(double)*GP_B2C_R[NN_B2C_R], NULL);
    firsttime = 0;
  }

  /******************************************************
             MPI: from the partitions B to C
  ******************************************************/

  request_send = malloc(sizeof(MPI_Request)*NN_B2C_S);
  request_recv = malloc(sizeof(MPI_Request)*NN_B2C_R);
  stat_send = malloc(sizeof(MPI_Status)*NN_B2C_S);
  stat_recv = malloc(sizeof(MPI_Status)*NN_B2C_R);

  NN_S = 0;
  NN_R = 0;

  /* MPI_Irecv */

  for (ID=0; ID<NN_B2C_R; ID++){

    IDR = ID_NN_B2C_R[ID];
    gp = GP_B2C_R[ID];

    if (IDR!=myid){ 
      MPI_Irecv( &Work_Array_Rcv_Grid_B2C[gp], Num_Rcv_Grid_B2C[IDR],
                 MPI_DOUBLE, IDR, tag, mpi_comm_level1, &request_recv[NN_R]);
      NN_R++;
    }
  }
 
  /* MPI_Isend */

  for (ID=0; ID<NN_B2C_S; ID++){

    IDS = ID_NN_B2C_S[ID];
    gp = GP_B2C_S[ID];

    /* copy Density_Grid_B0 to Work_Array_Snd_Grid_B2C */

    for (LN=0; LN<Num_Snd_Grid_B2C[IDS]; LN++){
      BN = Index_Snd_Grid_B2C[IDS][LN];
      Work_Array_Snd_Grid_B2C[gp+LN] = data_B[BN];
    } 

    if (IDS!=myid){
      MPI_Isend( &Work_Array_Snd_Grid_B2C[gp], Num_Snd_Grid_B2C[IDS], 
		 MPI_DOUBLE, IDS, tag, mpi_comm_level1, &request_send[NN_S]);
      NN_S++;
    }
  }

  /* MPI_Waitall */

  if (NN_S!=0) MPI_Waitall(NN_S,request_send,stat_send);
  if (NN_R!=0) MPI_Waitall(NN_R,request_recv,stat_recv);

  free(request_send);
  free(request_recv);
  free(stat_send);
  free(stat_recv);

  /* copy Work_Array_Rcv_Grid_B2C to data_C */

  for (ID=0; ID<NN_B2C_R; ID++){

    IDR = ID_NN_B2C_R[ID];

    if (IDR==myid){
      gp = GP_B2C_S[ID];
      for (LN=0; LN<Num_Rcv_Grid_B2C[IDR]; LN++){
	CN = Index_Rcv_Grid_B2C[IDR][LN];
	data_C[CN] = Work_Array_Snd_Grid_B2C[gp+LN];
      } 
    }
    else{

      gp = GP_B2C_R[ID];
      for (LN=0; LN<Num_Rcv_Grid_B2C[IDR]; LN++){
	CN = Index_Rcv_Grid_B2C[IDR][LN];
	data_C[CN] = Work_Array_Rcv_Grid_B2C[gp+LN];
      }
    }
  }

  /* freeing of arrays */
  free(Work_Array_Snd_Grid_B2C);
  free(Work_Array_Rcv_Grid_B2C);
}





void Density_Grid_Copy_B2D(double **Density_Grid_B0)
{
  static int firsttime=1;
  int DN,BN,LN,spin,i,gp,NN_S,NN_R;
  double *Work_Array_Snd_Grid_B2D;
  double *Work_Array_Rcv_Grid_B2D;
  int numprocs,myid,tag=999,ID,IDS,IDR;
  MPI_Status stat;
  MPI_Request request;
  MPI_Status *stat_send;
  MPI_Status *stat_recv;
  MPI_Request *request_send;
  MPI_Request *request_recv;

  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  /* allocation of arrays */
  
  Work_Array_Snd_Grid_B2D = (double*)malloc(sizeof(double)*GP_B2D_S[NN_B2D_S]*(SpinP_switch+1)); 
  Work_Array_Rcv_Grid_B2D = (double*)malloc(sizeof(double)*GP_B2D_R[NN_B2D_R]*(SpinP_switch+1)); 

  if (firsttime==1){
    PrintMemory("Set_Density_Grid: Work_Array_Snd_Grid_B2D",
		sizeof(double)*GP_B2D_S[NN_B2D_S]*(SpinP_switch+1), NULL);
    PrintMemory("Set_Density_Grid: Work_Array_Rcv_Grid_B2D",
		sizeof(double)*GP_B2D_R[NN_B2D_R]*(SpinP_switch+1), NULL);
    firsttime = 0;
  }

  /******************************************************
             MPI: from the partitions B to D
  ******************************************************/

  request_send = malloc(sizeof(MPI_Request)*NN_B2D_S);
  request_recv = malloc(sizeof(MPI_Request)*NN_B2D_R);
  stat_send = malloc(sizeof(MPI_Status)*NN_B2D_S);
  stat_recv = malloc(sizeof(MPI_Status)*NN_B2D_R);

  NN_S = 0;
  NN_R = 0;

  /* MPI_Irecv */

  for (ID=0; ID<NN_B2D_R; ID++){

    IDR = ID_NN_B2D_R[ID];
    gp = GP_B2D_R[ID];

    if (IDR!=myid){ 
      MPI_Irecv( &Work_Array_Rcv_Grid_B2D[(SpinP_switch+1)*gp], Num_Rcv_Grid_B2D[IDR]*(SpinP_switch+1),
                 MPI_DOUBLE, IDR, tag, mpi_comm_level1, &request_recv[NN_R]);
      NN_R++;
    }
  }

  /* MPI_Isend */

  for (ID=0; ID<NN_B2D_S; ID++){

    IDS = ID_NN_B2D_S[ID];
    gp = GP_B2D_S[ID];

    /* copy Density_Grid_B0 to Work_Array_Snd_Grid_B2D */

    for (LN=0; LN<Num_Snd_Grid_B2D[IDS]; LN++){

      BN = Index_Snd_Grid_B2D[IDS][LN];

      if (SpinP_switch==0){
        Work_Array_Snd_Grid_B2D[gp+LN]       = Density_Grid_B0[0][BN];
      }
      else if (SpinP_switch==1){
        Work_Array_Snd_Grid_B2D[2*gp+2*LN+0] = Density_Grid_B0[0][BN];
        Work_Array_Snd_Grid_B2D[2*gp+2*LN+1] = Density_Grid_B0[1][BN];
      }
      else if (SpinP_switch==3){
        Work_Array_Snd_Grid_B2D[4*gp+4*LN+0] = Density_Grid_B0[0][BN];
        Work_Array_Snd_Grid_B2D[4*gp+4*LN+1] = Density_Grid_B0[1][BN];
        Work_Array_Snd_Grid_B2D[4*gp+4*LN+2] = Density_Grid_B0[2][BN];
        Work_Array_Snd_Grid_B2D[4*gp+4*LN+3] = Density_Grid_B0[3][BN];
      }
    } /* LN */        

    if (IDS!=myid){
      MPI_Isend( &Work_Array_Snd_Grid_B2D[(SpinP_switch+1)*gp], Num_Snd_Grid_B2D[IDS]*(SpinP_switch+1), 
		 MPI_DOUBLE, IDS, tag, mpi_comm_level1, &request_send[NN_S]);
      NN_S++;
    }
  }

  /* MPI_Waitall */

  if (NN_S!=0) MPI_Waitall(NN_S,request_send,stat_send);
  if (NN_R!=0) MPI_Waitall(NN_R,request_recv,stat_recv);

  free(request_send);
  free(request_recv);
  free(stat_send);
  free(stat_recv);

  /* copy Work_Array_Rcv_Grid_B2D to Density_Grid_D */

  for (ID=0; ID<NN_B2D_R; ID++){

    IDR = ID_NN_B2D_R[ID];

    if (IDR==myid){

      gp = GP_B2D_S[ID];

      for (LN=0; LN<Num_Rcv_Grid_B2D[IDR]; LN++){

	DN = Index_Rcv_Grid_B2D[IDR][LN];

	if (SpinP_switch==0){
	  Density_Grid_D[0][DN] = Work_Array_Snd_Grid_B2D[gp+LN];
	}     
	else if (SpinP_switch==1){
	  Density_Grid_D[0][DN] = Work_Array_Snd_Grid_B2D[2*gp+2*LN+0];
	  Density_Grid_D[1][DN] = Work_Array_Snd_Grid_B2D[2*gp+2*LN+1];
	}     
	else if (SpinP_switch==3){
	  Density_Grid_D[0][DN] = Work_Array_Snd_Grid_B2D[4*gp+4*LN+0];
	  Density_Grid_D[1][DN] = Work_Array_Snd_Grid_B2D[4*gp+4*LN+1];
	  Density_Grid_D[2][DN] = Work_Array_Snd_Grid_B2D[4*gp+4*LN+2];
	  Density_Grid_D[3][DN] = Work_Array_Snd_Grid_B2D[4*gp+4*LN+3];
	}
      } /* LN */   

    }

    else{

      gp = GP_B2D_R[ID];

      for (LN=0; LN<Num_Rcv_Grid_B2D[IDR]; LN++){

	DN = Index_Rcv_Grid_B2D[IDR][LN];

	if (SpinP_switch==0){
	  Density_Grid_D[0][DN] = Work_Array_Rcv_Grid_B2D[gp+LN];
	}     
	else if (SpinP_switch==1){
	  Density_Grid_D[0][DN] = Work_Array_Rcv_Grid_B2D[2*gp+2*LN+0];
	  Density_Grid_D[1][DN] = Work_Array_Rcv_Grid_B2D[2*gp+2*LN+1];
	}     
	else if (SpinP_switch==3){
	  Density_Grid_D[0][DN] = Work_Array_Rcv_Grid_B2D[4*gp+4*LN+0];
	  Density_Grid_D[1][DN] = Work_Array_Rcv_Grid_B2D[4*gp+4*LN+1];
	  Density_Grid_D[2][DN] = Work_Array_Rcv_Grid_B2D[4*gp+4*LN+2];
	  Density_Grid_D[3][DN] = Work_Array_Rcv_Grid_B2D[4*gp+4*LN+3];
	}
      }

    }
  }

  /* if (SpinP_switch==0), copy Density_Grid to Density_Grid */

  if (SpinP_switch==0){
    for (BN=0; BN<My_NumGridB_AB; BN++){
      Density_Grid_B0[1][BN] = Density_Grid_B0[0][BN]; 
    }

    for (DN=0; DN<My_NumGridD; DN++){
      Density_Grid_D[1][DN] = Density_Grid_D[0][DN]; 
    }
  }

  /* freeing of arrays */
  free(Work_Array_Snd_Grid_B2D);
  free(Work_Array_Rcv_Grid_B2D);
}


void diagonalize_nc_density(double **Density_Grid_B0)
{
  int BN,DN,Mc_AN,Gc_AN,Nog,GRc;
  double Re11,Re22,Re12,Im12;
  double phi[2],theta[2],sit,cot,sip,cop;
  double d1,d2,d3,x,y,z,Cxyz[4];
  double Nup[2],Ndown[2];
  /* for OpenMP */
  int OMPID,Nthrds;

  /************************************
     Density_Grid in the partition B
  ************************************/

#pragma omp parallel shared(Density_Grid_B0,My_NumGridB_AB) private(OMPID,Nthrds,BN,Re11,Re22,Re12,Im12,Nup,Ndown,theta,phi) default(none)
  {

    /* get info. on OpenMP */ 

    OMPID = omp_get_thread_num();
    Nthrds = omp_get_num_threads();

    for (BN=OMPID; BN<My_NumGridB_AB; BN+=Nthrds){

      Re11 = Density_Grid_B0[0][BN];
      Re22 = Density_Grid_B0[1][BN];
      Re12 = Density_Grid_B0[2][BN];
      Im12 = Density_Grid_B0[3][BN];

      EulerAngle_Spin( 1, Re11, Re22, Re12, Im12, Re12, -Im12, Nup, Ndown, theta, phi );

      Density_Grid_B0[0][BN] = Nup[0];
      Density_Grid_B0[1][BN] = Ndown[0];
      Density_Grid_B0[2][BN] = theta[0];
      Density_Grid_B0[3][BN] = phi[0];
    }

#pragma omp flush(Density_Grid_B0)

  } /* #pragma omp parallel */

  /************************************
     Density_Grid in the partition D
  ************************************/

#pragma omp parallel shared(Density_Grid_D,My_NumGridD) private(OMPID,Nthrds,DN,Re11,Re22,Re12,Im12,Nup,Ndown,theta,phi) default(none)
  {

    /* get info. on OpenMP */ 

    OMPID = omp_get_thread_num();
    Nthrds = omp_get_num_threads();

    for (DN=OMPID; DN<My_NumGridD; DN+=Nthrds){

      Re11 = Density_Grid_D[0][DN];
      Re22 = Density_Grid_D[1][DN];
      Re12 = Density_Grid_D[2][DN];
      Im12 = Density_Grid_D[3][DN];

      EulerAngle_Spin( 1, Re11, Re22, Re12, Im12, Re12, -Im12, Nup, Ndown, theta, phi );

      Density_Grid_D[0][DN] = Nup[0];
      Density_Grid_D[1][DN] = Ndown[0];
      Density_Grid_D[2][DN] = theta[0];
      Density_Grid_D[3][DN] = phi[0];
    }

#pragma omp flush(Density_Grid_D)

  } /* #pragma omp parallel */

}
