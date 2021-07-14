/**********************************************************************
  Set_ProExpn_VNA.c:

     Set_ProExpn_VNA.c is a subroutine to calculate matrix elements and
     the derivatives of VNA projector expansion in the momentum space.

  Log of Set_ProExpn_VNA.c:

     8/Apr/2004  Released by T.Ozaki

***********************************************************************/

#define  measure_time   0

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <sys/types.h>
#include <sys/times.h>
#include <sys/time.h> 
#include <sys/stat.h>
#include <unistd.h>
#include "openmx_common.h"
#include "mpi.h"
#include <omp.h>


#ifdef _FTRACE
#include <ftrace.h>
#else
//dummy macro
#define ftrace_region_begin(x)   
#define ftrace_region_end(x)   
#endif

static double Set_ProExpn(double ****HVNA, Type_DS_VNA *****DS_VNA);
static double Set_VNA2(double ****HVNA, double *****HVNA2);
static double Set_VNA3(double *****HVNA3);

#ifdef kcomp
static void Spherical_Bessel2( double x, int lmax, double *sb, double *dsb );
#elif defined (_NEC)
inline static void Spherical_Bessel2( double x, int lmax, double *sb, double *dsb );
#else 
inline void Spherical_Bessel2( double x, int lmax, double *sb, double *dsb );
/*
inline static dcomplex**** Allocate4D_dcomplex(int size_1, int size_2, int size_3, int size_4);
inline static dcomplex** Allocate2D_dcomplex(int size_1, int size_2);
inline static void Free4D_dcomplex(dcomplex**** buffer);
inline static void Free2D_dcomplex(dcomplex** buffer);
*/
#endif     

#ifdef _NEC
#define _NEC3
#define _NEC2
#define _NEC1
#define _NEC_ARRAY
#define _NEC2_ARRAY
#define _NEC1_ARRAY
//#define _NEC_BESSEL_WITH_SAFTY  //this relates Spherical_Bessel2array, but just numerical safty, and this blocks vecterization//
static inline void Spherical_Bessel2array(int x_size, double* x_array, double coef, int lmax, double *sb, double *dsb );

#include "Gaunt.h"

#define _NEC_MG

#define AI_L0_UNROLL_SIZE  (7)      /*(3*2 + 1)*/
#define AI_L1_UNROLL_SIZE   (19)     /* 19 = ((3+6)*2 +1) */
#define _NEC_LIST_V4
//#define _NEC_LIST_V4NL   //this is not fast//
#endif

double Set_ProExpn_VNA(double ****HVNA, double *****HVNA2, Type_DS_VNA *****DS_VNA)
{
  double time1,time2,time3;

  /* separable form */

  time1 = Set_ProExpn(HVNA, DS_VNA);

  /* one-center (orbitals) but two-center integrals, <Phi_{LM,L'M'}|VNA> */

  time2 = Set_VNA2(HVNA, HVNA2);

  /* one-center (orbitals) but two-center integrals, <VNA|Phi_{LM,L'M'}> */

  time3 = Set_VNA3(HVNA3);

  if (measure_time){
    printf("Time Set_ProExpn=%7.3f Set_VNA2=%7.3f Set_VNA3=%7.3f\n",time1,time2,time3);
  }

  return time1+time2+time3;
}



double Set_ProExpn(double ****HVNA, Type_DS_VNA *****DS_VNA)
{
  /****************************************************
   evaluate matrix elements of neutral atom potential
   in the KB or Blochl separable form in momentum space.
  ****************************************************/
  static int firsttime=1;
  int L2,L3,L,i,j,Mc_AN,Gc_AN,h_AN,Cwan,i1,j1;
  int kk,Ls,n,j0,jg0,Mj_AN0,m,L1,GL,Mul1;
  int tno0,tno1,tno2,p,q,po1,po,Original_Mc_AN;
  int L0,Mul0,spe,Gh_AN,Hwan,fan,jg,k,kg,wakg,kl;
  int size_SumNL0,size_TmpNL;
  int Mc_AN2,Gc_AN2,num,size1,size2;
  int *Snd_DS_VNA_Size,*Rcv_DS_VNA_Size;  
  double kmin,kmax,Sk,Dk,Normk,dmp;
  Type_DS_VNA *tmp_array;
  Type_DS_VNA *tmp_array2;
  double tmp0,tmp3,tmp4,tmp5,tmp10;
  double ****Bessel_Pro00;
  double ****Bessel_Pro01;
  double ene,sum,h,coe0,sj,sy,sjp,syp;
  double rcutA,rcutB,rcut;
  double TStime,TEtime;
  double stime,etime;
#ifdef _NEC
  int mn;
#endif
  double Stime_atom,Etime_atom;
  int numprocs,myid,tag=999,ID,IDS,IDR;
  int *VNA_List;
  int *VNA_List2;
  int Num_RVNA;
  double **NLH;
  double time0,time1,time2,time3,time4,time5;
  dcomplex sum0,sum1,sum2; 

  MPI_Status stat;
  MPI_Request request;

  int OneD_Nloop,*OneD2Mc_AN,*OneD2h_AN;

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  if (measure_time){
    time1 = 0.0;
    time2 = 0.0;
    time3 = 0.0;
    time4 = 0.0;
    time5 = 0.0;
    dtime(&stime);
  }

  dtime(&TStime);

  /****************************************************
                 allocation of arrays:
  ****************************************************/
  Num_RVNA = List_YOUSO[34]*(List_YOUSO[35] + 1);

  Snd_DS_VNA_Size = (int*)malloc(sizeof(int)*numprocs);
  Rcv_DS_VNA_Size = (int*)malloc(sizeof(int)*numprocs);
#ifndef _NEC_LIST_V4
  VNA_List  = (int*)malloc(sizeof(int)*(List_YOUSO[34]*(List_YOUSO[35] + 1)+2) ); 
  VNA_List2 = (int*)malloc(sizeof(int)*(List_YOUSO[34]*(List_YOUSO[35] + 1)+2) ); 
#ifdef _NEC
  int* VNA_offset= (int*)malloc(sizeof(int)*(List_YOUSO[34]*(List_YOUSO[35] + 1)+2) ); 
#endif
#endif
#ifdef _NEC_LIST_LM
  /* size is (L+1)^2 * (List_YOUSO[34]) */
  int* ai_VNA_List_L1 = (int*)malloc(sizeof(int)*(List_YOUSO[35]+1)*(List_YOUSO[35]+1)*(List_YOUSO[34]));
  int* ai_VNA_List_Mul1 = (int*)malloc(sizeof(int)*(List_YOUSO[35]+1)*(List_YOUSO[35]+1)*(List_YOUSO[34]));
  int* ai_VNA_List_L1M1 = (int*)malloc(sizeof(int)*(List_YOUSO[35]+1)*(List_YOUSO[35]+1)*(List_YOUSO[34]));
#endif

  Bessel_Pro00 = (double****)malloc(sizeof(double***)*SpeciesNum); 
  for (spe=0; spe<SpeciesNum; spe++){
    Bessel_Pro00[spe] = (double***)malloc(sizeof(double**)*(Spe_MaxL_Basis[spe]+1)); 
    for (L0=0; L0<=Spe_MaxL_Basis[spe]; L0++){
      Bessel_Pro00[spe][L0] = (double**)malloc(sizeof(double*)*Spe_Num_Basis[spe][L0]); 
      for (Mul0=0; Mul0<Spe_Num_Basis[spe][L0]; Mul0++){
	Bessel_Pro00[spe][L0][Mul0] = (double*)malloc(sizeof(double)*GL_Mesh); 
      }
    }
  }

  Bessel_Pro01 = (double****)malloc(sizeof(double***)*SpeciesNum); 
  for (spe=0; spe<SpeciesNum; spe++){
    Bessel_Pro01[spe] = (double***)malloc(sizeof(double**)*(Spe_MaxL_Basis[spe]+1)); 
    for (L0=0; L0<=Spe_MaxL_Basis[spe]; L0++){
      Bessel_Pro01[spe][L0] = (double**)malloc(sizeof(double*)*Spe_Num_Basis[spe][L0]); 
      for (Mul0=0; Mul0<Spe_Num_Basis[spe][L0]; Mul0++){
        Bessel_Pro01[spe][L0][Mul0] = (double*)malloc(sizeof(double)*GL_Mesh); 
      }
    }
  }

  size_TmpNL = (List_YOUSO[25]+1)*List_YOUSO[24]*
               (2*(List_YOUSO[25]+1)+1)*Num_RVNA*(2*List_YOUSO[30]+1);
  size_SumNL0 = (List_YOUSO[25]+1)*List_YOUSO[24]*(Num_RVNA+1);

  /* PrintMemory */
  if (firsttime) {
    PrintMemory("Set_ProExpn_VNA: SumNL0",sizeof(double)*size_SumNL0,NULL);
    PrintMemory("Set_ProExpn_VNA: SumNLr0",sizeof(double)*size_SumNL0,NULL);
    PrintMemory("Set_ProExpn_VNA: TmpNL", sizeof(dcomplex)*size_TmpNL,NULL);
    PrintMemory("Set_ProExpn_VNA: TmpNLr",sizeof(dcomplex)*size_TmpNL,NULL);
    PrintMemory("Set_ProExpn_VNA: TmpNLt",sizeof(dcomplex)*size_TmpNL,NULL);
    PrintMemory("Set_ProExpn_VNA: TmpNLp",sizeof(dcomplex)*size_TmpNL,NULL);
    firsttime=0;
  }

  if (measure_time){
	  dtime(&etime);
	  time1 += etime - stime;
    dtime(&stime);
  }


  /************************************************************
   Spe_Num_RVPS[Hwan]     ->  Num_RVNA = List_YOUSO[34]*(List_YOUSO[35]+1) 
   Spe_VPS_List[Hwan][L]  ->  VNA_List[L] = 0,0,0,.., 1,1,1,.., 2,2,2,..,
   List_YOUSO[19]         ->  Num_RVNA
   List_YOUSO[30]         ->  List_YOUSO[35]
  *************************************************************/

#ifdef _NEC_LIST_V4
  const int width_Mul = List_YOUSO[34];
  const int Num_RVNA_with_M = List_YOUSO[34] * (List_YOUSO[35]+1)*(List_YOUSO[35]+1);//sizeof i,j,and 2i*1//
  const int L1_max = List_YOUSO[35];
#else
  L = 0;
  int offset = 0;
  for (i=0; i<=List_YOUSO[35]; i++){    /* max L */
    for (j=0; j<List_YOUSO[34]; j++){   /* # of radial projectors */
      VNA_List[L]  = i;
      VNA_List2[L] = j;
#ifdef _NEC      
      VNA_offset[L] = offset; //offset of L that is sum of (L1*2+1)//
      offset += 2*i + 1;
#endif      
      L++;
    }
  }
#endif
  /************************************************************
   pre-calculate the time consuming part in the 'time1' part
  *************************************************************/

  kmin = Radial_kmin;
  kmax = PAO_Nkmax;
  Sk = kmax + kmin;
  Dk = kmax - kmin;

  for (spe=0; spe<SpeciesNum; spe++){

    for (L0=0; L0<=Spe_MaxL_Basis[spe]; L0++){
      for (Mul0=0; Mul0<Spe_Num_Basis[spe][L0]; Mul0++){

//#pragma omp parallel shared(Mul0,L0,spe,Bessel_Pro01,Bessel_Pro00,GL_Weight,Dk,GL_NormK) private(i,Normk,tmp10)
	{         

	  int OMPID,Nthrds,Nprocs;

	  /* get info. on OpenMP */ 

	  OMPID = omp_get_thread_num();
	  Nthrds = omp_get_num_threads();
	  Nprocs = omp_get_num_procs();

	  for (i=OMPID*GL_Mesh/Nthrds; i<(OMPID+1)*GL_Mesh/Nthrds; i++){

	    Normk = GL_NormK[i];
	    tmp10 = 0.50*Dk*GL_Weight[i]*Normk*Normk;

	    Bessel_Pro00[spe][L0][Mul0][i] = RF_BesselF(spe,L0,Mul0,Normk)*tmp10;
	    Bessel_Pro01[spe][L0][Mul0][i] = Bessel_Pro00[spe][L0][Mul0][i]*Normk;
	  }

#pragma omp flush(Bessel_Pro00,Bessel_Pro01)

	} /* #pragma omp parallel */

      }
    }
  }
  /************************************************************
    start the main calculation
  *************************************************************/

  /* one-dimensionalize the Mc_AN and h_AN loops */

  OneD_Nloop = 0;
  for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
    Gc_AN = M2G[Mc_AN];    
    for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
      OneD_Nloop++;
    }
  }  

  OneD2Mc_AN = (int*)malloc(sizeof(int)*(OneD_Nloop+1));
  OneD2h_AN = (int*)malloc(sizeof(int)*(OneD_Nloop+1));

  OneD_Nloop = 0;
  for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
    Gc_AN = M2G[Mc_AN];    
    for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
      OneD2Mc_AN[OneD_Nloop] = Mc_AN; 
      OneD2h_AN[OneD_Nloop] = h_AN; 
      OneD_Nloop++;
    }
  }

#ifdef _NEC
ftrace_region_begin("MG1");

  int max_L0 = 0;
  for (int spe=0; spe<SpeciesNum; spe++){
    if(max_L0 < Spe_MaxL_Basis[spe]) max_L0 = Spe_MaxL_Basis[spe];
  }
  const int ai_L0_LIMIT = max_L0 + 1;
  const int ai_L0_LIMIT_SQ = ai_L0_LIMIT * ai_L0_LIMIT;
  const int ai_L1_LIMIT = List_YOUSO[35] + 1;
  const int ai_L1_LIMIT_SQ = ai_L1_LIMIT*ai_L1_LIMIT;
  const int ai_LL_LIMIT = (ai_L0_LIMIT > ai_L1_LIMIT) ? ai_L0_LIMIT*2-1 : ai_L1_LIMIT*2-1;
  const int ai_LL_LIMIT_SQ = ai_LL_LIMIT*ai_LL_LIMIT;
  double* Gaunt_list = MakeGaunt(ai_L0_LIMIT_SQ, ai_L1_LIMIT_SQ, ai_LL_LIMIT_SQ);
ftrace_region_end("MG1");
#endif

  if (measure_time){
	  dtime(&etime);
	  time2 += etime - stime;
    dtime(&stime);
  }

#pragma omp parallel shared(time_per_atom,DS_VNA,Comp2Real,Spe_VNA_Bessel,Bessel_Pro01,Bessel_Pro00,GL_NormK,List_YOUSO,VNA_List,Num_RVNA,Spe_Num_Basis,Spe_MaxL_Basis,PAO_Nkmax,atv,Gxyz,ncn,natn,Spe_Atom_Cut1,WhatSpecies,M2G,OneD2h_AN,OneD2Mc_AN,OneD_Nloop,time1,time2,time3) 

  {
    int i,j,k,l,m,num0,num1;
    int Mc_AN,h_AN,Gc_AN,Cwan,Gh_AN; 
    int Rnh,Hwan,L0,Mul0,L,L1,M0,M1,LL,Mul1,Ls;
    int OMPID,Nthrds,Nprocs,Nloop;
    int Lmax,Lmax_Four_Int;
    double Stime_atom,Etime_atom;
    double rcutA,rcutB,rcut;
    double dx,dy,dz;
    double kmin,kmax,Sk,Dk;
    double gant,SH[2],dSHt[2],dSHp[2];
    double S_coordinate[3];
    double r,theta,phi;
    double Normk,tmp0,tmp1,tmp2;
    double siT,coT,siP,coP;
    double *SumNL0;
    double *SumNLr0;
    double *Bes00;
    double *Bes01;
#ifndef _NEC3    
    double tmp_SphB[30];
    double tmp_SphBp[30];
    double SphB[30][GL_Mesh];
    double SphBp[30][GL_Mesh];
#endif    
    double stime,etime;

    dcomplex Ctmp1,Ctmp0,Ctmp2;
    dcomplex CsumNL0,CsumNLr,CsumNLt,CsumNLp;
    dcomplex CY,CYt,CYp,CY1,CYt1,CYp1,Cpow;
#ifdef _NEC_LIST_V4
    dcomplex **TmpNL;
    dcomplex **TmpNLr;
    dcomplex **TmpNLt;
    dcomplex **TmpNLp;
#else
    dcomplex ***TmpNL;
    dcomplex ***TmpNLr;
    dcomplex ***TmpNLt;
    dcomplex ***TmpNLp;
#endif    
    dcomplex *CmatNL0;
    dcomplex *CmatNLr;
    dcomplex *CmatNLt;
    dcomplex *CmatNLp;

    /* allocation of arrays */
#ifdef _NEC_LIST_V4
	TmpNL = (dcomplex**)malloc(sizeof(dcomplex*)*(2*(List_YOUSO[25]+1)+1));
  TmpNL[0] = (dcomplex*)malloc(sizeof(dcomplex)*(2*(List_YOUSO[25]+1)+1)*Num_RVNA_with_M);
  for (k=0; k<(2*(List_YOUSO[25]+1)+1); k++){
	  TmpNL[k] = TmpNL[0] + k * Num_RVNA_with_M;	  
	}

	TmpNLr = (dcomplex**)malloc(sizeof(dcomplex*)*(2*(List_YOUSO[25]+1)+1));
  TmpNLr[0] = (dcomplex*)malloc(sizeof(dcomplex)*(2*(List_YOUSO[25]+1)+1)*Num_RVNA_with_M);
  for (k=0; k<(2*(List_YOUSO[25]+1)+1); k++){
	  TmpNLr[k] = TmpNLr[0] + k * Num_RVNA_with_M;	  
	}
  
	TmpNLt = (dcomplex**)malloc(sizeof(dcomplex*)*(2*(List_YOUSO[25]+1)+1));
  TmpNLt[0] = (dcomplex*)malloc(sizeof(dcomplex)*(2*(List_YOUSO[25]+1)+1)*Num_RVNA_with_M);
  for (k=0; k<(2*(List_YOUSO[25]+1)+1); k++){
	  TmpNLt[k] = TmpNLt[0] + k * Num_RVNA_with_M;	  
	}
  
	TmpNLp = (dcomplex**)malloc(sizeof(dcomplex*)*(2*(List_YOUSO[25]+1)+1));
  TmpNLp[0] = (dcomplex*)malloc(sizeof(dcomplex)*(2*(List_YOUSO[25]+1)+1)*Num_RVNA_with_M);
  for (k=0; k<(2*(List_YOUSO[25]+1)+1); k++){
	  TmpNLp[k] = TmpNLp[0] + k * Num_RVNA_with_M;	  
	}
#else
	TmpNL = (dcomplex***)malloc(sizeof(dcomplex**)*(2*(List_YOUSO[25]+1)+1));
	for (k=0; k<(2*(List_YOUSO[25]+1)+1); k++){
	  TmpNL[k] = (dcomplex**)malloc(sizeof(dcomplex*)*Num_RVNA);
	  for (l=0; l<Num_RVNA; l++){
	    TmpNL[k][l] = (dcomplex*)malloc(sizeof(dcomplex)*(2*List_YOUSO[35]+1));
	  }
	}
    
	TmpNLr = (dcomplex***)malloc(sizeof(dcomplex**)*(2*(List_YOUSO[25]+1)+1));
	for (k=0; k<(2*(List_YOUSO[25]+1)+1); k++){
	  TmpNLr[k] = (dcomplex**)malloc(sizeof(dcomplex*)*Num_RVNA);
	  for (l=0; l<Num_RVNA; l++){
	    TmpNLr[k][l] = (dcomplex*)malloc(sizeof(dcomplex)*(2*List_YOUSO[35]+1));
	  }
	}

	TmpNLt = (dcomplex***)malloc(sizeof(dcomplex**)*(2*(List_YOUSO[25]+1)+1));
	for (k=0; k<(2*(List_YOUSO[25]+1)+1); k++){
	  TmpNLt[k] = (dcomplex**)malloc(sizeof(dcomplex*)*Num_RVNA);
	  for (l=0; l<Num_RVNA; l++){
	    TmpNLt[k][l] = (dcomplex*)malloc(sizeof(dcomplex)*(2*List_YOUSO[35]+1));
	  }
	}

	TmpNLp = (dcomplex***)malloc(sizeof(dcomplex**)*(2*(List_YOUSO[25]+1)+1));
	for (k=0; k<(2*(List_YOUSO[25]+1)+1); k++){
	  TmpNLp[k] = (dcomplex**)malloc(sizeof(dcomplex*)*Num_RVNA);
	  for (l=0; l<Num_RVNA; l++){
	    TmpNLp[k][l] = (dcomplex*)malloc(sizeof(dcomplex)*(2*List_YOUSO[35]+1));
	  }
	}
#endif
#ifdef _NEC_LIST_V4NL    
	  SumNL0 = (double*)malloc(sizeof(double)*Num_RVNA_with_M);
    
    SumNLr0 = (double*)malloc(sizeof(double)*Num_RVNA_with_M);
#else
	  SumNL0 = (double*)malloc(sizeof(double)*Num_RVNA);
    
    SumNLr0 = (double*)malloc(sizeof(double)*Num_RVNA);
#endif    
#ifdef _NEC_LIST_V4    
    CmatNL0 = (dcomplex*)malloc(sizeof(dcomplex)*Num_RVNA_with_M);

    CmatNLr = (dcomplex*)malloc(sizeof(dcomplex)*Num_RVNA_with_M);
    
    CmatNLt = (dcomplex*)malloc(sizeof(dcomplex)*Num_RVNA_with_M);

    CmatNLp = (dcomplex*)malloc(sizeof(dcomplex)*Num_RVNA_with_M);
#else
    CmatNL0 = (dcomplex*)malloc(sizeof(dcomplex)*(2*List_YOUSO[35]+1));

    CmatNLr = (dcomplex*)malloc(sizeof(dcomplex)*(2*List_YOUSO[35]+1));
    
    CmatNLt = (dcomplex*)malloc(sizeof(dcomplex)*(2*List_YOUSO[35]+1));

    CmatNLp = (dcomplex*)malloc(sizeof(dcomplex)*(2*List_YOUSO[35]+1));
#endif


    Bes00 = (double*)malloc(sizeof(double)*GL_Mesh); 
    Bes01 = (double*)malloc(sizeof(double)*GL_Mesh); 
#ifdef _NEC3
    double* tmp_SphB_D  = (double*)malloc(sizeof(double)*30*GL_Mesh);
    double* tmp_SphBp_D = (double*)malloc(sizeof(double)*30*GL_Mesh);
#endif
    /* get info. on OpenMP */ 

    OMPID = omp_get_thread_num();
    Nthrds = omp_get_num_threads();
    Nprocs = omp_get_num_procs();

    /* one-dimensionalized loop */

    for (Nloop=OMPID*OneD_Nloop/Nthrds; Nloop<(OMPID+1)*OneD_Nloop/Nthrds; Nloop++){

      dtime(&Stime_atom);
ftrace_region_begin("R1");
      /* get Mc_AN and h_AN */

      Mc_AN = OneD2Mc_AN[Nloop];
      h_AN  = OneD2h_AN[Nloop];

      /* set data on Mc_AN */

      Gc_AN = M2G[Mc_AN];
      Cwan = WhatSpecies[Gc_AN];
      rcutA = Spe_Atom_Cut1[Cwan];

      /* set data on h_AN */

      Gh_AN = natn[Gc_AN][h_AN];        
      Rnh = ncn[Gc_AN][h_AN];
      Hwan = WhatSpecies[Gh_AN];
      rcutB = Spe_Atom_Cut1[Hwan];
      rcut = rcutA + rcutB;

      dx = Gxyz[Gh_AN][1] + atv[Rnh][1] - Gxyz[Gc_AN][1]; 
      dy = Gxyz[Gh_AN][2] + atv[Rnh][2] - Gxyz[Gc_AN][2]; 
      dz = Gxyz[Gh_AN][3] + atv[Rnh][3] - Gxyz[Gc_AN][3];

      xyz2spherical(dx,dy,dz,0.0,0.0,0.0,S_coordinate); 
      r     = S_coordinate[0];
      theta = S_coordinate[1];
      phi   = S_coordinate[2];

      /* for empty atoms or finite elemens basis */

      if (r<1.0e-10) r = 1.0e-10;

      /* precalculation of sin and cos */

      siT = sin(theta);
      coT = cos(theta);
      siP = sin(phi);
      coP = cos(phi);

      /****************************************************
         evaluate ovelap integrals <chi0|P> between PAOs
         and progectors of nonlocal VNA potentials. 
      ****************************************************/
      /****************************************************
              \int RL(k)*RL'(k)*jl(k*R) k^2 dk^3 
      ****************************************************/

      kmin = Radial_kmin;
      kmax = PAO_Nkmax;
      Sk = kmax + kmin;
      Dk = kmax - kmin;

ftrace_region_end("R1");
ftrace_region_begin("R2");


      Lmax = -10;
      for (L=0; L<Num_RVNA; L++){
#ifdef _NEC_LIST_V4
        const int L1 = L / width_Mul;
        if (Lmax<L1) Lmax = L1;
#else
	      if (Lmax<VNA_List[L]) Lmax = VNA_List[L];
#endif
      }
      if (Spe_MaxL_Basis[Cwan]<Lmax)
	      Lmax_Four_Int = 2*Lmax;
      else 
        Lmax_Four_Int = 2*Spe_MaxL_Basis[Cwan];

/* List_YOUSO[35] and Lmax expect 7 /*
/* printf("AITUNE: Lmax = %d, %d, %d\n", Lmax, Spe_MaxL_Basis[Cwan], List_YOUSO[35]); */


      /* calculate SphB and SphBp */
#ifdef _NEC3
#ifdef _NEC_ARRAY
      Spherical_Bessel2array(GL_Mesh, GL_NormK, r, Lmax_Four_Int, &(tmp_SphB_D[0]),&(tmp_SphBp_D[0]));
#else
      #pragma _NEC vector
      for (i=0; i<GL_Mesh; i++){
        const double Normk = GL_NormK[i];
        Spherical_Bessel2(Normk*r,Lmax_Four_Int,&(tmp_SphB_D[i*(Lmax_Four_Int+3)]),&(tmp_SphBp_D[i*(Lmax_Four_Int+3)]));
      }
#endif    
#else      
#ifdef kcomp
#else 
#pragma forceinline recursive
#endif      
      for (i=0; i<GL_Mesh; i++){
        Normk = GL_NormK[i];
        Spherical_Bessel2(Normk*r,Lmax_Four_Int,tmp_SphB,tmp_SphBp);
        for(LL=0; LL<=Lmax_Four_Int; LL++){ 
          SphB[LL][i]  = tmp_SphB[LL]; 
          SphBp[LL][i] = tmp_SphBp[LL]; 
	      }
      }
#endif        

ftrace_region_end("R2");

      /* LL loop */
      
      

      num0 = 0;      
      for (L0=0; L0<=Spe_MaxL_Basis[Cwan]; L0++){
        for (Mul0=0; Mul0<Spe_Num_Basis[Cwan][L0]; Mul0++){


ftrace_region_begin("R3");
          for (M0=-L0; M0<=L0; M0++){
#ifdef _NEC_LIST_V4
            for (int ilm=0; ilm<Num_RVNA_with_M; ilm++){
              //const int L1 = (int)sqrt(0.1 + (L / width_Mul));
              //const int offset = VNA_offset[L];
              TmpNL[L0+M0][ilm] = Complex(0.0,0.0); 
              TmpNLr[L0+M0][ilm] = Complex(0.0,0.0); 
              TmpNLt[L0+M0][ilm] = Complex(0.0,0.0); 
              TmpNLp[L0+M0][ilm] = Complex(0.0,0.0); 
            }
#else

            for (L=0; L<Num_RVNA; L++){
            const int L1 = VNA_List[L];
            const int offset = VNA_offset[L];
              for (M1=-L1; M1<=L1; M1++){
                TmpNL[L0+M0][L][L1+M1] = Complex(0.0,0.0); 
                TmpNLr[L0+M0][L][L1+M1] = Complex(0.0,0.0); 
                TmpNLt[L0+M0][L][L1+M1] = Complex(0.0,0.0); 
                TmpNLp[L0+M0][L][L1+M1] = Complex(0.0,0.0); 
              }
            }
#endif
          }

ftrace_region_end("R3");
          
          for(LL=0; LL<=Lmax_Four_Int; LL++){ 
          
ftrace_region_begin("R4.1");
#ifdef _NEC_LIST_V4NL
            for (int iml=0; iml<Num_RVNA_with_M; iml++){
              SumNL0[iml] = 0.0;
              SumNLr0[iml] = 0.0;
            }
#else
            for (L=0; L<Num_RVNA; L++){
              SumNL0[L] = 0.0;
              SumNLr0[L] = 0.0;
            }
#endif
ftrace_region_end("R4.1");
ftrace_region_begin("R4.2");

        /* Gauss-Legendre quadrature */
#ifdef _NEC3
#ifdef _NEC_ARRAY
            #pragma _NEC vector
            for (i=0; i<GL_Mesh; i++){
              Bes00[i] = Bessel_Pro00[Cwan][L0][Mul0][i]*tmp_SphB_D[i + LL*GL_Mesh]; 
              Bes01[i] = Bessel_Pro01[Cwan][L0][Mul0][i]*tmp_SphBp_D[i + LL*GL_Mesh]; 
            }
#else
            #pragma _NEC vector
            for (i=0; i<GL_Mesh; i++){
              Bes00[i] = Bessel_Pro00[Cwan][L0][Mul0][i]*tmp_SphB_D[i*(Lmax_Four_Int+3) + LL]; 
              Bes01[i] = Bessel_Pro01[Cwan][L0][Mul0][i]*tmp_SphBp_D[i*(Lmax_Four_Int+3) + LL]; 
            }
#endif            
#else
            for (i=0; i<GL_Mesh; i++){
              Bes00[i] = Bessel_Pro00[Cwan][L0][Mul0][i]*SphB[LL][i]; 
              Bes01[i] = Bessel_Pro01[Cwan][L0][Mul0][i]*SphBp[LL][i]; 
            }
#endif   

ftrace_region_end("R4.2");
ftrace_region_begin("R4.3");
            L = 0;
            for (L1=0; L1<=List_YOUSO[35]; L1++){ 
              for (Mul1=0; Mul1<List_YOUSO[34]; Mul1++){         

                tmp1 = 0.0;
                tmp2 = 0.0;

                #pragma _NEC vector
                for (i=0; i<GL_Mesh; i++){
                  tmp1 += Bes00[i]*Spe_VNA_Bessel[Hwan][L1][Mul1][i];
                  tmp2 += Bes01[i]*Spe_VNA_Bessel[Hwan][L1][Mul1][i];
            		}
#ifdef _NEC_LIST_V4NL
                const int offset = L1*L1*width_Mul + Mul1*(2*L1+1);
                #pragma _NEC vector
                for (int k=0; k<=2*L1; k++){
                  SumNL0[offset+k]  = tmp1;
                  SumNLr0[offset+k] = tmp2;
                }
#else
                SumNL0[L]  = tmp1;
                SumNLr0[L] = tmp2;

                L++; 
#endif                
	            }
	          }

ftrace_region_end("R4.3");
ftrace_region_begin("R4.4");	  
	          /* derivatives of "on site" */
        
	          if (h_AN==0){ 
	            for (L=0; L<Num_RVNA; L++){
		            SumNLr0[L] = 0.0;
	            }
	          }
          
          
ftrace_region_end("R4.4");
ftrace_region_begin("R4.5");
	/****************************************************
         for "overlap",
	    sum_m 8*(-i)^{-L0+L1+l}*
	          C_{L0,-M0,L1,M1,l,m}*
	          \int RL(k)*RL'(k)*jl(k*R) k^2 dk^3,
	****************************************************/

          	for(m=-LL; m<=LL; m++){ 

//ftrace_region_begin("R4.5.1");
              ComplexSH(LL,m,theta,phi,SH,dSHt,dSHp);
//ftrace_region_end("R4.5.1");
              SH[1]   = -SH[1];
              dSHt[1] = -dSHt[1];
              dSHp[1] = -dSHp[1];

              CY  = Complex(SH[0],SH[1]);
              CYt = Complex(dSHt[0],dSHt[1]);
              CYp = Complex(dSHp[0],dSHp[1]);
#ifdef _NEC_LIST_V4
              #pragma _NEC vector
              for (int iml=0; iml<Num_RVNA_with_M; iml++){
                
              /* note:
                iml = offset + L1 + M1;
                offset = (2*L1+1) * Mul1 + L1*L1*width_Mul;
                L1 + M1 <= 2*L1+1;
                iml <= (2*L1+1) * (Mul1+1) + L1*L1*width_Mul;
                iml/width_Mul >= L1*L1  ->  sqrt(iml/width_Mul) >= L1;
                iml/width_Mul = L1*L1 + (2*L1+1) * (Mul1+1)/width_Mul < L1*L1 + (2*L1+1) = (L1+1)*(L1+1)
                  -> sqrt(iml/width_Mul) < L1+1;
                Consequently, we can get L1 from iml
                L1 = sqrt(iml/width_Mul)
              */   
                const int L1 = (int)sqrt(iml/width_Mul);
                const int M1 = ((iml - L1*L1*width_Mul) % (2*L1+1)) - L1;
                const int L = ((iml - L1*L1*width_Mul) / (2*L1+1)) + L1*width_Mul;
                //const int offset = (2*L1+1) * Mul1 + L1*L1*width_Mul;
                //int M1 = iml - offset - L1;
                int M0 = M1 + m;
                int Ls = -L0 + L1 + LL;
      	        if ( ((L1-LL)<=L0) && ((LL-L1)<=L0) && L0<=(L1+LL) ){
                  if ((M0<=L0) && (-M0<=L0)){
                    Cpow = Im_pow(-1,Ls);
                    CY1  = Cmul(Cpow,CY);
                    CYt1 = Cmul(Cpow,CYt);
                    CYp1 = Cmul(Cpow,CYp);
                  
#ifdef _NEC_MG
                      gant = GetGaunt(L0,M0,L1,M1,LL,m, Gaunt_list, ai_L1_LIMIT_SQ, ai_LL_LIMIT_SQ);
#else
                      gant = Gaunt(L0,M0,L1,M1,LL,m);
#endif
#ifdef _NEC_LIST_V4NL
                      double tmp2 = gant*SumNL0[iml];
                      double tmp0 = gant*SumNLr0[iml];
#else                      
                      double tmp2 = gant*SumNL0[L];
                      double tmp0 = gant*SumNLr0[L];
#endif
                      /* S */ 
                      TmpNL[L0+M0][iml].r += CY1.r*tmp2;
                      TmpNL[L0+M0][iml].i += CY1.i*tmp2;
                      /* dS/dr */ 
                      
                      TmpNLr[L0+M0][iml].r += CY1.r*tmp0;
                      TmpNLr[L0+M0][iml].i += CY1.i*tmp0;
                      /* dS/dt */ 
                      TmpNLt[L0+M0][iml].r += CYt1.r*tmp2;
                      TmpNLt[L0+M0][iml].i += CYt1.i*tmp2;
                      /* dS/dp */ 
                      TmpNLp[L0+M0][iml].r += CYp1.r*tmp2;
                      TmpNLp[L0+M0][iml].i += CYp1.i*tmp2;
                    }
                  }
                }

#else              
	            for (L=0; L<Num_RVNA; L++){
                const int L1 = VNA_List[L];
                const int offset = VNA_offset[L];
                Ls = -L0 + L1 + LL;

      	        if ( ((L1-LL)<=L0) && ((LL-L1)<=L0) && L0<=(L1+LL) ){
                  Cpow = Im_pow(-1,Ls);
                  CY1  = Cmul(Cpow,CY);
                  CYt1 = Cmul(Cpow,CYt);
                  CYp1 = Cmul(Cpow,CYp);

                  for (M0=-L0; M0<=L0; M0++){

                    M1 = M0 - m;

                    if ((M1<=L1) && (-M1<=L1)){

//ftrace_region_begin("R4.5.2");
#ifdef _NEC_MG
                      gant = GetGaunt(L0,M0,L1,M1,LL,m, Gaunt_list, ai_L1_LIMIT_SQ, ai_LL_LIMIT_SQ);
#else
                      gant = Gaunt(L0,M0,L1,M1,LL,m);
#endif
//ftrace_region_end("R4.5.2");                      
                      tmp2 = gant*SumNL0[L];

                      /* S */ 
                      TmpNL[L0+M0][L][L1+M1].r += CY1.r*tmp2;
                      TmpNL[L0+M0][L][L1+M1].i += CY1.i*tmp2;

                      /* dS/dr */ 

                      tmp0 = gant*SumNLr0[L];

                      TmpNLr[L0+M0][L][L1+M1].r += CY1.r*tmp0;
                      TmpNLr[L0+M0][L][L1+M1].i += CY1.i*tmp0;

                      /* dS/dt */ 

                      TmpNLt[L0+M0][L][L1+M1].r += CYt1.r*tmp2;
                      TmpNLt[L0+M0][L][L1+M1].i += CYt1.i*tmp2;

                      /* dS/dp */ 

                      TmpNLp[L0+M0][L][L1+M1].r += CYp1.r*tmp2;
                      TmpNLp[L0+M0][L][L1+M1].i += CYp1.i*tmp2;

                    }
                  }
                }

              }/* L */
#endif
            }/* m */
                    
ftrace_region_end("R4.5");
          } /* LL */
     //   }
     // }
      /****************************************************
                        Complex to Real
      ****************************************************/


ftrace_region_begin("R5");
       
//      num0 = 0;      
  //    for (L0=0; L0<=Spe_MaxL_Basis[Cwan]; L0++){
    //    for (Mul0=0; Mul0<Spe_Num_Basis[Cwan][L0]; Mul0++){
     

#ifdef _NEC_LIST_V4


          for (M0=-L0; M0<=L0; M0++){
            if( 2*L0 < AI_L0_UNROLL_SIZE){
              #pragma _NEC vector                
              for (int iml=0; iml<Num_RVNA_with_M; iml++){


                dcomplex CsumNL0 = Complex(0.0,0.0);
                dcomplex CsumNLr = Complex(0.0,0.0);
                dcomplex CsumNLt = Complex(0.0,0.0);
                dcomplex CsumNLp = Complex(0.0,0.0);

                #pragma _NEC loop_count(AI_L0_UNROLL_SIZE)
                #pragma _NEC unroll(AI_L0_UNROLL_SIZE)
                for (int L0k=0; L0k<AI_L0_UNROLL_SIZE; L0k++){
                if(L0k<=2*L0){
                  Ctmp1 = Conjg(Comp2Real[L0][L0+M0][L0k]);

                  /* S */
                  Ctmp0 = TmpNL[L0k][iml];
                  Ctmp2 = Cmul(Ctmp1,Ctmp0);
                  CsumNL0 = Cadd(CsumNL0,Ctmp2);

                  /* dS/dr */
                  Ctmp0 = TmpNLr[L0k][iml];
                  Ctmp2 = Cmul(Ctmp1,Ctmp0);
                  CsumNLr = Cadd(CsumNLr,Ctmp2);

                  /* dS/dt */
                  Ctmp0 = TmpNLt[L0k][iml];
                  Ctmp2 = Cmul(Ctmp1,Ctmp0);
                  CsumNLt = Cadd(CsumNLt,Ctmp2);

                  /* dS/dp */
                  Ctmp0 = TmpNLp[L0k][iml];
                  Ctmp2 = Cmul(Ctmp1,Ctmp0);
                  CsumNLp = Cadd(CsumNLp,Ctmp2);
                }
                }

                CmatNL0[iml] = CsumNL0;
                CmatNLr[iml] = CsumNLr;
                CmatNLt[iml] = CsumNLt;
                CmatNLp[iml] = CsumNLp;
              }
            }
            else{
              /*printf("AITUNE:L0*2 = %d >= AI_L0_UNROLL_SIZE=%d \n", L0*2, AI_L0_UNROLL_SIZE);*/

              #pragma _NEC vector                
              for (int iml=0; iml<Num_RVNA_with_M; iml++){


                dcomplex CsumNL0 = Complex(0.0,0.0);
                dcomplex CsumNLr = Complex(0.0,0.0);
                dcomplex CsumNLt = Complex(0.0,0.0);
                dcomplex CsumNLp = Complex(0.0,0.0);

                for (int L0k=0; L0k<=2*L0; L0k++){
                  Ctmp1 = Conjg(Comp2Real[L0][L0+M0][L0k]);

                  /* S */
                  Ctmp0 = TmpNL[L0k][iml];
                  Ctmp2 = Cmul(Ctmp1,Ctmp0);
                  CsumNL0 = Cadd(CsumNL0,Ctmp2);

                  /* dS/dr */
                  Ctmp0 = TmpNLr[L0k][iml];
                  Ctmp2 = Cmul(Ctmp1,Ctmp0);
                  CsumNLr = Cadd(CsumNLr,Ctmp2);

                  /* dS/dt */
                  Ctmp0 = TmpNLt[L0k][iml];
                  Ctmp2 = Cmul(Ctmp1,Ctmp0);
                  CsumNLt = Cadd(CsumNLt,Ctmp2);

                  /* dS/dp */
                  Ctmp0 = TmpNLp[L0k][iml];
                  Ctmp2 = Cmul(Ctmp1,Ctmp0);
                  CsumNLp = Cadd(CsumNLp,Ctmp2);
                
                }

                CmatNL0[iml] = CsumNL0;
                CmatNLr[iml] = CsumNLr;
                CmatNLt[iml] = CsumNLt;
                CmatNLp[iml] = CsumNLp;

              }
            }

            if(L1_max*2 < AI_L1_UNROLL_SIZE){
              #pragma _NEC vector         
              for (int iml=0; iml<Num_RVNA_with_M; iml++){

                CsumNL0 = Complex(0.0,0.0);
                CsumNLr = Complex(0.0,0.0);
                CsumNLt = Complex(0.0,0.0);
                CsumNLp = Complex(0.0,0.0);
                

                const int L1 = (int)sqrt(iml/width_Mul);
                const int L1M1 = ((iml - L1*L1*width_Mul) % (2*L1+1));
                //const int Mul1 = ((iml - L1*L1*width_Mul) / (2*L1+1));
                #pragma _NEC loop_count(AI_L1_UNROLL_SIZE)
                #pragma _NEC unroll(AI_L1_UNROLL_SIZE)
                for (int L1k=0; L1k<AI_L1_UNROLL_SIZE; L1k++){
                  if(L1k<=2*L1){
                  /* S */ 

                  CsumNL0.r += (CmatNL0[iml-L1M1+L1k].r*Comp2Real[L1][L1M1][L1k].r
                    - CmatNL0[iml-L1M1+L1k].i*Comp2Real[L1][L1M1][L1k].i);

                  CsumNL0.i += (CmatNL0[iml-L1M1+L1k].r*Comp2Real[L1][L1M1][L1k].i
                    + CmatNL0[iml-L1M1+L1k].i*Comp2Real[L1][L1M1][L1k].r);

                  /* dS/dr */ 

                  CsumNLr.r += (CmatNLr[iml-L1M1+L1k].r*Comp2Real[L1][L1M1][L1k].r
                    - CmatNLr[iml-L1M1+L1k].i*Comp2Real[L1][L1M1][L1k].i);

                  CsumNLr.i += (CmatNLr[iml-L1M1+L1k].r*Comp2Real[L1][L1M1][L1k].i
                    + CmatNLr[iml-L1M1+L1k].i*Comp2Real[L1][L1M1][L1k].r);

                  /* dS/dt */ 

                  CsumNLt.r += (CmatNLt[iml-L1M1+L1k].r*Comp2Real[L1][L1M1][L1k].r
                    - CmatNLt[iml-L1M1+L1k].i*Comp2Real[L1][L1M1][L1k].i);

                  CsumNLt.i += (CmatNLt[iml-L1M1+L1k].r*Comp2Real[L1][L1M1][L1k].i
                    + CmatNLt[iml-L1M1+L1k].i*Comp2Real[L1][L1M1][L1k].r);

                  /* dS/dp */ 

                  CsumNLp.r += (CmatNLp[iml-L1M1+L1k].r*Comp2Real[L1][L1M1][L1k].r
                    - CmatNLp[iml-L1M1+L1k].i*Comp2Real[L1][L1M1][L1k].i);

                  CsumNLp.i += (CmatNLp[iml-L1M1+L1k].r*Comp2Real[L1][L1M1][L1k].i
                    + CmatNLp[iml-L1M1+L1k].i*Comp2Real[L1][L1M1][L1k].r);
                  }
                }



                DS_VNA[0][Mc_AN][h_AN][num0+L0+M0][iml] = 8.0*(Type_DS_VNA)CsumNL0.r;

                if (h_AN!=0){

                  if (fabs(siT)<10e-14){

                    DS_VNA[1][Mc_AN][h_AN][num0+L0+M0][iml] =
                      -8.0*(Type_DS_VNA)(siT*coP*CsumNLr.r + coT*coP/r*CsumNLt.r);

                    DS_VNA[2][Mc_AN][h_AN][num0+L0+M0][iml] =
                      -8.0*(Type_DS_VNA)(siT*siP*CsumNLr.r + coT*siP/r*CsumNLt.r);

                    DS_VNA[3][Mc_AN][h_AN][num0+L0+M0][iml] =
                      -8.0*(Type_DS_VNA)(coT*CsumNLr.r - siT/r*CsumNLt.r);
                  }

                  else{

                    DS_VNA[1][Mc_AN][h_AN][num0+L0+M0][iml] =
                      -8.0*(Type_DS_VNA)(siT*coP*CsumNLr.r + coT*coP/r*CsumNLt.r
                      - siP/siT/r*CsumNLp.r);

                    DS_VNA[2][Mc_AN][h_AN][num0+L0+M0][iml] =
                      -8.0*(Type_DS_VNA)(siT*siP*CsumNLr.r + coT*siP/r*CsumNLt.r
                      + coP/siT/r*CsumNLp.r);

                    DS_VNA[3][Mc_AN][h_AN][num0+L0+M0][iml] =
                      -8.0*(Type_DS_VNA)(coT*CsumNLr.r - siT/r*CsumNLt.r);
                  }

                }

                else{
                  DS_VNA[1][Mc_AN][h_AN][num0+L0+M0][iml] = 0.0;
                  DS_VNA[2][Mc_AN][h_AN][num0+L0+M0][iml] = 0.0;
                  DS_VNA[3][Mc_AN][h_AN][num0+L0+M0][iml] = 0.0;
                } 

              }
            }else{
              /*printf("AITUNE:L1_max*2 = %d >= AI_L1_UNROLL_SIZE=%d \n", L1_max*2, AI_L1_UNROLL_SIZE);*/
              #pragma _NEC vector         
              for (int iml=0; iml<Num_RVNA_with_M; iml++){

                CsumNL0 = Complex(0.0,0.0);
                CsumNLr = Complex(0.0,0.0);
                CsumNLt = Complex(0.0,0.0);
                CsumNLp = Complex(0.0,0.0);
                

                const int L1 = (int)sqrt(iml/width_Mul);
                const int L1M1 = ((iml - L1*L1*width_Mul) % (2*L1+1));
                for (int L1k=0; L1k<=2*L1; L1k++){
                  
                  /* S */ 

                  CsumNL0.r += (CmatNL0[iml-L1M1+L1k].r*Comp2Real[L1][L1M1][L1k].r
                    - CmatNL0[iml-L1M1+L1k].i*Comp2Real[L1][L1M1][L1k].i);

                  CsumNL0.i += (CmatNL0[iml-L1M1+L1k].r*Comp2Real[L1][L1M1][L1k].i
                    + CmatNL0[iml-L1M1+L1k].i*Comp2Real[L1][L1M1][L1k].r);

                  /* dS/dr */ 

                  CsumNLr.r += (CmatNLr[iml-L1M1+L1k].r*Comp2Real[L1][L1M1][L1k].r
                    - CmatNLr[iml-L1M1+L1k].i*Comp2Real[L1][L1M1][L1k].i);

                  CsumNLr.i += (CmatNLr[iml-L1M1+L1k].r*Comp2Real[L1][L1M1][L1k].i
                    + CmatNLr[iml-L1M1+L1k].i*Comp2Real[L1][L1M1][L1k].r);

                  /* dS/dt */ 

                  CsumNLt.r += (CmatNLt[iml-L1M1+L1k].r*Comp2Real[L1][L1M1][L1k].r
                    - CmatNLt[iml-L1M1+L1k].i*Comp2Real[L1][L1M1][L1k].i);

                  CsumNLt.i += (CmatNLt[iml-L1M1+L1k].r*Comp2Real[L1][L1M1][L1k].i
                    + CmatNLt[iml-L1M1+L1k].i*Comp2Real[L1][L1M1][L1k].r);

                  /* dS/dp */ 

                  CsumNLp.r += (CmatNLp[iml-L1M1+L1k].r*Comp2Real[L1][L1M1][L1k].r
                    - CmatNLp[iml-L1M1+L1k].i*Comp2Real[L1][L1M1][L1k].i);

                  CsumNLp.i += (CmatNLp[iml-L1M1+L1k].r*Comp2Real[L1][L1M1][L1k].i
                    + CmatNLp[iml-L1M1+L1k].i*Comp2Real[L1][L1M1][L1k].r);
                  
                }



                DS_VNA[0][Mc_AN][h_AN][num0+L0+M0][iml] = 8.0*(Type_DS_VNA)CsumNL0.r;

                if (h_AN!=0){

                  if (fabs(siT)<10e-14){

                    DS_VNA[1][Mc_AN][h_AN][num0+L0+M0][iml] =
                      -8.0*(Type_DS_VNA)(siT*coP*CsumNLr.r + coT*coP/r*CsumNLt.r);

                    DS_VNA[2][Mc_AN][h_AN][num0+L0+M0][iml] =
                      -8.0*(Type_DS_VNA)(siT*siP*CsumNLr.r + coT*siP/r*CsumNLt.r);

                    DS_VNA[3][Mc_AN][h_AN][num0+L0+M0][iml] =
                      -8.0*(Type_DS_VNA)(coT*CsumNLr.r - siT/r*CsumNLt.r);
                  }

                  else{

                    DS_VNA[1][Mc_AN][h_AN][num0+L0+M0][iml] =
                      -8.0*(Type_DS_VNA)(siT*coP*CsumNLr.r + coT*coP/r*CsumNLt.r
                      - siP/siT/r*CsumNLp.r);

                    DS_VNA[2][Mc_AN][h_AN][num0+L0+M0][iml] =
                      -8.0*(Type_DS_VNA)(siT*siP*CsumNLr.r + coT*siP/r*CsumNLt.r
                      + coP/siT/r*CsumNLp.r);

                    DS_VNA[3][Mc_AN][h_AN][num0+L0+M0][iml] =
                      -8.0*(Type_DS_VNA)(coT*CsumNLr.r - siT/r*CsumNLt.r);
                  }

                }

                else{
                  DS_VNA[1][Mc_AN][h_AN][num0+L0+M0][iml] = 0.0;
                  DS_VNA[2][Mc_AN][h_AN][num0+L0+M0][iml] = 0.0;
                  DS_VNA[3][Mc_AN][h_AN][num0+L0+M0][iml] = 0.0;
                } 

              }
            }
          }
#else

          num1 = 0;
          for (L=0; L<Num_RVNA; L++){
            L1 = VNA_List[L];

            for (M0=-L0; M0<=L0; M0++){
              for (M1=-L1; M1<=L1; M1++){

                dcomplex CsumNL0 = Complex(0.0,0.0);
                dcomplex CsumNLr = Complex(0.0,0.0);
                dcomplex CsumNLt = Complex(0.0,0.0);
                dcomplex CsumNLp = Complex(0.0,0.0);

                for (k=-L0; k<=L0; k++){

                  Ctmp1 = Conjg(Comp2Real[L0][L0+M0][L0+k]);

                  /* S */

                  Ctmp0 = TmpNL[L0+k][L][L1+M1];

                  Ctmp2 = Cmul(Ctmp1,Ctmp0);
                  CsumNL0 = Cadd(CsumNL0,Ctmp2);

                  /* dS/dr */

                  Ctmp0 = TmpNLr[L0+k][L][L1+M1];

                  Ctmp2 = Cmul(Ctmp1,Ctmp0);
                  CsumNLr = Cadd(CsumNLr,Ctmp2);

                  /* dS/dt */


                  Ctmp0 = TmpNLt[L0+k][L][L1+M1];

                  Ctmp2 = Cmul(Ctmp1,Ctmp0);
                  CsumNLt = Cadd(CsumNLt,Ctmp2);

                  /* dS/dp */


                  Ctmp0 = TmpNLp[L0+k][L][L1+M1];

                  Ctmp2 = Cmul(Ctmp1,Ctmp0);
                  CsumNLp = Cadd(CsumNLp,Ctmp2);

                }

                CmatNL0[L1+M1] = CsumNL0;
                CmatNLr[L1+M1] = CsumNLr;
                CmatNLt[L1+M1] = CsumNLt;
                CmatNLp[L1+M1] = CsumNLp;

              }
              

              for (M1=-L1; M1<=L1; M1++){

                CsumNL0 = Complex(0.0,0.0);
                CsumNLr = Complex(0.0,0.0);
                CsumNLt = Complex(0.0,0.0);
                CsumNLp = Complex(0.0,0.0);

                for (k=-L1; k<=L1; k++){

                  /* S */ 

                  CsumNL0.r += (CmatNL0[L1+k].r*Comp2Real[L1][L1+M1][L1+k].r
                    - CmatNL0[L1+k].i*Comp2Real[L1][L1+M1][L1+k].i);

                  CsumNL0.i += (CmatNL0[L1+k].r*Comp2Real[L1][L1+M1][L1+k].i
                    + CmatNL0[L1+k].i*Comp2Real[L1][L1+M1][L1+k].r);

                  /* dS/dr */ 

                  CsumNLr.r += (CmatNLr[L1+k].r*Comp2Real[L1][L1+M1][L1+k].r
                    - CmatNLr[L1+k].i*Comp2Real[L1][L1+M1][L1+k].i);

                  CsumNLr.i += (CmatNLr[L1+k].r*Comp2Real[L1][L1+M1][L1+k].i
                    + CmatNLr[L1+k].i*Comp2Real[L1][L1+M1][L1+k].r);

                  /* dS/dt */ 

                  CsumNLt.r += (CmatNLt[L1+k].r*Comp2Real[L1][L1+M1][L1+k].r
                    - CmatNLt[L1+k].i*Comp2Real[L1][L1+M1][L1+k].i);

                  CsumNLt.i += (CmatNLt[L1+k].r*Comp2Real[L1][L1+M1][L1+k].i
                    + CmatNLt[L1+k].i*Comp2Real[L1][L1+M1][L1+k].r);

                  /* dS/dp */ 

                  CsumNLp.r += (CmatNLp[L1+k].r*Comp2Real[L1][L1+M1][L1+k].r
                    - CmatNLp[L1+k].i*Comp2Real[L1][L1+M1][L1+k].i);

                  CsumNLp.i += (CmatNLp[L1+k].r*Comp2Real[L1][L1+M1][L1+k].i
                    + CmatNLp[L1+k].i*Comp2Real[L1][L1+M1][L1+k].r);

                }

                DS_VNA[0][Mc_AN][h_AN][num0+L0+M0][num1+L1+M1] = 8.0*(Type_DS_VNA)CsumNL0.r;

                if (h_AN!=0){

                  if (fabs(siT)<10e-14){

                    DS_VNA[1][Mc_AN][h_AN][num0+L0+M0][num1+L1+M1] =
                      -8.0*(Type_DS_VNA)(siT*coP*CsumNLr.r + coT*coP/r*CsumNLt.r);

                    DS_VNA[2][Mc_AN][h_AN][num0+L0+M0][num1+L1+M1] =
                      -8.0*(Type_DS_VNA)(siT*siP*CsumNLr.r + coT*siP/r*CsumNLt.r);

                    DS_VNA[3][Mc_AN][h_AN][num0+L0+M0][num1+L1+M1] =
                      -8.0*(Type_DS_VNA)(coT*CsumNLr.r - siT/r*CsumNLt.r);
                  }

                  else{

                    DS_VNA[1][Mc_AN][h_AN][num0+L0+M0][num1+L1+M1] =
                      -8.0*(Type_DS_VNA)(siT*coP*CsumNLr.r + coT*coP/r*CsumNLt.r
                      - siP/siT/r*CsumNLp.r);

                    DS_VNA[2][Mc_AN][h_AN][num0+L0+M0][num1+L1+M1] =
                      -8.0*(Type_DS_VNA)(siT*siP*CsumNLr.r + coT*siP/r*CsumNLt.r
                      + coP/siT/r*CsumNLp.r);

                    DS_VNA[3][Mc_AN][h_AN][num0+L0+M0][num1+L1+M1] =
                      -8.0*(Type_DS_VNA)(coT*CsumNLr.r - siT/r*CsumNLt.r);
                  }

                }

                else{
                  DS_VNA[1][Mc_AN][h_AN][num0+L0+M0][num1+L1+M1] = 0.0;
                  DS_VNA[2][Mc_AN][h_AN][num0+L0+M0][num1+L1+M1] = 0.0;
                  DS_VNA[3][Mc_AN][h_AN][num0+L0+M0][num1+L1+M1] = 0.0;
                } 

              }
            }
            num1 = num1 + 2*L1 + 1; 
          }/* L, that is, L1 */
#endif

ftrace_region_end("R5");
          num0 = num0 + 2*L0 + 1; 
        } /* Mul0 */
      } /* L0 */

      


      dtime(&Etime_atom);
      time_per_atom[Gc_AN] += Etime_atom - Stime_atom;

    } /* Mc_AN */

    /* freeing of arrays */
#ifdef _NEC3
  	free(tmp_SphB_D);
	  free(tmp_SphBp_D);
#endif

    free(CmatNL0);
    free(CmatNLr);
    free(CmatNLt);
    free(CmatNLp);

    free(Bes00);
    free(Bes01);

    free(SumNLr0);

    free(SumNL0);

#ifdef _NEC_LIST_V4
  free(TmpNLp[0]);
  free(TmpNLp);
  free(TmpNLr[0]);
  free(TmpNLr);
  free(TmpNLt[0]);
  free(TmpNLt);
  free(TmpNL[0]);
  free(TmpNL);
#else
	for (k=0; k<(2*(List_YOUSO[25]+1)+1); k++){
	  for (l=0; l<Num_RVNA; l++){
	    free(TmpNLp[k][l]);
	  }
          free(TmpNLp[k]);
	}
    free(TmpNLp);

	for (k=0; k<(2*(List_YOUSO[25]+1)+1); k++){
	  for (l=0; l<Num_RVNA; l++){
	    free(TmpNLt[k][l]);
	  }
	  free(TmpNLt[k]);
	}
    free(TmpNLt);

	for (k=0; k<(2*(List_YOUSO[25]+1)+1); k++){
	  for (l=0; l<Num_RVNA; l++){
	    free(TmpNLr[k][l]);
	  }
	  free(TmpNLr[k]);
	}
	free(TmpNLr);

  for (k=0; k<(2*(List_YOUSO[25]+1)+1); k++){
	  for (l=0; l<Num_RVNA; l++){
	    free(TmpNL[k][l]);
	  }
	  free(TmpNL[k]);
	}
	free(TmpNL);
#endif
#pragma omp flush(DS_VNA)

  } /* #pragma omp parallel */
  
  if (measure_time){
	  dtime(&etime);
	  time3 += etime - stime;
  }

  /*******************************************************
   *******************************************************
     multiplying overlap integrals WITH COMMUNICATION

     MPI: communicate only for k=0
     DS_VNA
  *******************************************************
  *******************************************************/

  MPI_Barrier(mpi_comm_level1);
  if (measure_time) dtime(&stime);

  /* allocation of array */

  NLH = (double**)malloc(sizeof(double*)*List_YOUSO[7]); 
  for (i=0; i<List_YOUSO[7]; i++){
    NLH[i] = (double*)malloc(sizeof(double)*List_YOUSO[7]); 
  }

  for (ID=0; ID<numprocs; ID++){
    F_Snd_Num_WK[ID] = 0;
    F_Rcv_Num_WK[ID] = 0;
  }

  do {

    /***********************************
          set the size of data
    ************************************/
ftrace_region_begin("DO-1");
    for (ID=0; ID<numprocs; ID++){

      IDS = (myid + ID) % numprocs;
      IDR = (myid - ID + numprocs) % numprocs;

      /* find the data size to send the block data */

      if ( 0<(F_Snd_Num[IDS]-F_Snd_Num_WK[IDS]) ){

	size1 = 0;
	n = F_Snd_Num_WK[IDS];

	Mc_AN = Snd_MAN[IDS][n];
	Gc_AN = Snd_GAN[IDS][n];
	Cwan = WhatSpecies[Gc_AN]; 
	tno1 = Spe_Total_NO[Cwan];
	for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
	  Gh_AN = natn[Gc_AN][h_AN];        
	  Hwan = WhatSpecies[Gh_AN];
	  tno2 = (List_YOUSO[35]+1)*(List_YOUSO[35]+1)*List_YOUSO[34];
	  size1 += tno1*tno2; 
	}
 
	Snd_DS_VNA_Size[IDS] = size1;
	MPI_Isend(&size1, 1, MPI_INT, IDS, tag, mpi_comm_level1, &request);
      }
      else{
	Snd_DS_VNA_Size[IDS] = 0;
      }

      /* receiving of the size of the data */

      if ( 0<(F_Rcv_Num[IDR]-F_Rcv_Num_WK[IDR]) ){
	MPI_Recv(&size2, 1, MPI_INT, IDR, tag, mpi_comm_level1, &stat);
	Rcv_DS_VNA_Size[IDR] = size2;
      }
      else{
	Rcv_DS_VNA_Size[IDR] = 0;
      }

      if ( 0<(F_Snd_Num[IDS]-F_Snd_Num_WK[IDS]) )  MPI_Wait(&request,&stat);

    } /* ID */
ftrace_region_end("DO-1");
    /***********************************
               data transfer
    ************************************/

    for (ID=0; ID<numprocs; ID++){

      IDS = (myid + ID) % numprocs;
      IDR = (myid - ID + numprocs) % numprocs;

      /******************************
         sending of the data 
      ******************************/
ftrace_region_begin("DO-2-1");
      if ( 0<(F_Snd_Num[IDS]-F_Snd_Num_WK[IDS]) ){

	size1 = Snd_DS_VNA_Size[IDS];

	/*
	printf("S myid=%2d ID=%2d IDS=%2d IDR=%2d  size1=%2d\n",myid,ID,IDS,IDR,size1);fflush(stdout);
	*/

	/* allocation of the array */

	tmp_array = (Type_DS_VNA*)malloc(sizeof(Type_DS_VNA)*size1);

	/* multidimentional array to the vector array */

	num = 0;
	n = F_Snd_Num_WK[IDS];

	Mc_AN = Snd_MAN[IDS][n];
	Gc_AN = Snd_GAN[IDS][n];
	Cwan = WhatSpecies[Gc_AN]; 
	tno1 = Spe_Total_NO[Cwan];
	for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
	  Gh_AN = natn[Gc_AN][h_AN];        
	  Hwan = WhatSpecies[Gh_AN];
	  tno2 = (List_YOUSO[35]+1)*(List_YOUSO[35]+1)*List_YOUSO[34];

	  for (i=0; i<tno1; i++){
	    for (j=0; j<tno2; j++){
	      tmp_array[num] = DS_VNA[0][Mc_AN][h_AN][i][j];
	      num++;
	    } 
	  } 
	}

	MPI_Isend(&tmp_array[0], size1, MPI_Type_DS_VNA, IDS, tag, mpi_comm_level1, &request);

      }
ftrace_region_end("DO-2-1");
      /******************************
        receiving of the block data
      ******************************/

      if ( 0<(F_Rcv_Num[IDR]-F_Rcv_Num_WK[IDR]) ){
ftrace_region_begin("DO-2-2");  //second worst//
	size2 = Rcv_DS_VNA_Size[IDR];

	/*
	printf("R myid=%2d ID=%2d IDS=%2d IDR=%2d  size2=%2d\n",myid,ID,IDS,IDR,size2);fflush(stdout);
	*/

	tmp_array2 = (Type_DS_VNA*)malloc(sizeof(Type_DS_VNA)*size2);
	MPI_Recv(&tmp_array2[0], size2, MPI_Type_DS_VNA, IDR, tag, mpi_comm_level1, &stat);

ftrace_region_end("DO-2-2");
ftrace_region_begin("DO-2-3");
	/* store */

	num = 0;
	n = F_Rcv_Num_WK[IDR];
	Original_Mc_AN = F_TopMAN[IDR] + n;

	Gc_AN = Rcv_GAN[IDR][n];
	Cwan = WhatSpecies[Gc_AN];
	tno1 = Spe_Total_NO[Cwan];
	for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
	  Gh_AN = natn[Gc_AN][h_AN];
	  Hwan = WhatSpecies[Gh_AN];
	  tno2 = (List_YOUSO[35]+1)*(List_YOUSO[35]+1)*List_YOUSO[34];
	  for (i=0; i<tno1; i++){
	    for (j=0; j<tno2; j++){
	      DS_VNA[0][Matomnum+1][h_AN][i][j] = tmp_array2[num];
	      num++;
	    }
	  }
	}

	/* free tmp_array2 */
	free(tmp_array2);

ftrace_region_end("DO-2-3");
ftrace_region_begin("DO-2-4"); //worst//
	/*****************************************
              multiplying overlap integrals
	*****************************************/

	for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){

	  dtime(&Stime_atom);
    
	  Gc_AN = M2G[Mc_AN];
	  Cwan = WhatSpecies[Gc_AN];
	  fan = FNAN[Gc_AN];
	  rcutA = Spe_Atom_Cut1[Cwan];

	  n = F_Rcv_Num_WK[IDR];
	  jg = Rcv_GAN[IDR][n];

	  for (j0=0; j0<=fan; j0++){

	    jg0 = natn[Gc_AN][j0];
	    Mj_AN0 = F_G2M[jg0];

  	    po = 0;
	    if (Original_Mc_AN==Mj_AN0){
	      po = 1;
	      j = j0;
#ifndef _NEC        
	    }
  
	    if (po==1){
#endif
	      Hwan = WhatSpecies[jg];
	      rcutB = Spe_Atom_Cut1[Hwan];
	      rcut = rcutA + rcutB;

	      for (k=0; k<=fan; k++){

		kg = natn[Gc_AN][k];
		wakg = WhatSpecies[kg];
		kl = RMI1[Mc_AN][j][k];

		if (0<=kl){

//#pragma omp parallel for private(m, n, sum, L ,L1, GL, Mul1, ene, L2 ,tmp0, L3) default(shared) if(Spe_Total_NO[Cwan] > 8)
		  for (int m=0; m<Spe_Total_NO[Cwan]; m++){
        #pragma _NEC novector
		    for (int n=0; n<Spe_Total_NO[Hwan]; n++){

		      double sum = 0.0;



#ifdef _NEC_LIST_V4
      #pragma _NEC ivdep
      #pragma _NEC vector
		  for (int iml=0; iml<Num_RVNA_with_M; iml++){
        const int GL = (int)sqrt(0.01 + iml/width_Mul);        
        const int Mul1 = ((iml - GL*GL*width_Mul) / (2*GL+1));
		    double ene = VNA_proj_ene[wakg][GL][Mul1];
		    sum += ene*(DS_VNA[0][Mc_AN][k][m][iml]*DS_VNA[0][Matomnum+1][kl][n][iml]); /*bracket is important because DS_VNA is float (not double))*/
      }
#else

		      L = 0;
          

		    for (int L1=0; L1<Num_RVNA; L1++){    

			const int GL = VNA_List[L1];
			const int Mul1 = VNA_List2[L1];
			double ene = VNA_proj_ene[wakg][GL][Mul1];
			const int L2 = 2*VNA_List[L1];
      const int offset = VNA_offset[L1]; //offset of L that is sum of (L1*2+1)//

			double tmp0 = 0.0;

			for (L3=0; L3<=L2; L3++){
			  tmp0 += DS_VNA[0][Mc_AN][k][m][L]*DS_VNA[0][Matomnum+1][kl][n][L];
			  L++;
			}

			sum += ene*tmp0; 

		      }
#endif
		      if (k==0)  NLH[m][n]  = sum;  
		      else       NLH[m][n] += sum; 

		    } /* n */
		  } /* m */

		} /* if (0<=kl)*/
	      } /* k */ 
       
	      /****************************************************
                            NLH to HVNA
	      ****************************************************/

	      dmp = dampingF(rcut,Dis[Gc_AN][j]);

	      for (i1=0; i1<Spe_Total_NO[Cwan]; i1++){
		for (j1=0; j1<Spe_Total_NO[Hwan]; j1++){
		  HVNA[Mc_AN][j][i1][j1] = dmp*NLH[i1][j1];
		}
	      }

	    } /* if (po==1) */

	  } /* j */            

	  dtime(&Etime_atom);
	  time_per_atom[Gc_AN] += Etime_atom - Stime_atom;

	} /* Mc_AN */

	/********************************************
            increment of F_Rcv_Num_WK[IDR] 
	********************************************/

	F_Rcv_Num_WK[IDR]++;

ftrace_region_end("DO-2-4");
      }
ftrace_region_begin("DO-2-5");

      if ( 0<(F_Snd_Num[IDS]-F_Snd_Num_WK[IDS]) ) {
	MPI_Wait(&request,&stat);
	free(tmp_array);  /* freeing of array */

	/********************************************
             increment of F_Snd_Num_WK[IDS]
	********************************************/

	F_Snd_Num_WK[IDS]++;
      } 
ftrace_region_end("DO-2-5");
    } /* ID */
    /*****************************************************
      check whether all the communications have finished
    *****************************************************/

ftrace_region_begin("DO-3");
    po = 0;
    for (ID=0; ID<numprocs; ID++){

      IDS = (myid + ID) % numprocs;
      IDR = (myid - ID + numprocs) % numprocs;

      if ( 0<(F_Snd_Num[IDS]-F_Snd_Num_WK[IDS]) ) po += F_Snd_Num[IDS]-F_Snd_Num_WK[IDS];
      if ( 0<(F_Rcv_Num[IDR]-F_Rcv_Num_WK[IDR]) ) po += F_Rcv_Num[IDR]-F_Rcv_Num_WK[IDR];
    }
  
ftrace_region_end("DO-3");
  } while (po!=0);
  
  /* freeing of arrays */
  
  for (i=0; i<List_YOUSO[7]; i++){
    free(NLH[i]);
  }
  free(NLH);
  
  if (measure_time){
    dtime(&etime);
    time4 += etime - stime;
  }
  /*******************************************************
   *******************************************************
    multiplying overlap integrals WITHOUT COMMUNICATION
  *******************************************************
  *******************************************************/
  if (measure_time) dtime(&stime);
  
#pragma omp parallel shared(List_YOUSO,time_per_atom,HVNA,Dis,DS_VNA,VNA_proj_ene,VNA_List2,VNA_List,Num_RVNA,Spe_Total_NO,RMI1,F_G2M,natn,Spe_Atom_Cut1,FNAN,WhatSpecies,M2G,Matomnum) 
  {         
  
    int OMPID,Nthrds,Nprocs;
    int Mc_AN,Gc_AN,Cwan,fan,Mj_AN,Hwan;
    int L,L1,GL,Mul1,L2,L3,i1,j1;
    int k,kg,wakg,kl,i,j,jg,m,n;
    double **NLH;
    double rcutA,rcutB,rcut,sum,ene,tmp0,dmp;
    double Stime_atom,Etime_atom;
#ifdef _NEC
    double tmp01[120];   /*NUM_RVNA=120*/
#endif

    /* allocation of array */

    NLH = (double**)malloc(sizeof(double*)*List_YOUSO[7]); 
    for (i=0; i<List_YOUSO[7]; i++){
      NLH[i] = (double*)malloc(sizeof(double)*List_YOUSO[7]); 
    }

    /* get info. on OpenMP */ 

    OMPID = omp_get_thread_num();
    Nthrds = omp_get_num_threads();
    Nprocs = omp_get_num_procs();

    for (Mc_AN=(OMPID*Matomnum/Nthrds+1); Mc_AN<((OMPID+1)*Matomnum/Nthrds+1); Mc_AN++){

      dtime(&Stime_atom);
    
      Gc_AN = M2G[Mc_AN];
      Cwan = WhatSpecies[Gc_AN];
      fan = FNAN[Gc_AN];
      rcutA = Spe_Atom_Cut1[Cwan];
  
      for (j=0; j<=fan; j++){
  
	jg = natn[Gc_AN][j];
	Mj_AN = F_G2M[jg];

	if (Mj_AN<=Matomnum){

	  Hwan = WhatSpecies[jg];
	  rcutB = Spe_Atom_Cut1[Hwan];
	  rcut = rcutA + rcutB;

	  for (k=0; k<=fan; k++){

	    kg = natn[Gc_AN][k];
	    wakg = WhatSpecies[kg];
	    kl = RMI1[Mc_AN][j][k];

	    if (0<=kl){


	    for (m=0; m<Spe_Total_NO[Cwan]; m++){
        #pragma _NEC novector
		    for (n=0; n<Spe_Total_NO[Hwan]; n++){

		  sum = 0.0;

#ifdef _NEC_LIST_V4
      #pragma _NEC ivdep
      #pragma _NEC vector
		  for (int iml=0; iml<Num_RVNA_with_M; iml++){
        const int GL = (int)sqrt(iml/width_Mul);        
        const int Mul1 = ((iml - GL*GL*width_Mul) / (2*GL+1));
		    double ene = VNA_proj_ene[wakg][GL][Mul1];
		    
        sum += ene* (DS_VNA[0][Mc_AN][k][m][iml]*DS_VNA[0][Mj_AN][kl][n][iml]);/*bracket is important because DS_VNA is float (not double))*/ 
		    
		  }
#else      
		  L = 0;

      #pragma _NEC ivdep
      #pragma _NEC vector
		  for (int L1=0; L1<Num_RVNA; L1++){

		    int GL = VNA_List[L1]; 
		    int Mul1 = VNA_List2[L1];
		    double ene = VNA_proj_ene[wakg][GL][Mul1];
		    int L2 = 2*VNA_List[L1];
    

		    double tmp0 = 0.0;
		    for (int L3=0; L3<=L2; L3++){
		      tmp0 += DS_VNA[0][Mc_AN][k][m][L]*DS_VNA[0][Mj_AN][kl][n][L]; 
          L++;
        }
		    
		    sum += ene*tmp0; 
		  }
#endif

		  if (k==0)  NLH[m][n]  = sum;  
		  else       NLH[m][n] += sum; 

    
	      }
        }
	    }
	  } /* k */ 
       
	  /****************************************************
                            NLH to HVNA
	  ****************************************************/

	  dmp = dampingF(rcut,Dis[Gc_AN][j]);

	  for (i1=0; i1<Spe_Total_NO[Cwan]; i1++){
	    for (j1=0; j1<Spe_Total_NO[Hwan]; j1++){
	      HVNA[Mc_AN][j][i1][j1] = dmp*NLH[i1][j1];
	    }
	  }

	} /* if (Mj_AN<=Matomnum) */
      } /* j */

      dtime(&Etime_atom);
      time_per_atom[Gc_AN] += Etime_atom - Stime_atom;

    } /* Mc_AN */

    /* freeing of arrays */
 
    for (i=0; i<List_YOUSO[7]; i++){
      free(NLH[i]);
    }
    free(NLH);

#pragma omp flush(HVNA)

  } /* #pragma omp parallel */
  if (measure_time){
    dtime(&etime);
    time5 += etime - stime;
  }

  if (measure_time){
    printf("Set_ProExpn_VNA myid=%2d time1=%7.3f time2=%7.3f time3=%7.3f time4=%7.3f time5=%7.3f\n",
            myid,time1,time2,time3,time4,time5);fflush(stdout); 
  }

  /****************************************************
                   freeing of arrays:
  ****************************************************/

  free(Snd_DS_VNA_Size);
  free(Rcv_DS_VNA_Size);
#ifndef _NEC_LIST_V4
  free(VNA_List);
  free(VNA_List2);
  free(VNA_offset);
#endif
  for (spe=0; spe<SpeciesNum; spe++){
    for (L0=0; L0<=Spe_MaxL_Basis[spe]; L0++){
      for (Mul0=0; Mul0<Spe_Num_Basis[spe][L0]; Mul0++){
        free(Bessel_Pro00[spe][L0][Mul0]);
      }
      free(Bessel_Pro00[spe][L0]);
    }
    free(Bessel_Pro00[spe]);
  }
  free(Bessel_Pro00);

  for (spe=0; spe<SpeciesNum; spe++){
    for (L0=0; L0<=Spe_MaxL_Basis[spe]; L0++){
      for (Mul0=0; Mul0<Spe_Num_Basis[spe][L0]; Mul0++){
        free(Bessel_Pro01[spe][L0][Mul0]);
      }
      free(Bessel_Pro01[spe][L0]);
    }
    free(Bessel_Pro01[spe]);
  }
  free(Bessel_Pro01);

  free(OneD2Mc_AN);
  free(OneD2h_AN);
  free(Gaunt_list);

  /* for time */
  dtime(&TEtime);
  time0 = TEtime - TStime;

  return time0;
} 




double Set_VNA2(double ****HVNA, double *****HVNA2)
{
  /****************************************************
   Evaluate matrix elements of one-center (orbitals)
   but two-center integrals for neutral atom potentials
   in the momentum space.

   <Phi_{LM,L'M'}|VNA>
  ****************************************************/
  static int firsttime=1;
  int L2,L3,L,GL,i,j;
  int k,kl,h_AN,Gh_AN,Rnh,Ls,n;
  int tno0,tno1,tno2,i1,j1,p;
  int Cwan,Hwan,fan,jg,kg,wakg,Lmax;
  int size_TmpHNA;
  int Mc_AN,Gc_AN,Mj_AN,num,size1,size2;
  double time0;
  double tmp2,tmp3,tmp4,tmp5;
  double TStime,TEtime;
  int Num_RVNA;
  /* for OpenMP */
  int OneD_Nloop,*OneD2Mc_AN,*OneD2h_AN;

  dtime(&TStime);

  /****************************************************
  allocation of arrays:
  ****************************************************/

  size_TmpHNA = (List_YOUSO[25]+1)*List_YOUSO[24]*
                (2*(List_YOUSO[25]+1)+1)*
                (List_YOUSO[25]+1)*List_YOUSO[24]*
                (2*(List_YOUSO[25]+1)+1);

  /* PrintMemory */
  if (firsttime) {
    PrintMemory("Set_ProExpn_VNA: TmpHNA",sizeof(double)*size_TmpHNA,NULL);
    PrintMemory("Set_ProExpn_VNA: TmpHNAr",sizeof(double)*size_TmpHNA,NULL);
    PrintMemory("Set_ProExpn_VNA: TmpHNAt",sizeof(double)*size_TmpHNA,NULL);
    PrintMemory("Set_ProExpn_VNA: TmpHNAp",sizeof(double)*size_TmpHNA,NULL);

    firsttime=0;
  }

  /* one-dimensionalize the Mc_AN and h_AN loops */

  OneD_Nloop = 0;
  for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
    Gc_AN = M2G[Mc_AN];    
    for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
      OneD_Nloop++;
    }
  }  

  OneD2Mc_AN = (int*)malloc(sizeof(int)*(OneD_Nloop+1));
  OneD2h_AN = (int*)malloc(sizeof(int)*(OneD_Nloop+1));

  OneD_Nloop = 0;
  for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
    Gc_AN = M2G[Mc_AN];    
    for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
      OneD2Mc_AN[OneD_Nloop] = Mc_AN; 
      OneD2h_AN[OneD_Nloop] = h_AN; 
      OneD_Nloop++;
    }
  }

ftrace_region_begin("MG2");

  int max_L0 = 0;
  for (int spe=0; spe<SpeciesNum; spe++){
    if(max_L0 < Spe_MaxL_Basis[spe]) max_L0 = Spe_MaxL_Basis[spe];
  }
  const int ai_L0_LIMIT = max_L0 + 1;
  const int ai_L0_LIMIT_SQ = ai_L0_LIMIT * ai_L0_LIMIT;
  const int ai_L1_LIMIT = ai_L0_LIMIT;
  const int ai_L1_LIMIT_SQ = ai_L0_LIMIT_SQ;
  const int ai_LL_LIMIT = max_L0*2 + 1;
  const int ai_LL_LIMIT_SQ = ai_LL_LIMIT*ai_LL_LIMIT;
  double* Gaunt_list = MakeGaunt(ai_L0_LIMIT_SQ, ai_L1_LIMIT_SQ, ai_LL_LIMIT_SQ);
ftrace_region_end("MG2");


#pragma omp parallel shared(List_YOUSO,time_per_atom,HVNA2,Comp2Real,GL_Weight,Spe_ProductRF_Bessel,Spe_CrudeVNA_Bessel,GL_NormK,Spe_Num_Basis,Spe_MaxL_Basis,PAO_Nkmax,atv,Gxyz,WhatSpecies,ncn,natn,M2G,OneD2h_AN,OneD2Mc_AN,OneD_Nloop) 
  {         
    int OMPID,Nthrds,Nprocs,Nloop;
    int Mc_AN,h_AN,Gc_AN,Cwan,Gh_AN;
    int Rnh,Hwan,L0,Mul0,L1,Mul1,M0,M1,LL;
    int Lmax_Four_Int,i,j,k,l,m,num0,num1;
    double dx,dy,dz;
    double siT,coT,siP,coP;
    double Normk,kmin,kmax,Sk,Dk;
    double gant,r,theta,phi;
    double SH[2],dSHt[2],dSHp[2];
    double S_coordinate[3];
    double Stime_atom,Etime_atom;
    double sum,sumr,sj,sjp,tmp0,tmp1,tmp10;
#ifndef _NEC2
    double **SphB,**SphBp;
    double *tmp_SphB,*tmp_SphBp;
#endif
    dcomplex Ctmp0,Ctmp1,Ctmp2,Cpow;
    dcomplex Csum,Csumr,Csumt,Csump;
    dcomplex CsumS0,CsumSr,CsumSt,CsumSp;
    dcomplex CY,CYt,CYp,CY1,CYt1,CYp1;
    dcomplex **CmatS0;
    dcomplex **CmatSr;
    dcomplex **CmatSt;
    dcomplex **CmatSp;
    dcomplex ******TmpHNA;
    dcomplex ******TmpHNAr;
    dcomplex ******TmpHNAt;
    dcomplex ******TmpHNAp;

    /* allocation of arrays */

    TmpHNA = (dcomplex******)malloc(sizeof(dcomplex*****)*(List_YOUSO[25]+1));
    for (i=0; i<(List_YOUSO[25]+1); i++){
      TmpHNA[i] = (dcomplex*****)malloc(sizeof(dcomplex****)*List_YOUSO[24]);
      for (j=0; j<List_YOUSO[24]; j++){
	TmpHNA[i][j] = (dcomplex****)malloc(sizeof(dcomplex***)*(2*(List_YOUSO[25]+1)+1));
	for (k=0; k<(2*(List_YOUSO[25]+1)+1); k++){
	  TmpHNA[i][j][k] = (dcomplex***)malloc(sizeof(dcomplex**)*(List_YOUSO[25]+1));
	  for (l=0; l<(List_YOUSO[25]+1); l++){
	    TmpHNA[i][j][k][l] = (dcomplex**)malloc(sizeof(dcomplex*)*List_YOUSO[24]);
	    for (m=0; m<List_YOUSO[24]; m++){
	      TmpHNA[i][j][k][l][m] = (dcomplex*)malloc(sizeof(dcomplex)*(2*(List_YOUSO[25]+1)+1));
	    }
	  }
	}
      }
    }

    TmpHNAr = (dcomplex******)malloc(sizeof(dcomplex*****)*(List_YOUSO[25]+1));
    for (i=0; i<(List_YOUSO[25]+1); i++){
      TmpHNAr[i] = (dcomplex*****)malloc(sizeof(dcomplex****)*List_YOUSO[24]);
      for (j=0; j<List_YOUSO[24]; j++){
	TmpHNAr[i][j] = (dcomplex****)malloc(sizeof(dcomplex***)*(2*(List_YOUSO[25]+1)+1));
	for (k=0; k<(2*(List_YOUSO[25]+1)+1); k++){
	  TmpHNAr[i][j][k] = (dcomplex***)malloc(sizeof(dcomplex**)*(List_YOUSO[25]+1));
	  for (l=0; l<(List_YOUSO[25]+1); l++){
	    TmpHNAr[i][j][k][l] = (dcomplex**)malloc(sizeof(dcomplex*)*List_YOUSO[24]);
	    for (m=0; m<List_YOUSO[24]; m++){
	      TmpHNAr[i][j][k][l][m] = (dcomplex*)malloc(sizeof(dcomplex)*(2*(List_YOUSO[25]+1)+1));
	    }
	  }
	}
      }
    }

    TmpHNAt = (dcomplex******)malloc(sizeof(dcomplex*****)*(List_YOUSO[25]+1));
    for (i=0; i<(List_YOUSO[25]+1); i++){
      TmpHNAt[i] = (dcomplex*****)malloc(sizeof(dcomplex****)*List_YOUSO[24]);
      for (j=0; j<List_YOUSO[24]; j++){
	TmpHNAt[i][j] = (dcomplex****)malloc(sizeof(dcomplex***)*(2*(List_YOUSO[25]+1)+1));
	for (k=0; k<(2*(List_YOUSO[25]+1)+1); k++){
	  TmpHNAt[i][j][k] = (dcomplex***)malloc(sizeof(dcomplex**)*(List_YOUSO[25]+1));
	  for (l=0; l<(List_YOUSO[25]+1); l++){
	    TmpHNAt[i][j][k][l] = (dcomplex**)malloc(sizeof(dcomplex*)*List_YOUSO[24]);
	    for (m=0; m<List_YOUSO[24]; m++){
	      TmpHNAt[i][j][k][l][m] = (dcomplex*)malloc(sizeof(dcomplex)*(2*(List_YOUSO[25]+1)+1));
	    }
	  }
	}
      }
    }

    TmpHNAp = (dcomplex******)malloc(sizeof(dcomplex*****)*(List_YOUSO[25]+1));
    for (i=0; i<(List_YOUSO[25]+1); i++){
      TmpHNAp[i] = (dcomplex*****)malloc(sizeof(dcomplex****)*List_YOUSO[24]);
      for (j=0; j<List_YOUSO[24]; j++){
	TmpHNAp[i][j] = (dcomplex****)malloc(sizeof(dcomplex***)*(2*(List_YOUSO[25]+1)+1));
	for (k=0; k<(2*(List_YOUSO[25]+1)+1); k++){
	  TmpHNAp[i][j][k] = (dcomplex***)malloc(sizeof(dcomplex**)*(List_YOUSO[25]+1));
	  for (l=0; l<(List_YOUSO[25]+1); l++){
	    TmpHNAp[i][j][k][l] = (dcomplex**)malloc(sizeof(dcomplex*)*List_YOUSO[24]);
	    for (m=0; m<List_YOUSO[24]; m++){
	      TmpHNAp[i][j][k][l][m] = (dcomplex*)malloc(sizeof(dcomplex)*(2*(List_YOUSO[25]+1)+1));
	    }
	  }
	}
      }
    }

    CmatS0 = (dcomplex**)malloc(sizeof(dcomplex*)*(2*(List_YOUSO[25]+1)+1));
    for (i=0; i<(2*(List_YOUSO[25]+1)+1); i++){
      CmatS0[i] = (dcomplex*)malloc(sizeof(dcomplex)*(2*(List_YOUSO[25]+1)+1));
    }

    CmatSr = (dcomplex**)malloc(sizeof(dcomplex*)*(2*(List_YOUSO[25]+1)+1));
    for (i=0; i<(2*(List_YOUSO[25]+1)+1); i++){
      CmatSr[i] = (dcomplex*)malloc(sizeof(dcomplex)*(2*(List_YOUSO[25]+1)+1));
    }

    CmatSt = (dcomplex**)malloc(sizeof(dcomplex*)*(2*(List_YOUSO[25]+1)+1));
    for (i=0; i<(2*(List_YOUSO[25]+1)+1); i++){
      CmatSt[i] = (dcomplex*)malloc(sizeof(dcomplex)*(2*(List_YOUSO[25]+1)+1));
    }

    CmatSp = (dcomplex**)malloc(sizeof(dcomplex*)*(2*(List_YOUSO[25]+1)+1));
    for (i=0; i<(2*(List_YOUSO[25]+1)+1); i++){
      CmatSp[i] = (dcomplex*)malloc(sizeof(dcomplex)*(2*(List_YOUSO[25]+1)+1));
    }

#ifdef _NEC2
    double* tmp_SphB_D  = (double*)malloc(sizeof(double)*30*GL_Mesh);
    double* tmp_SphBp_D = (double*)malloc(sizeof(double)*30*GL_Mesh);
#endif

    /* get info. on OpenMP */ 

    OMPID = omp_get_thread_num();
    Nthrds = omp_get_num_threads();
    Nprocs = omp_get_num_procs();

    /* one-dimensionalized loop */

    for (Nloop=OMPID*OneD_Nloop/Nthrds; Nloop<(OMPID+1)*OneD_Nloop/Nthrds; Nloop++){

      dtime(&Stime_atom);

      /* get Mc_AN and h_AN */

      Mc_AN = OneD2Mc_AN[Nloop];
      h_AN  = OneD2h_AN[Nloop];

      /* set data on Mc_AN */

      Gc_AN = M2G[Mc_AN];
      Cwan = WhatSpecies[Gc_AN];

      /* set data on h_AN */

      Gh_AN = natn[Gc_AN][h_AN];        
      Rnh = ncn[Gc_AN][h_AN];
      Hwan = WhatSpecies[Gh_AN];

      dx = Gxyz[Gh_AN][1] + atv[Rnh][1] - Gxyz[Gc_AN][1]; 
      dy = Gxyz[Gh_AN][2] + atv[Rnh][2] - Gxyz[Gc_AN][2]; 
      dz = Gxyz[Gh_AN][3] + atv[Rnh][3] - Gxyz[Gc_AN][3];

      xyz2spherical(dx,dy,dz,0.0,0.0,0.0,S_coordinate); 
      r     = S_coordinate[0];
      theta = S_coordinate[1];
      phi   = S_coordinate[2];

      /* for empty atoms or finite elemens basis */

      if (r<1.0e-10) r = 1.0e-10;

      /* precalculation of sin and cos */

      siT = sin(theta);
      coT = cos(theta);
      siP = sin(phi);
      coP = cos(phi);

      /****************************************************
         evaluate ovelap integrals <Phi_{LM,L'M'}|VNA>
         between VNA and PAO.
      ****************************************************/
      /****************************************************
       \int VNA(k)*Spe_ProductRF_Bessel(k)*jl(k*R) k^2 dk^3 
      ****************************************************/

      kmin = Radial_kmin;
      kmax = PAO_Nkmax;
      Sk = kmax + kmin;
      Dk = kmax - kmin;

      for (L0=0; L0<=Spe_MaxL_Basis[Cwan]; L0++){
	for (Mul0=0; Mul0<Spe_Num_Basis[Cwan][L0]; Mul0++){
	  for (L1=0; L1<=Spe_MaxL_Basis[Cwan]; L1++){
	    for (Mul1=0; Mul1<Spe_Num_Basis[Cwan][L1]; Mul1++){
	      for (M0=-L0; M0<=L0; M0++){
		for (M1=-L1; M1<=L1; M1++){
		  TmpHNA[L0][Mul0][L0+M0][L1][Mul1][L1+M1] = Complex(0.0,0.0); 
		  TmpHNAr[L0][Mul0][L0+M0][L1][Mul1][L1+M1] = Complex(0.0,0.0); 
		  TmpHNAt[L0][Mul0][L0+M0][L1][Mul1][L1+M1] = Complex(0.0,0.0); 
		  TmpHNAp[L0][Mul0][L0+M0][L1][Mul1][L1+M1] = Complex(0.0,0.0); 
		}
	      }
	    }
	  }
	}
      }


      /* calculate SphB and SphBp */
#ifdef _NEC2
      
      Lmax_Four_Int = 2*Spe_MaxL_Basis[Cwan];
#ifdef _NEC2_ARRAY
      Spherical_Bessel2array(GL_Mesh, GL_NormK, r, Lmax_Four_Int, &(tmp_SphB_D[0]),&(tmp_SphBp_D[0]));
#else
      #pragma _NEC vector
      for (i=0; i<GL_Mesh; i++){
	      const double Normk = GL_NormK[i];
        Spherical_Bessel2(Normk*r,Lmax_Four_Int,&(tmp_SphB_D[i*(Lmax_Four_Int+3)]),&(tmp_SphBp_D[i*(Lmax_Four_Int+3)]));
      }
#endif      
#else
      
      /* allocate SphB and SphBp */

      Lmax_Four_Int = 2*Spe_MaxL_Basis[Cwan];

      SphB = (double**)malloc(sizeof(double*)*(Lmax_Four_Int+3));
      for(LL=0; LL<(Lmax_Four_Int+3); LL++){ 
	SphB[LL] = (double*)malloc(sizeof(double)*GL_Mesh);
      }

      SphBp = (double**)malloc(sizeof(double*)*(Lmax_Four_Int+3));
      for(LL=0; LL<(Lmax_Four_Int+3); LL++){ 
	SphBp[LL] = (double*)malloc(sizeof(double)*GL_Mesh);
      }

      tmp_SphB  = (double*)malloc(sizeof(double)*(Lmax_Four_Int+3));
      tmp_SphBp = (double*)malloc(sizeof(double)*(Lmax_Four_Int+3));

      /* calculate SphB and SphBp */
#ifdef kcomp
#else 
#pragma forceinline recursive
#endif      
      for (i=0; i<GL_Mesh; i++){
	Normk = GL_NormK[i];
	Spherical_Bessel2(Normk*r,Lmax_Four_Int,tmp_SphB,tmp_SphBp);
	for(LL=0; LL<=Lmax_Four_Int; LL++){ 
	  SphB[LL][i]  = tmp_SphB[LL]; 
	  SphBp[LL][i] = tmp_SphBp[LL]; 
	}
      }

      free(tmp_SphB);
      free(tmp_SphBp);
#endif  

      /* loops for L0, Mul0, L1, and Mul1 */

      for (L0=0; L0<=Spe_MaxL_Basis[Cwan]; L0++){
	for (Mul0=0; Mul0<Spe_Num_Basis[Cwan][L0]; Mul0++){
	  for (L1=0; L1<=Spe_MaxL_Basis[Cwan]; L1++){
	    for (Mul1=0; Mul1<Spe_Num_Basis[Cwan][L1]; Mul1++){

	      if (L0<=L1){
#ifdef _NEC2
    /* sum over LL */

    for(LL=0; LL<=2*L1; LL++){ 
#else
		Lmax_Four_Int = 2*L1;//AITUNE: this causes memory leak at free(SphB)

		/* sum over LL */

		for(LL=0; LL<=Lmax_Four_Int; LL++){ 
#endif
		  if (abs(L1-LL)<=L0 && L0<=(L1+LL) ){

		    sum  = 0.0;
		    sumr = 0.0;

		    /* Gauss-Legendre quadrature */

		    for (i=0; i<GL_Mesh; i++){

		      Normk = GL_NormK[i];
#ifdef _NEC2
#ifdef _NEC2_ARRAY
          sj  =  tmp_SphB_D[i + LL*GL_Mesh];
          sjp =  tmp_SphBp_D[i + LL*GL_Mesh];
#else
          sj  =  tmp_SphB_D[i*(Lmax_Four_Int+3) + LL];
          sjp =  tmp_SphBp_D[i*(Lmax_Four_Int+3) + LL];
#endif
#else
		      sj  =  SphB[LL][i];
		      sjp = SphBp[LL][i];
#endif
		      tmp0 = Spe_CrudeVNA_Bessel[Hwan][i];
		      tmp1 = Spe_ProductRF_Bessel[Cwan][L0][Mul0][L1][Mul1][LL][i]; 
		      tmp10 = 0.50*Dk*tmp0*tmp1*GL_Weight[i]*Normk*Normk;

		      sum  += tmp10*sj;
		      sumr += tmp10*sjp*Normk;
		    }

		    /* sum over m */

		    for (M0=-L0; M0<=L0; M0++){
		      for (M1=-L1; M1<=L1; M1++){

			Csum = Complex(0.0,0.0);
			Csumr = Complex(0.0,0.0);
			Csumt = Complex(0.0,0.0);
			Csump = Complex(0.0,0.0);

			for(m=-LL; m<=LL; m++){ 

			  if ( (M1-m)==M0){

			    ComplexSH(LL,m,theta,phi,SH,dSHt,dSHp);
			    CY  = Complex(SH[0],SH[1]);
			    CYt = Complex(dSHt[0],dSHt[1]);
			    CYp = Complex(dSHp[0],dSHp[1]);
#ifdef _NEC_MG
			    gant = pow(-1.0,(double)abs(m))*GetGaunt(L0,M0,L1,M1,LL,-m, Gaunt_list, ai_L1_LIMIT_SQ, ai_LL_LIMIT_SQ);
#else
			    gant = pow(-1.0,(double)abs(m))*Gaunt(L0,M0,L1,M1,LL,-m);
#endif 
			    /* S */

			    Ctmp2 = CRmul(CY,gant);          
			    Csum = Cadd(Csum,Ctmp2);             

			    /* dS/dt */
                      
			    Ctmp2 = CRmul(CYt,gant);          
			    Csumt = Cadd(Csumt,Ctmp2);             
                      
			    /* dS/dp */
                      
			    Ctmp2 = CRmul(CYp,gant);          
			    Csump = Cadd(Csump,Ctmp2);             
			  }

			}

			/* S */
			Ctmp1 = CRmul(Csum,sum);
			TmpHNA[L0][Mul0][L0+M0][L1][Mul1][L1+M1] 
			  = Cadd(TmpHNA[L0][Mul0][L0+M0][L1][Mul1][L1+M1],Ctmp1);

			/* dS/dr */
			Ctmp1 = CRmul(Csum,sumr);
			TmpHNAr[L0][Mul0][L0+M0][L1][Mul1][L1+M1] 
			  = Cadd(TmpHNAr[L0][Mul0][L0+M0][L1][Mul1][L1+M1],Ctmp1);

			/* dS/dt */
			Ctmp1 = CRmul(Csumt,sum);
			TmpHNAt[L0][Mul0][L0+M0][L1][Mul1][L1+M1] 
			  = Cadd(TmpHNAt[L0][Mul0][L0+M0][L1][Mul1][L1+M1],Ctmp1);

			/* dS/dp */
			Ctmp1 = CRmul(Csump,sum);
			TmpHNAp[L0][Mul0][L0+M0][L1][Mul1][L1+M1] 
			  = Cadd(TmpHNAp[L0][Mul0][L0+M0][L1][Mul1][L1+M1],Ctmp1);
		      }
		    }
		  }
		} /* LL */
	      } /* if (L0<=L1) */
	    }
	  }
	}
      }

      /* free SphB and SphBp */
#ifndef _NEC2      
      for(LL=0; LL<(Lmax_Four_Int+3); LL++){ 
	free(SphB[LL]);
      }
      free(SphB);

      for(LL=0; LL<(Lmax_Four_Int+3); LL++){ 
	free(SphBp[LL]);
      }
      free(SphBp);
#endif
      /* copy the upper part to the lower part */

      for (L0=0; L0<=Spe_MaxL_Basis[Cwan]; L0++){
	for (Mul0=0; Mul0<Spe_Num_Basis[Cwan][L0]; Mul0++){
	  for (L1=0; L1<=Spe_MaxL_Basis[Cwan]; L1++){
	    for (Mul1=0; Mul1<Spe_Num_Basis[Cwan][L1]; Mul1++){
	      if (L0<=L1){
		for (M0=-L0; M0<=L0; M0++){
		  for (M1=-L1; M1<=L1; M1++){

		    TmpHNA[L1][Mul1][L1+M1][L0][Mul0][L0+M0] 
		      = Conjg(TmpHNA[L0][Mul0][L0+M0][L1][Mul1][L1+M1]); 

		    TmpHNAr[L1][Mul1][L1+M1][L0][Mul0][L0+M0]
		      = Conjg(TmpHNAr[L0][Mul0][L0+M0][L1][Mul1][L1+M1]);
 
		    TmpHNAt[L1][Mul1][L1+M1][L0][Mul0][L0+M0]
		      = Conjg(TmpHNAt[L0][Mul0][L0+M0][L1][Mul1][L1+M1]); 

		    TmpHNAp[L1][Mul1][L1+M1][L0][Mul0][L0+M0]
		      = Conjg(TmpHNAp[L0][Mul0][L0+M0][L1][Mul1][L1+M1]); 

		  }
		}
	      }
	    }
	  }
	}
      }

      /****************************************************
 	                 Complex to Real
      ****************************************************/

      num0 = 0;
      for (L0=0; L0<=Spe_MaxL_Basis[Cwan]; L0++){
	for (Mul0=0; Mul0<Spe_Num_Basis[Cwan][L0]; Mul0++){
       
	  num1 = 0;
	  for (L1=0; L1<=Spe_MaxL_Basis[Cwan]; L1++){
	    for (Mul1=0; Mul1<Spe_Num_Basis[Cwan][L1]; Mul1++){

	      for (M0=-L0; M0<=L0; M0++){
		for (M1=-L1; M1<=L1; M1++){

		  CsumS0 = Complex(0.0,0.0);
		  CsumSr = Complex(0.0,0.0);
		  CsumSt = Complex(0.0,0.0);
		  CsumSp = Complex(0.0,0.0);

		  for (k=-L0; k<=L0; k++){

		    Ctmp1 = Conjg(Comp2Real[L0][L0+M0][L0+k]);

		    /* S */

		    Ctmp0 = TmpHNA[L0][Mul0][L0+k][L1][Mul1][L1+M1];
		    Ctmp2 = Cmul(Ctmp1,Ctmp0);
		    CsumS0 = Cadd(CsumS0,Ctmp2);

		    /* dS/dr */

		    Ctmp0 = TmpHNAr[L0][Mul0][L0+k][L1][Mul1][L1+M1];
		    Ctmp2 = Cmul(Ctmp1,Ctmp0);
		    CsumSr = Cadd(CsumSr,Ctmp2);

		    /* dS/dt */

		    Ctmp0 = TmpHNAt[L0][Mul0][L0+k][L1][Mul1][L1+M1];
		    Ctmp2 = Cmul(Ctmp1,Ctmp0);
		    CsumSt = Cadd(CsumSt,Ctmp2);

		    /* dS/dp */

		    Ctmp0 = TmpHNAp[L0][Mul0][L0+k][L1][Mul1][L1+M1];
		    Ctmp2 = Cmul(Ctmp1,Ctmp0);
		    CsumSp = Cadd(CsumSp,Ctmp2);

		  }

		  CmatS0[L0+M0][L1+M1] = CsumS0;
		  CmatSr[L0+M0][L1+M1] = CsumSr;
		  CmatSt[L0+M0][L1+M1] = CsumSt;
		  CmatSp[L0+M0][L1+M1] = CsumSp;

		}
	      }

	      for (M0=-L0; M0<=L0; M0++){
		for (M1=-L1; M1<=L1; M1++){

		  CsumS0 = Complex(0.0,0.0);
		  CsumSr = Complex(0.0,0.0);
		  CsumSt = Complex(0.0,0.0);
		  CsumSp = Complex(0.0,0.0);

		  for (k=-L1; k<=L1; k++){

		    /* S */ 

		    Ctmp1 = Cmul(CmatS0[L0+M0][L1+k],Comp2Real[L1][L1+M1][L1+k]);
		    CsumS0 = Cadd(CsumS0,Ctmp1);

		    /* dS/dr */ 

		    Ctmp1 = Cmul(CmatSr[L0+M0][L1+k],Comp2Real[L1][L1+M1][L1+k]);
		    CsumSr = Cadd(CsumSr,Ctmp1);

		    /* dS/dt */ 

		    Ctmp1 = Cmul(CmatSt[L0+M0][L1+k],Comp2Real[L1][L1+M1][L1+k]);
		    CsumSt = Cadd(CsumSt,Ctmp1);

		    /* dS/dp */ 

		    Ctmp1 = Cmul(CmatSp[L0+M0][L1+k],Comp2Real[L1][L1+M1][L1+k]);
		    CsumSp = Cadd(CsumSp,Ctmp1);

		  }

		  HVNA2[0][Mc_AN][h_AN][num0+L0+M0][num1+L1+M1] = 8.0*CsumS0.r;

		  if (h_AN!=0){

		    if (fabs(siT)<10e-14){

		      HVNA2[1][Mc_AN][h_AN][num0+L0+M0][num1+L1+M1] =
			-8.0*(siT*coP*CsumSr.r + coT*coP/r*CsumSt.r);

		      HVNA2[2][Mc_AN][h_AN][num0+L0+M0][num1+L1+M1] =
			-8.0*(siT*siP*CsumSr.r + coT*siP/r*CsumSt.r);

		      HVNA2[3][Mc_AN][h_AN][num0+L0+M0][num1+L1+M1] =
			-8.0*(coT*CsumSr.r - siT/r*CsumSt.r);
		    }

		    else{

		      HVNA2[1][Mc_AN][h_AN][num0+L0+M0][num1+L1+M1] =
			-8.0*(siT*coP*CsumSr.r + coT*coP/r*CsumSt.r
			     - siP/siT/r*CsumSp.r);

		      HVNA2[2][Mc_AN][h_AN][num0+L0+M0][num1+L1+M1] =
			-8.0*(siT*siP*CsumSr.r + coT*siP/r*CsumSt.r
			     + coP/siT/r*CsumSp.r);

		      HVNA2[3][Mc_AN][h_AN][num0+L0+M0][num1+L1+M1] =
			-8.0*(coT*CsumSr.r - siT/r*CsumSt.r);
		    }
		  }
		  else{
		    HVNA2[1][Mc_AN][h_AN][num0+L0+M0][num1+L1+M1] = 0.0;
		    HVNA2[2][Mc_AN][h_AN][num0+L0+M0][num1+L1+M1] = 0.0;
		    HVNA2[3][Mc_AN][h_AN][num0+L0+M0][num1+L1+M1] = 0.0;
		  }

		}
	      }

	      num1 = num1 + 2*L1 + 1; 
	    }
	  }

	  num0 = num0 + 2*L0 + 1; 
	}
      }

      dtime(&Etime_atom);
      time_per_atom[Gc_AN] += Etime_atom - Stime_atom;

    } /* Mc_AN */

    /* freeing of arrays */
    
#ifdef _NEC2
  	free(tmp_SphB_D);
	  free(tmp_SphBp_D);
#endif

    for (i=0; i<(List_YOUSO[25]+1); i++){
      for (j=0; j<List_YOUSO[24]; j++){
	for (k=0; k<(2*(List_YOUSO[25]+1)+1); k++){
	  for (l=0; l<(List_YOUSO[25]+1); l++){
	    for (m=0; m<List_YOUSO[24]; m++){
	      free(TmpHNA[i][j][k][l][m]);
	    }
	    free(TmpHNA[i][j][k][l]);
	  }
	  free(TmpHNA[i][j][k]);
	}
	free(TmpHNA[i][j]);
      }
      free(TmpHNA[i]);
    }
    free(TmpHNA);

    for (i=0; i<(List_YOUSO[25]+1); i++){
      for (j=0; j<List_YOUSO[24]; j++){
	for (k=0; k<(2*(List_YOUSO[25]+1)+1); k++){
	  for (l=0; l<(List_YOUSO[25]+1); l++){
	    for (m=0; m<List_YOUSO[24]; m++){
	      free(TmpHNAr[i][j][k][l][m]);
	    }
	    free(TmpHNAr[i][j][k][l]);
	  }
	  free(TmpHNAr[i][j][k]);
	}
	free(TmpHNAr[i][j]);
      }
      free(TmpHNAr[i]);
    }
    free(TmpHNAr);

    for (i=0; i<(List_YOUSO[25]+1); i++){
      for (j=0; j<List_YOUSO[24]; j++){
	for (k=0; k<(2*(List_YOUSO[25]+1)+1); k++){
	  for (l=0; l<(List_YOUSO[25]+1); l++){
	    for (m=0; m<List_YOUSO[24]; m++){
	      free(TmpHNAt[i][j][k][l][m]);
	    }
	    free(TmpHNAt[i][j][k][l]);
	  }
	  free(TmpHNAt[i][j][k]);
	}
	free(TmpHNAt[i][j]);
      }
      free(TmpHNAt[i]);
    }
    free(TmpHNAt);

    for (i=0; i<(List_YOUSO[25]+1); i++){
      for (j=0; j<List_YOUSO[24]; j++){
	for (k=0; k<(2*(List_YOUSO[25]+1)+1); k++){
	  for (l=0; l<(List_YOUSO[25]+1); l++){
	    for (m=0; m<List_YOUSO[24]; m++){
	      free(TmpHNAp[i][j][k][l][m]);
	    }
	    free(TmpHNAp[i][j][k][l]);
	  }
	  free(TmpHNAp[i][j][k]);
	}
	free(TmpHNAp[i][j]);
      }
      free(TmpHNAp[i]);
    }
    free(TmpHNAp);

    for (i=0; i<(2*(List_YOUSO[25]+1)+1); i++){
      free(CmatS0[i]);
    }
    free(CmatS0);

    for (i=0; i<(2*(List_YOUSO[25]+1)+1); i++){
      free(CmatSr[i]);
    }
    free(CmatSr);

    for (i=0; i<(2*(List_YOUSO[25]+1)+1); i++){
      free(CmatSt[i]);
    }
    free(CmatSt);

    for (i=0; i<(2*(List_YOUSO[25]+1)+1); i++){
      free(CmatSp[i]);
    }
    free(CmatSp);

#pragma omp flush(HVNA2)

  } /* #pragma omp parallel */

  /****************************************************
    HVNA[Mc_AN][0] = sum_{h_AN} HVNA2
  ****************************************************/

  for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){

    Gc_AN = M2G[Mc_AN];
    Cwan = WhatSpecies[Gc_AN];
    for (i=0; i<Spe_Total_NO[Cwan]; i++){
      for (j=0; j<Spe_Total_NO[Cwan]; j++){
	HVNA[Mc_AN][0][i][j] = 0.0;
      }
    }

    for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
      for (i=0; i<Spe_Total_NO[Cwan]; i++){
        for (j=0; j<Spe_Total_NO[Cwan]; j++){
          HVNA[Mc_AN][0][i][j] += HVNA2[0][Mc_AN][h_AN][i][j];
	}
      }
    }
  } 

  /****************************************************
  freeing of arrays:
  ****************************************************/

  free(OneD2Mc_AN);
  free(OneD2h_AN);
  free(Gaunt_list);

  /* for time */
  dtime(&TEtime);
  time0 = TEtime - TStime;
  return time0;
} 



double Set_VNA3(double *****HVNA3)
{
  /****************************************************
   Evaluate matrix elements of one-center (orbitals)
   but two-center integrals for neutral atom potentials
   in the momentum space.

   <VNA|Phi_{LM,L'M'}>
  ****************************************************/
  static int firsttime=0; /* due to the same as in Set_VNA2 */
  int L2,L3,L,GL,i,j;
  int k,kl,h_AN,Gh_AN,Rnh,Ls,n;
  int tno0,tno1,tno2,i1,j1,p;
  int Cwan,Hwan,fan,jg,kg,wakg,Lmax;
  int size_TmpHNA;
  int Mc_AN,Gc_AN,Mj_AN,num,size1,size2;
  double time0;
  double tmp2,tmp3,tmp4,tmp5;
  double TStime,TEtime;
  int Num_RVNA;

  /* for OpenMP */
  int OneD_Nloop,*OneD2Mc_AN,*OneD2h_AN;

  dtime(&TStime);

  /****************************************************
  allocation of arrays:
  ****************************************************/

  size_TmpHNA = (List_YOUSO[25]+1)*List_YOUSO[24]*
                (2*(List_YOUSO[25]+1)+1)*
                (List_YOUSO[25]+1)*List_YOUSO[24]*
                (2*(List_YOUSO[25]+1)+1);

  /* PrintMemory */
  if (firsttime) {
    PrintMemory("Set_ProExpn_VNA: TmpHNA",sizeof(double)*size_TmpHNA,NULL);
    PrintMemory("Set_ProExpn_VNA: TmpHNAr",sizeof(double)*size_TmpHNA,NULL);
    PrintMemory("Set_ProExpn_VNA: TmpHNAt",sizeof(double)*size_TmpHNA,NULL);
    PrintMemory("Set_ProExpn_VNA: TmpHNAp",sizeof(double)*size_TmpHNA,NULL);

    firsttime=0;
  }

  /* one-dimensionalize the Mc_AN and h_AN loops */

  OneD_Nloop = 0;
  for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
    Gc_AN = M2G[Mc_AN];    
    for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
      OneD_Nloop++;
    }
  }  

  OneD2Mc_AN = (int*)malloc(sizeof(int)*(OneD_Nloop+1));
  OneD2h_AN = (int*)malloc(sizeof(int)*(OneD_Nloop+1));

  OneD_Nloop = 0;
  for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
    Gc_AN = M2G[Mc_AN];    
    for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
      OneD2Mc_AN[OneD_Nloop] = Mc_AN; 
      OneD2h_AN[OneD_Nloop] = h_AN; 
      OneD_Nloop++;
    }
  }

ftrace_region_begin("MG3");
  int max_L0 = 0;
  for (int spe=0; spe<SpeciesNum; spe++){
    if(max_L0 < Spe_MaxL_Basis[spe]) max_L0 = Spe_MaxL_Basis[spe];
  }
  const int ai_L0_LIMIT = max_L0 + 1;
  const int ai_L0_LIMIT_SQ = ai_L0_LIMIT * ai_L0_LIMIT;
  const int ai_L1_LIMIT = ai_L0_LIMIT;
  const int ai_L1_LIMIT_SQ = ai_L0_LIMIT_SQ;
  const int ai_LL_LIMIT = max_L0*2 + 1;
  const int ai_LL_LIMIT_SQ = ai_LL_LIMIT*ai_LL_LIMIT;
  double* Gaunt_list = MakeGaunt(ai_L0_LIMIT_SQ, ai_L1_LIMIT_SQ, ai_LL_LIMIT_SQ);
ftrace_region_end("MG3");

#pragma omp parallel shared(List_YOUSO,time_per_atom,HVNA3,Comp2Real,GL_Weight,Spe_ProductRF_Bessel,Spe_CrudeVNA_Bessel,GL_NormK,Spe_Num_Basis,Spe_MaxL_Basis,PAO_Nkmax,atv,Gxyz,WhatSpecies,ncn,natn,M2G,OneD2h_AN,OneD2Mc_AN,OneD_Nloop) 
  {         
    int OMPID,Nthrds,Nprocs,Nloop;
    int Mc_AN,h_AN,Gc_AN,Cwan,Gh_AN;
    int Rnh,Hwan,L0,Mul0,L1,Mul1,M0,M1,LL;
    int Lmax_Four_Int,i,j,k,l,m,num0,num1;
    double dx,dy,dz;
    double siT,coT,siP,coP;
    double Normk,kmin,kmax,Sk,Dk;
    double gant,r,theta,phi;
    double SH[2],dSHt[2],dSHp[2];
    double S_coordinate[3];
    double Stime_atom,Etime_atom;
    double sum,sumr,sj,sjp,tmp0,tmp1,tmp10;
#ifndef _NEC1
    double **SphB,**SphBp;
    double *tmp_SphB,*tmp_SphBp;
#endif
    dcomplex Ctmp0,Ctmp1,Ctmp2,Cpow;
    dcomplex Csum,Csumr,Csumt,Csump;
    dcomplex CsumS0,CsumSr,CsumSt,CsumSp;
    dcomplex CY,CYt,CYp,CY1,CYt1,CYp1;
    dcomplex **CmatS0;
    dcomplex **CmatSr;
    dcomplex **CmatSt;
    dcomplex **CmatSp;
    dcomplex ******TmpHNA;
    dcomplex ******TmpHNAr;
    dcomplex ******TmpHNAt;
    dcomplex ******TmpHNAp;

    /* allocation of arrays */

    TmpHNA = (dcomplex******)malloc(sizeof(dcomplex*****)*(List_YOUSO[25]+1));
    for (i=0; i<(List_YOUSO[25]+1); i++){
      TmpHNA[i] = (dcomplex*****)malloc(sizeof(dcomplex****)*List_YOUSO[24]);
      for (j=0; j<List_YOUSO[24]; j++){
	TmpHNA[i][j] = (dcomplex****)malloc(sizeof(dcomplex***)*(2*(List_YOUSO[25]+1)+1));
	for (k=0; k<(2*(List_YOUSO[25]+1)+1); k++){
	  TmpHNA[i][j][k] = (dcomplex***)malloc(sizeof(dcomplex**)*(List_YOUSO[25]+1));
	  for (l=0; l<(List_YOUSO[25]+1); l++){
	    TmpHNA[i][j][k][l] = (dcomplex**)malloc(sizeof(dcomplex*)*List_YOUSO[24]);
	    for (m=0; m<List_YOUSO[24]; m++){
	      TmpHNA[i][j][k][l][m] = (dcomplex*)malloc(sizeof(dcomplex)*(2*(List_YOUSO[25]+1)+1));
	    }
	  }
	}
      }
    }

    TmpHNAr = (dcomplex******)malloc(sizeof(dcomplex*****)*(List_YOUSO[25]+1));
    for (i=0; i<(List_YOUSO[25]+1); i++){
      TmpHNAr[i] = (dcomplex*****)malloc(sizeof(dcomplex****)*List_YOUSO[24]);
      for (j=0; j<List_YOUSO[24]; j++){
	TmpHNAr[i][j] = (dcomplex****)malloc(sizeof(dcomplex***)*(2*(List_YOUSO[25]+1)+1));
	for (k=0; k<(2*(List_YOUSO[25]+1)+1); k++){
	  TmpHNAr[i][j][k] = (dcomplex***)malloc(sizeof(dcomplex**)*(List_YOUSO[25]+1));
	  for (l=0; l<(List_YOUSO[25]+1); l++){
	    TmpHNAr[i][j][k][l] = (dcomplex**)malloc(sizeof(dcomplex*)*List_YOUSO[24]);
	    for (m=0; m<List_YOUSO[24]; m++){
	      TmpHNAr[i][j][k][l][m] = (dcomplex*)malloc(sizeof(dcomplex)*(2*(List_YOUSO[25]+1)+1));
	    }
	  }
	}
      }
    }

    TmpHNAt = (dcomplex******)malloc(sizeof(dcomplex*****)*(List_YOUSO[25]+1));
    for (i=0; i<(List_YOUSO[25]+1); i++){
      TmpHNAt[i] = (dcomplex*****)malloc(sizeof(dcomplex****)*List_YOUSO[24]);
      for (j=0; j<List_YOUSO[24]; j++){
	TmpHNAt[i][j] = (dcomplex****)malloc(sizeof(dcomplex***)*(2*(List_YOUSO[25]+1)+1));
	for (k=0; k<(2*(List_YOUSO[25]+1)+1); k++){
	  TmpHNAt[i][j][k] = (dcomplex***)malloc(sizeof(dcomplex**)*(List_YOUSO[25]+1));
	  for (l=0; l<(List_YOUSO[25]+1); l++){
	    TmpHNAt[i][j][k][l] = (dcomplex**)malloc(sizeof(dcomplex*)*List_YOUSO[24]);
	    for (m=0; m<List_YOUSO[24]; m++){
	      TmpHNAt[i][j][k][l][m] = (dcomplex*)malloc(sizeof(dcomplex)*(2*(List_YOUSO[25]+1)+1));
	    }
	  }
	}
      }
    }

    TmpHNAp = (dcomplex******)malloc(sizeof(dcomplex*****)*(List_YOUSO[25]+1));
    for (i=0; i<(List_YOUSO[25]+1); i++){
      TmpHNAp[i] = (dcomplex*****)malloc(sizeof(dcomplex****)*List_YOUSO[24]);
      for (j=0; j<List_YOUSO[24]; j++){
	TmpHNAp[i][j] = (dcomplex****)malloc(sizeof(dcomplex***)*(2*(List_YOUSO[25]+1)+1));
	for (k=0; k<(2*(List_YOUSO[25]+1)+1); k++){
	  TmpHNAp[i][j][k] = (dcomplex***)malloc(sizeof(dcomplex**)*(List_YOUSO[25]+1));
	  for (l=0; l<(List_YOUSO[25]+1); l++){
	    TmpHNAp[i][j][k][l] = (dcomplex**)malloc(sizeof(dcomplex*)*List_YOUSO[24]);
	    for (m=0; m<List_YOUSO[24]; m++){
	      TmpHNAp[i][j][k][l][m] = (dcomplex*)malloc(sizeof(dcomplex)*(2*(List_YOUSO[25]+1)+1));
	    }
	  }
	}
      }
    }

    CmatS0 = (dcomplex**)malloc(sizeof(dcomplex*)*(2*(List_YOUSO[25]+1)+1));
    for (i=0; i<(2*(List_YOUSO[25]+1)+1); i++){
      CmatS0[i] = (dcomplex*)malloc(sizeof(dcomplex)*(2*(List_YOUSO[25]+1)+1));
    }

    CmatSr = (dcomplex**)malloc(sizeof(dcomplex*)*(2*(List_YOUSO[25]+1)+1));
    for (i=0; i<(2*(List_YOUSO[25]+1)+1); i++){
      CmatSr[i] = (dcomplex*)malloc(sizeof(dcomplex)*(2*(List_YOUSO[25]+1)+1));
    }

    CmatSt = (dcomplex**)malloc(sizeof(dcomplex*)*(2*(List_YOUSO[25]+1)+1));
    for (i=0; i<(2*(List_YOUSO[25]+1)+1); i++){
      CmatSt[i] = (dcomplex*)malloc(sizeof(dcomplex)*(2*(List_YOUSO[25]+1)+1));
    }

    CmatSp = (dcomplex**)malloc(sizeof(dcomplex*)*(2*(List_YOUSO[25]+1)+1));
    for (i=0; i<(2*(List_YOUSO[25]+1)+1); i++){
      CmatSp[i] = (dcomplex*)malloc(sizeof(dcomplex)*(2*(List_YOUSO[25]+1)+1));
    }

#ifdef _NEC1
    double* tmp_SphB_D  = (double*)malloc(sizeof(double)*30*GL_Mesh);
    double* tmp_SphBp_D = (double*)malloc(sizeof(double)*30*GL_Mesh);
#endif

    /* get info. on OpenMP */ 

    OMPID = omp_get_thread_num();
    Nthrds = omp_get_num_threads();
    Nprocs = omp_get_num_procs();

    /* one-dimensionalized loop */

    for (Nloop=OMPID*OneD_Nloop/Nthrds; Nloop<(OMPID+1)*OneD_Nloop/Nthrds; Nloop++){

      dtime(&Stime_atom);

      /* get Mc_AN and h_AN */

      Mc_AN = OneD2Mc_AN[Nloop];
      h_AN  = OneD2h_AN[Nloop];

      /* set data on Mc_AN */

      Gc_AN = M2G[Mc_AN];
      Cwan = WhatSpecies[Gc_AN];

      /* set data on h_AN */

      Gh_AN = natn[Gc_AN][h_AN];        
      Rnh = ncn[Gc_AN][h_AN];
      Hwan = WhatSpecies[Gh_AN];

      dx = Gxyz[Gc_AN][1] - Gxyz[Gh_AN][1] - atv[Rnh][1];
      dy = Gxyz[Gc_AN][2] - Gxyz[Gh_AN][2] - atv[Rnh][2];
      dz = Gxyz[Gc_AN][3] - Gxyz[Gh_AN][3] - atv[Rnh][3];

      xyz2spherical(dx,dy,dz,0.0,0.0,0.0,S_coordinate); 
      r     = S_coordinate[0];
      theta = S_coordinate[1];
      phi   = S_coordinate[2];

      /* for empty atoms or finite elemens basis */

      if (r<1.0e-10) r = 1.0e-10;

      /* precalculation of sin and cos */

      siT = sin(theta);
      coT = cos(theta);
      siP = sin(phi);
      coP = cos(phi);

      /****************************************************
         evaluate ovelap integrals <Phi_{LM,L'M'}|VNA>
         between VNA and PAO.
      ****************************************************/
      /****************************************************
       \int VNA(k)*Spe_ProductRF_Bessel(k)*jl(k*R) k^2 dk^3 
      ****************************************************/

      kmin = Radial_kmin;
      kmax = PAO_Nkmax;
      Sk = kmax + kmin;
      Dk = kmax - kmin;

      for (L0=0; L0<=Spe_MaxL_Basis[Hwan]; L0++){
	for (Mul0=0; Mul0<Spe_Num_Basis[Hwan][L0]; Mul0++){
	  for (L1=0; L1<=Spe_MaxL_Basis[Hwan]; L1++){
	    for (Mul1=0; Mul1<Spe_Num_Basis[Hwan][L1]; Mul1++){
	      for (M0=-L0; M0<=L0; M0++){
		for (M1=-L1; M1<=L1; M1++){
		  TmpHNA[L0][Mul0][L0+M0][L1][Mul1][L1+M1] = Complex(0.0,0.0); 
		  TmpHNAr[L0][Mul0][L0+M0][L1][Mul1][L1+M1] = Complex(0.0,0.0); 
		  TmpHNAt[L0][Mul0][L0+M0][L1][Mul1][L1+M1] = Complex(0.0,0.0); 
		  TmpHNAp[L0][Mul0][L0+M0][L1][Mul1][L1+M1] = Complex(0.0,0.0); 
		}
	      }
	    }
	  }
	}
      }

#ifdef _NEC1
      /* calculate SphB and SphBp */
      Lmax_Four_Int = 2*Spe_MaxL_Basis[Hwan];
#ifdef _NEC1_ARRAY
      Spherical_Bessel2array(GL_Mesh, GL_NormK, r, Lmax_Four_Int, &(tmp_SphB_D[0]),&(tmp_SphBp_D[0]));
#else
      #pragma _NEC_vector
      for (i=0; i<GL_Mesh; i++){
	      const double Normk = GL_NormK[i];
        Spherical_Bessel2(Normk*r,Lmax_Four_Int,&(tmp_SphB_D[i*(Lmax_Four_Int+3)]),&(tmp_SphBp_D[i*(Lmax_Four_Int+3)]));
      }
#endif
#else
      /* allocate SphB and SphBp */

      Lmax_Four_Int = 2*Spe_MaxL_Basis[Hwan];

      SphB = (double**)malloc(sizeof(double*)*(Lmax_Four_Int+3));
      for(LL=0; LL<(Lmax_Four_Int+3); LL++){ 
	SphB[LL] = (double*)malloc(sizeof(double)*GL_Mesh);
      }

      SphBp = (double**)malloc(sizeof(double*)*(Lmax_Four_Int+3));
      for(LL=0; LL<(Lmax_Four_Int+3); LL++){ 
	SphBp[LL] = (double*)malloc(sizeof(double)*GL_Mesh);
      }

      tmp_SphB  = (double*)malloc(sizeof(double)*(Lmax_Four_Int+3));
      tmp_SphBp = (double*)malloc(sizeof(double)*(Lmax_Four_Int+3));

      /* calculate SphB and SphBp */
#ifdef kcomp
#else 
#pragma forceinline recursive
#endif      
      for (i=0; i<GL_Mesh; i++){
	Normk = GL_NormK[i];
	Spherical_Bessel2(Normk*r,Lmax_Four_Int,tmp_SphB,tmp_SphBp);
	for(LL=0; LL<=Lmax_Four_Int; LL++){ 
	  SphB[LL][i]  = tmp_SphB[LL]; 
	  SphBp[LL][i] = tmp_SphBp[LL]; 
	}
      }

      free(tmp_SphB);
      free(tmp_SphBp);
#endif
      /* loops for L0, Mul0, L1, and Mul1 */

      for (L0=0; L0<=Spe_MaxL_Basis[Hwan]; L0++){
	for (Mul0=0; Mul0<Spe_Num_Basis[Hwan][L0]; Mul0++){
	  for (L1=0; L1<=Spe_MaxL_Basis[Hwan]; L1++){
	    for (Mul1=0; Mul1<Spe_Num_Basis[Hwan][L1]; Mul1++){

	      if (L0<=L1){
#ifdef _NEC1
		/* sum over LL */

		for(LL=0; LL<=2*L1; LL++){ 
#else
		Lmax_Four_Int = 2*L1;  //AITUNE: this causes memory leak

		/* sum over LL */

		for(LL=0; LL<=Lmax_Four_Int; LL++){ 
#endif
		  if (abs(L1-LL)<=L0 && L0<=(L1+LL) ){

		    sum  = 0.0;
		    sumr = 0.0;

		    /* Gauss-Legendre quadrature */

		    for (i=0; i<GL_Mesh; i++){

		      Normk = GL_NormK[i];
#ifdef _NEC1
#ifdef _NEC1_ARRAY
          sj = tmp_SphB_D[i + LL*GL_Mesh];
          sjp = tmp_SphBp_D[i + LL*GL_Mesh];
#else
          sj = tmp_SphB_D[i*(Lmax_Four_Int+3) + LL];
          sjp = tmp_SphBp_D[i*(Lmax_Four_Int+3) + LL];
#endif          
#else          
		      sj  =  SphB[LL][i];
		      sjp = SphBp[LL][i];
#endif
		      tmp0 = Spe_CrudeVNA_Bessel[Cwan][i];
		      tmp1 = Spe_ProductRF_Bessel[Hwan][L0][Mul0][L1][Mul1][LL][i]; 
		      tmp10 = 0.50*Dk*tmp0*tmp1*GL_Weight[i]*Normk*Normk;

		      sum  += tmp10*sj;
		      sumr += tmp10*sjp*Normk;
		    }

		    /* sum over m */

		    for (M0=-L0; M0<=L0; M0++){
		      for (M1=-L1; M1<=L1; M1++){

			Csum = Complex(0.0,0.0);
			Csumr = Complex(0.0,0.0);
			Csumt = Complex(0.0,0.0);
			Csump = Complex(0.0,0.0);

			for(m=-LL; m<=LL; m++){ 

			  if ( (M1-m)==M0){

			    ComplexSH(LL,m,theta,phi,SH,dSHt,dSHp);
			    CY  = Complex(SH[0],SH[1]);
			    CYt = Complex(dSHt[0],dSHt[1]);
			    CYp = Complex(dSHp[0],dSHp[1]);
#ifdef _NEC_MG
			    gant = pow(-1.0,(double)abs(m))*GetGaunt(L0,M0,L1,M1,LL,-m, Gaunt_list, ai_L1_LIMIT_SQ, ai_LL_LIMIT_SQ);
#else
			    gant = pow(-1.0,(double)abs(m))*Gaunt(L0,M0,L1,M1,LL,-m);
#endif
			    /* S */

			    Ctmp2 = CRmul(CY,gant);          
			    Csum = Cadd(Csum,Ctmp2);             

			    /* dS/dt */
                      
			    Ctmp2 = CRmul(CYt,gant);          
			    Csumt = Cadd(Csumt,Ctmp2);             
                      
			    /* dS/dp */
                      
			    Ctmp2 = CRmul(CYp,gant);          
			    Csump = Cadd(Csump,Ctmp2);             
			  }

			}

			/* S */
			Ctmp1 = CRmul(Csum,sum);
			TmpHNA[L0][Mul0][L0+M0][L1][Mul1][L1+M1] 
			  = Cadd(TmpHNA[L0][Mul0][L0+M0][L1][Mul1][L1+M1],Ctmp1);

			/* dS/dr */
			Ctmp1 = CRmul(Csum,sumr);
			TmpHNAr[L0][Mul0][L0+M0][L1][Mul1][L1+M1] 
			  = Cadd(TmpHNAr[L0][Mul0][L0+M0][L1][Mul1][L1+M1],Ctmp1);

			/* dS/dt */
			Ctmp1 = CRmul(Csumt,sum);
			TmpHNAt[L0][Mul0][L0+M0][L1][Mul1][L1+M1] 
			  = Cadd(TmpHNAt[L0][Mul0][L0+M0][L1][Mul1][L1+M1],Ctmp1);

			/* dS/dp */
			Ctmp1 = CRmul(Csump,sum);
			TmpHNAp[L0][Mul0][L0+M0][L1][Mul1][L1+M1] 
			  = Cadd(TmpHNAp[L0][Mul0][L0+M0][L1][Mul1][L1+M1],Ctmp1);
		      }
		    }
		  }
		} /* LL */
	      } /* if (L0<=L1) */
	    }
	  }
	}
      }

      /* free SphB and SphBp */
#ifndef _NEC1    
      for(LL=0; LL<(Lmax_Four_Int+3); LL++){ 
	free(SphB[LL]);
      }
      free(SphB);

      for(LL=0; LL<(Lmax_Four_Int+3); LL++){ 
	free(SphBp[LL]);
      }
      free(SphBp);
#endif
      /* copy the upper part to the lower part */

      for (L0=0; L0<=Spe_MaxL_Basis[Hwan]; L0++){
	for (Mul0=0; Mul0<Spe_Num_Basis[Hwan][L0]; Mul0++){
	  for (L1=0; L1<=Spe_MaxL_Basis[Hwan]; L1++){
	    for (Mul1=0; Mul1<Spe_Num_Basis[Hwan][L1]; Mul1++){
	      if (L0<=L1){
		for (M0=-L0; M0<=L0; M0++){
		  for (M1=-L1; M1<=L1; M1++){

		    TmpHNA[L1][Mul1][L1+M1][L0][Mul0][L0+M0] 
		      = Conjg(TmpHNA[L0][Mul0][L0+M0][L1][Mul1][L1+M1]); 

		    TmpHNAr[L1][Mul1][L1+M1][L0][Mul0][L0+M0]
		      = Conjg(TmpHNAr[L0][Mul0][L0+M0][L1][Mul1][L1+M1]);
 
		    TmpHNAt[L1][Mul1][L1+M1][L0][Mul0][L0+M0]
		      = Conjg(TmpHNAt[L0][Mul0][L0+M0][L1][Mul1][L1+M1]); 

		    TmpHNAp[L1][Mul1][L1+M1][L0][Mul0][L0+M0]
		      = Conjg(TmpHNAp[L0][Mul0][L0+M0][L1][Mul1][L1+M1]); 

		  }
		}
	      }
	    }
	  }
	}
      }

      /****************************************************
 	                 Complex to Real
      ****************************************************/

      num0 = 0;
      for (L0=0; L0<=Spe_MaxL_Basis[Hwan]; L0++){
	for (Mul0=0; Mul0<Spe_Num_Basis[Hwan][L0]; Mul0++){
       
	  num1 = 0;
	  for (L1=0; L1<=Spe_MaxL_Basis[Hwan]; L1++){
	    for (Mul1=0; Mul1<Spe_Num_Basis[Hwan][L1]; Mul1++){

	      for (M0=-L0; M0<=L0; M0++){
		for (M1=-L1; M1<=L1; M1++){

		  CsumS0 = Complex(0.0,0.0);
		  CsumSr = Complex(0.0,0.0);
		  CsumSt = Complex(0.0,0.0);
		  CsumSp = Complex(0.0,0.0);

		  for (k=-L0; k<=L0; k++){

		    Ctmp1 = Conjg(Comp2Real[L0][L0+M0][L0+k]);

		    /* S */

		    Ctmp0 = TmpHNA[L0][Mul0][L0+k][L1][Mul1][L1+M1];
		    Ctmp2 = Cmul(Ctmp1,Ctmp0);
		    CsumS0 = Cadd(CsumS0,Ctmp2);

		    /* dS/dr */

		    Ctmp0 = TmpHNAr[L0][Mul0][L0+k][L1][Mul1][L1+M1];
		    Ctmp2 = Cmul(Ctmp1,Ctmp0);
		    CsumSr = Cadd(CsumSr,Ctmp2);

		    /* dS/dt */

		    Ctmp0 = TmpHNAt[L0][Mul0][L0+k][L1][Mul1][L1+M1];
		    Ctmp2 = Cmul(Ctmp1,Ctmp0);
		    CsumSt = Cadd(CsumSt,Ctmp2);

		    /* dS/dp */

		    Ctmp0 = TmpHNAp[L0][Mul0][L0+k][L1][Mul1][L1+M1];
		    Ctmp2 = Cmul(Ctmp1,Ctmp0);
		    CsumSp = Cadd(CsumSp,Ctmp2);

		  }

		  CmatS0[L0+M0][L1+M1] = CsumS0;
		  CmatSr[L0+M0][L1+M1] = CsumSr;
		  CmatSt[L0+M0][L1+M1] = CsumSt;
		  CmatSp[L0+M0][L1+M1] = CsumSp;

		}
	      }

	      for (M0=-L0; M0<=L0; M0++){
		for (M1=-L1; M1<=L1; M1++){

		  CsumS0 = Complex(0.0,0.0);
		  CsumSr = Complex(0.0,0.0);
		  CsumSt = Complex(0.0,0.0);
		  CsumSp = Complex(0.0,0.0);

		  for (k=-L1; k<=L1; k++){

		    /* S */ 


		    Ctmp1 = Cmul(CmatS0[L0+M0][L1+k],Comp2Real[L1][L1+M1][L1+k]);
		    CsumS0 = Cadd(CsumS0,Ctmp1);

		    /* dS/dr */ 

		    Ctmp1 = Cmul(CmatSr[L0+M0][L1+k],Comp2Real[L1][L1+M1][L1+k]);
		    CsumSr = Cadd(CsumSr,Ctmp1);

		    /* dS/dt */ 

		    Ctmp1 = Cmul(CmatSt[L0+M0][L1+k],Comp2Real[L1][L1+M1][L1+k]);
		    CsumSt = Cadd(CsumSt,Ctmp1);

		    /* dS/dp */ 

		    Ctmp1 = Cmul(CmatSp[L0+M0][L1+k],Comp2Real[L1][L1+M1][L1+k]);
		    CsumSp = Cadd(CsumSp,Ctmp1);

		  }

		  HVNA3[0][Mc_AN][h_AN][num0+L0+M0][num1+L1+M1] = 8.0*CsumS0.r;

		  if (h_AN!=0){

		    if (fabs(siT)<10e-14){

		      HVNA3[1][Mc_AN][h_AN][num0+L0+M0][num1+L1+M1] =
			-8.0*(siT*coP*CsumSr.r + coT*coP/r*CsumSt.r);

		      HVNA3[2][Mc_AN][h_AN][num0+L0+M0][num1+L1+M1] =
			-8.0*(siT*siP*CsumSr.r + coT*siP/r*CsumSt.r);

		      HVNA3[3][Mc_AN][h_AN][num0+L0+M0][num1+L1+M1] =
			-8.0*(coT*CsumSr.r - siT/r*CsumSt.r);
		    }

		    else{

		      HVNA3[1][Mc_AN][h_AN][num0+L0+M0][num1+L1+M1] =
			-8.0*(siT*coP*CsumSr.r + coT*coP/r*CsumSt.r
			     - siP/siT/r*CsumSp.r);

		      HVNA3[2][Mc_AN][h_AN][num0+L0+M0][num1+L1+M1] =
			-8.0*(siT*siP*CsumSr.r + coT*siP/r*CsumSt.r
			     + coP/siT/r*CsumSp.r);

		      HVNA3[3][Mc_AN][h_AN][num0+L0+M0][num1+L1+M1] =
			-8.0*(coT*CsumSr.r - siT/r*CsumSt.r);
		    }
		  }
		  else{
		    HVNA3[1][Mc_AN][h_AN][num0+L0+M0][num1+L1+M1] = 0.0;
		    HVNA3[2][Mc_AN][h_AN][num0+L0+M0][num1+L1+M1] = 0.0;
		    HVNA3[3][Mc_AN][h_AN][num0+L0+M0][num1+L1+M1] = 0.0;
		  }

		}
	      }

	      num1 = num1 + 2*L1 + 1; 
	    }
	  }

	  num0 = num0 + 2*L0 + 1; 
	}
      }

      dtime(&Etime_atom);
      time_per_atom[Gc_AN] += Etime_atom - Stime_atom;

    } /* Mc_AN */

    /* freeing of arrays */
    
#ifdef _NEC1
  	free(tmp_SphB_D);
	  free(tmp_SphBp_D);
#endif

    for (i=0; i<(List_YOUSO[25]+1); i++){
      for (j=0; j<List_YOUSO[24]; j++){
	for (k=0; k<(2*(List_YOUSO[25]+1)+1); k++){
	  for (l=0; l<(List_YOUSO[25]+1); l++){
	    for (m=0; m<List_YOUSO[24]; m++){
	      free(TmpHNA[i][j][k][l][m]);
	    }
	    free(TmpHNA[i][j][k][l]);
	  }
	  free(TmpHNA[i][j][k]);
	}
	free(TmpHNA[i][j]);
      }
      free(TmpHNA[i]);
    }
    free(TmpHNA);

    for (i=0; i<(List_YOUSO[25]+1); i++){
      for (j=0; j<List_YOUSO[24]; j++){
	for (k=0; k<(2*(List_YOUSO[25]+1)+1); k++){
	  for (l=0; l<(List_YOUSO[25]+1); l++){
	    for (m=0; m<List_YOUSO[24]; m++){
	      free(TmpHNAr[i][j][k][l][m]);
	    }
	    free(TmpHNAr[i][j][k][l]);
	  }
	  free(TmpHNAr[i][j][k]);
	}
	free(TmpHNAr[i][j]);
      }
      free(TmpHNAr[i]);
    }
    free(TmpHNAr);

    for (i=0; i<(List_YOUSO[25]+1); i++){
      for (j=0; j<List_YOUSO[24]; j++){
	for (k=0; k<(2*(List_YOUSO[25]+1)+1); k++){
	  for (l=0; l<(List_YOUSO[25]+1); l++){
	    for (m=0; m<List_YOUSO[24]; m++){
	      free(TmpHNAt[i][j][k][l][m]);
	    }
	    free(TmpHNAt[i][j][k][l]);
	  }
	  free(TmpHNAt[i][j][k]);
	}
	free(TmpHNAt[i][j]);
      }
      free(TmpHNAt[i]);
    }
    free(TmpHNAt);

    for (i=0; i<(List_YOUSO[25]+1); i++){
      for (j=0; j<List_YOUSO[24]; j++){
	for (k=0; k<(2*(List_YOUSO[25]+1)+1); k++){
	  for (l=0; l<(List_YOUSO[25]+1); l++){
	    for (m=0; m<List_YOUSO[24]; m++){
	      free(TmpHNAp[i][j][k][l][m]);
	    }
	    free(TmpHNAp[i][j][k][l]);
	  }
	  free(TmpHNAp[i][j][k]);
	}
	free(TmpHNAp[i][j]);
      }
      free(TmpHNAp[i]);
    }
    free(TmpHNAp);

    for (i=0; i<(2*(List_YOUSO[25]+1)+1); i++){
      free(CmatS0[i]);
    }
    free(CmatS0);

    for (i=0; i<(2*(List_YOUSO[25]+1)+1); i++){
      free(CmatSr[i]);
    }
    free(CmatSr);

    for (i=0; i<(2*(List_YOUSO[25]+1)+1); i++){
      free(CmatSt[i]);
    }
    free(CmatSt);

    for (i=0; i<(2*(List_YOUSO[25]+1)+1); i++){
      free(CmatSp[i]);
    }
    free(CmatSp);

#pragma omp flush(HVNA3)

  } /* #pragma omp parallel */

  /****************************************************
  freeing of arrays:
  ****************************************************/

  free(OneD2Mc_AN);
  free(OneD2h_AN);
  free(Gaunt_list);

  /* for time */
  dtime(&TEtime);
  time0 = TEtime - TStime;
  return time0;
} 


#define xmin  0.0
#define asize_lmax 30

#ifdef kcomp
static void Spherical_Bessel2( double x, int lmax, double *sb, double *dsb )
#elif defined (_NEC)
inline static void Spherical_Bessel2( double x, int lmax, double *sb, double *dsb )
#else 
inline void Spherical_Bessel2( double x, int lmax, double *sb, double *dsb )
#endif      
{
  int m,n,nmax;
  double tsb[asize_lmax+10];
  double invx,vsb0,vsb1,vsb2,vsbi;
  double j0,j1,j0p,j1p,sf,tmp,si,co,ix,ix2;
#ifndef _NEC
  if (x<0.0){
    printf("minus x is invalid for Spherical_Bessel\n");
    exit(0);
  }
#endif
  /* find an appropriate nmax */

  nmax = lmax + 3*x + 20;
  if (nmax<100) nmax = 100;
#ifndef _NEC
  if (asize_lmax<(lmax+1)){
    printf("asize_lmax should be larger than %d in Spherical_Bessel.c\n",lmax+1);
    exit(0);
  }
#endif
  /* if x is larger than xmin */

  if ( xmin < x ){

    invx = 1.0/x;

    /* initial values */

    vsb0 = 0.0;
    vsb1 = 1.0e-14;

    /* downward recurrence from nmax-2 to lmax+2 */
#ifdef _NEC
    #pragma _NEC novector
#endif
    for ( n=nmax-1; (lmax+2)<n; n-- ){

      vsb2 = (2.0*n + 1.0)*invx*vsb1 - vsb0;

      if (1.0e+250<vsb2){
        tmp = 1.0/vsb2;
        vsb2 *= tmp;
        vsb1 *= tmp;
      }

      vsbi = vsb0;
      vsb0 = vsb1;
      vsb1 = vsb2;
    }

    /* downward recurrence from lmax+1 to 0 */

    n = lmax + 3;
    tsb[n-1] = vsb1;
    tsb[n  ] = vsb0;
    tsb[n+1] = vsbi;

    tmp = tsb[n-1];
    tsb[n-1] /= tmp;
    tsb[n  ] /= tmp;
#ifdef _NEC
    #pragma _NEC novector
#endif
    for ( n=lmax+2; 0<n; n-- ){

      tsb[n-1] = (2.0*n + 1.0)*invx*tsb[n] - tsb[n+1];

      if (1.0e+250<tsb[n-1]){
        tmp = tsb[n-1];
#ifdef _NEC
    #pragma _NEC novector
#endif
        for (m=n-1; m<=lmax+1; m++){
          tsb[m] /= tmp;
        }
      }
    }

    /* normalization */

    si = sin(x);
    co = cos(x);
    ix = 1.0/x;
    ix2 = ix*ix;
    j0 = si*ix;
    j1 = si*ix*ix - co*ix;

    if (fabs(tsb[1])<fabs(tsb[0])) sf = j0/tsb[0];
    else                           sf = j1/tsb[1];

    /* tsb to sb */
#ifdef _NEC
    #pragma _NEC novector
#endif
    for ( n=0; n<=lmax+1; n++ ){
      sb[n] = tsb[n]*sf;
    }

    /* derivative of sb */

    dsb[0] = co*ix - si*ix*ix;
#ifdef _NEC
    #pragma _NEC novector
#endif
    for ( n=1; n<=lmax; n++ ){
      dsb[n] = ( (double)n*sb[n-1] - (double)(n+1.0)*sb[n+1] )/(2.0*(double)n + 1.0);
    }

  }

  /* if x is smaller than xmin */

  else {

    /* sb */
#ifdef _NEC
    #pragma _NEC novector
#endif
    for ( n=0; n<=lmax; n++ ){
      sb[n] = 0.0;
    }
    sb[0] = 1.0;

    /* derivative of sb */

    dsb[0] = 0.0;
#ifdef _NEC
    #pragma _NEC novector
#endif
    for ( n=1; n<=lmax; n++ ){
      dsb[n] = ( (double)n*sb[n-1] - (double)(n+1.0)*sb[n+1] )/(2.0*(double)n + 1.0);
    }
  }
}

#ifdef _NEC
#if 1
static inline void Spherical_Bessel2array(int x_size, double* x_array, double r, int Lmax_Four_Int, double *sb, double *dsb )
{
  double * B = (double*)malloc(sizeof(double)*(Lmax_Four_Int+3)*x_size);
  double * dB = (double*)malloc(sizeof(double)*(Lmax_Four_Int+3)*x_size);

  const int lmax = Lmax_Four_Int;

#define _NEC_NMAX_V1
#ifdef _NEC_NMAX_V1
  int nmax = 100;
  
  #pragma _NEC ivdep
  #pragma _NEC vector
  for (int i=0; i<x_size; i++){
    const double x = GL_NormK[i] * r;
    const int nn = lmax + 3*x + 20;
    if( nmax < nn) {
      nmax = nn;
    }
  }
#endif  

  
  double* vsb0 = (double*)malloc(sizeof(double)*x_size);
  double* vsb1 = (double*)malloc(sizeof(double)*x_size);
  double* vsb2 = (double*)malloc(sizeof(double)*x_size);
  double* vsbi = (double*)malloc(sizeof(double)*x_size);
  double* sf = (double*)malloc(sizeof(double)*x_size);
  
  double* tsb = (double*)malloc(sizeof(double)*(asize_lmax+10)*x_size);

  /* initial values */
  #pragma _NEC vector
  for (int i=0; i<x_size; i++){
      vsb0[i] = 0.0;
      vsb1[i] = 1.0e-14;
  }  

  /* downward recurrence from nmax-2 to lmax+2 */
#ifdef _NEC_NMAX_V1
  #pragma _NEC novector
  for ( int n=nmax-1; (lmax+2)<n; n-- ){

    #pragma _NEC vector
    for (int i=0; i<x_size; i++){
      const double x = GL_NormK[i] * r;
      const double invx = 1.0/x;

      vsb2[i] = (2.0*n + 1.0)*invx*vsb1[i] - vsb0[i];

      if (1.0e+250<vsb2[i]){
        const double tmp = 1.0/vsb2[i];
        vsb2[i] *= tmp;
        vsb1[i] *= tmp;
      }

      vsbi[i] = vsb0[i];
      vsb0[i] = vsb1[i];
      vsb1[i] = vsb2[i];
    }
  }
#else

  #pragma _NEC novector
  for (int i=0; i<x_size; i++){
    const double x = GL_NormK[i] * r;
    const double invx = 1.0/x;
    const int nn = lmax + 3*x + 20;
    const int nmax = (nn > 100) ? nn : 100;

    #pragma _NEC ivdep
    #pragma _NEC vector
    for ( int n=nmax-1; (lmax+2)<n; n-- ){
      vsb2[i] = (2.0*n + 1.0)*invx*vsb1[i] - vsb0[i];

      if (1.0e+250<vsb2[i]){
        const double tmp = 1.0/vsb2[i];
        vsb2[i] *= tmp;
        vsb1[i] *= tmp;
      }

      vsbi[i] = vsb0[i];
      vsb0[i] = vsb1[i];
      vsb1[i] = vsb2[i];
    }
  }
#endif  

  
  #pragma _NEC vector
  for (int i=0; i<x_size; i++){
    const double x = GL_NormK[i] * r;
    const double invx = 1.0/x;

      /* downward recurrence from lmax+1 to 0 */
      {
        const int n = lmax + 3;
        tsb[(n-1)*x_size + i] = vsb1[i];
        tsb[n*x_size + i  ] = vsb0[i];
        tsb[(n+1)*x_size + i] = vsbi[i];

        //const double tmp = tsb[n-1];
        tsb[(n-1)*x_size + i] /= vsb1[i];
        tsb[n*x_size + i  ] /= vsb1[i];
      }
  }


  #pragma _NEC novector
  for ( int n=lmax+2; 0<n; n-- ){

    #pragma _NEC vector
    for (int i=0; i<x_size; i++){
      const double x = GL_NormK[i] * r;
      const double invx = 1.0/x;

      tsb[(n-1)*x_size + i] = (2.0*n + 1.0)*invx*tsb[n*x_size + i] - tsb[(n+1)*x_size + i];

#ifdef _NEC_BESSEL_WITH_SAFTY 
      if (1.0e+250<tsb[(n-1)*x_size + i]){
        const double tmp = tsb[(n-1)*x_size + i];
    #pragma _NEC novector
        for (int m=n-1; m<=lmax+1; m++){
          tsb[m*x_size + i] /= tmp;
        }
      }
#endif      
    }
  }

  /* normalization */
  #pragma _NEC vector
  for (int i=0; i<x_size; i++){
    const double x = GL_NormK[i] * r;

      const double si = sin(x);
      const double co = cos(x);
      const double ix = 1.0/x;
      const double ix2 = ix*ix;
      const double j0 = si*ix;
      const double j1 = si*ix*ix - co*ix;
      
      if(xmin < x){
        if (fabs(tsb[1*x_size + i])<fabs(tsb[0*x_size + i])) sf[i] = j0/tsb[0*x_size + i];
        else                           sf[i] = j1/tsb[1*x_size + i];

        dB[0*x_size+i] = co*ix - si*ix*ix;

      }else{
        sf[i] = 0.0;
        dB[0*x_size+i] = 0.0;
      }
  }

  /* tsb to B */
  #pragma _NEC novector
  for (int n=0; n<=lmax+1; n++ ){
    #pragma _NEC vector
    for (int i=0; i<x_size; i++){
      B[n*x_size + i] = tsb[n*x_size + i]*sf[i];
    }
  }

  #pragma _NEC vector
  for (int i=0; i<x_size; i++){
    const double x = GL_NormK[i] * r;
    if(!(xmin < x)){
      B[0*x_size+i] = 1.0; 
    }
  }

    /* derivative of B */
  
  #pragma _NEC novector
  for (int n=1; n<=lmax; n++ ){
    #pragma _NEC vector
    for (int i=0; i<x_size; i++){  
      dB[n*x_size+i] = ( (double)n*B[(n-1)*x_size+i] - (double)(n+1.0)*B[(n+1)*x_size+i] )/(2.0*(double)n + 1.0);
    }
  }

 
  
  for(int k = 0; k < (Lmax_Four_Int+3); ++k){
    #pragma _NEC vector
    for (int i=0; i<x_size; i++){
    /****************************************************************/    
      sb[i + k * x_size] = B[i+k*x_size];
      dsb[i + k * x_size] = dB[i+k*x_size];
    }
  }

  free(B);
  free(dB);
  free(vsb0);
  free(vsb1);
  free(vsb2);
  free(vsbi);
  free(sf);
  free(tsb);
}
#else
static inline void Spherical_Bessel2array(int x_size, double* x_array, double coef, int lmax, double *sb, double *dsb )
{
  /*
  int m,n,nmax;
  double tsb[asize_lmax+10];
  double invx,vsb0,vsb1,vsb2,vsbi;
  double j0,j1,j0p,j1p,sf,tmp,si,co,ix,ix2;
  */
  double* tsb = (double*)malloc(sizeof(double)*(asize_lmax+10) * x_size);
  double* vbs = (double*)malloc(sizeof(double)*x_size * 4);
  /* find an appropriate nmax */

  //int nmax = lmax + 3*x + 20;
  //if (nmax<100) nmax = 100;
  int nmax = 100;
#ifndef _NEC
  if (asize_lmax<(lmax+1)){
    printf("asize_lmax should be larger than %d in Spherical_Bessel.c\n",lmax+1);
    exit(0);
  }
#endif
  /* if x is larger than xmin */
  
    //if ( xmin < x )
    {


    /* initial values */
    #pragma _NEC vector
    for(int i = 0; i < x_size; ++i ){

      vbs[i*4 + 0] = 0.0;
      vbs[i*4 + 1] = 1.0e-14;
    }

    /* downward recurrence from nmax-2 to lmax+2 */
    for (int n=nmax-1; (lmax+2)<n; n-- ){

      #pragma _NEC vector
      for(int i = 0; i < x_size ;++i){

        double x = x_array[i] * coef;
        double invx = 1.0/x;

        vbs[i*4 + 2] = (2.0*n + 1.0)*invx*vbs[i*4 + 1] - vbs[i*4 + 0];

        if (1.0e+250<vbs[i*4 + 2]){
          double tmp = 1.0/vbs[i*4 + 2];
          vbs[i*4 + 2] *= tmp;
          vbs[i*4 + 1] *= tmp;
        }

        vbs[i*4 + 3] = vbs[i*4 + 0];
        vbs[i*4 + 0] = vbs[i*4 + 1];
        vbs[i*4 + 1] = vbs[i*4 + 2];
      }
    }

    /* downward recurrence from lmax+1 to 0 */

    int n = lmax + 3;
    #pragma _NEC vector
    for(int i = 0; i < x_size ; ++i){
      tsb[(n-1)*x_size + i] = vbs[i*4 + 1];
      tsb[(n  )*x_size + i] = vbs[i*4 + 0];
      tsb[(n+1)*x_size + i] = vbs[i*4 + 3];
      
      tsb[(n-1)*x_size + i] /= vbs[i*4 + 1];
      tsb[(n  )*x_size + i] /= vbs[i*4 + 1];
    }

    for (int n=lmax+2; 0<n; n-- ){

#pragma _NEC vector
      for(int i = 0; i < x_size ; ++i){
        double x = x_array[i] * coef;
        double invx = 1.0/x;
        tsb[(n-1)*x_size + i] = (2.0*n + 1.0)*invx*tsb[(n)*x_size + i] - tsb[(n+1)*x_size + i];
      }
#if 0
/* this routine is probably for error, but this is difficult for vectorization*/
      if (1.0e+250<tsb[(n-1)*x_size + i]){
        tmp = tsb[(n-1)*x_size + i];
#ifdef _NEC
    #pragma _NEC novector
#endif
        for (int m=n-1; m<=lmax+1; m++){
          tsb[m] /= tmp;
        }
      }
#endif      
    }

    /* normalization */

#pragma _NEC vector
    for(int i = 0; i < x_size ; ++i){
      double x = x_array[i] * coef;
      double si = sin(x);
      double co = cos(x);
      double ix = 1.0/x;
      double ix2 = ix*ix;
      double j0 = si*ix;
      double j1 = si*ix*ix - co*ix;

      double sf;
      if (x <= xmin) sf = 0.0;
      else if (fabs(tsb[x_size + i])<fabs(tsb[i])) sf = j0/tsb[i];
      else                           sf = j1/tsb[x_size + i];

      vbs[i*4 + 0] = sf; //temporary//

      dsb[i] = (x > xmin) ? co*ix - si*ix*ix : 0.0;
    }

    /* tsb to sb */
    //for (int n=0; n<=lmax+1; n++ ){ //this is probably bug.
    for (int n=0; n<=lmax; n++ ){      
      #pragma _NEC vector
      for(int i = 0; i < x_size ; ++i){
        double sf = vbs[i*4 + 0];
        sb[n*x_size + i] = tsb[n*x_size + i]*sf;
      }
    }

    /* derivative of sb */

    for ( n=1; n<=lmax; n++ ){
      #pragma _NEC vector
      for(int i = 0; i < x_size ; ++i){
        dsb[n*x_size + i] = ( (double)n*sb[(n-1)*x_size + i] - (double)(n+1.0)*sb[(n+1)*x_size + i] )/(2.0*(double)n + 1.0);
      }
    }

  }

  free(tsb);
  free(vbs);
}
#endif
#endif

/*
inline dcomplex** Allocate2D_dcomplex(int size_1, int size_2)
{ 
  int i, j;

  dcomplex** buffer = (dcomplex**)malloc(sizeof(dcomplex*)*size_1);
  buffer[0] = (dcomplex*)malloc(sizeof(dcomplex)*size_1*size_2);

  for (i=0; i<size_1; i++){
    buffer[i] = buffer[0] + i * size_2;
    for (j=0; j<size_2; j++){
      buffer[i][j] = Complex(0.0,0.0);
    }
  }

  return buffer;
}

inline void Free2D_dcomplex(dcomplex** buffer)
{ 
  free(buffer[0]);
  free(buffer);
}


inline dcomplex**** Allocate4D_dcomplex(int size_1, int size_2, int size_3, int size_4)
{ 
  int i, j, k, l;

  dcomplex**** buffer = (dcomplex****)malloc(sizeof(dcomplex***)*size_1);
  buffer[0] = (dcomplex***)malloc(sizeof(dcomplex**)*size_1*size_2);
  buffer[0][0] = (dcomplex**)malloc(sizeof(dcomplex*)*size_1*size_2*size_3);
  buffer[0][0][0] = (dcomplex*)malloc(sizeof(dcomplex)*size_1*size_2*size_3*size_4);
		
  for (i=0; i<size_1; i++){
    buffer[i] = buffer[0] + i * size_2;
    for (j=0; j<size_2; j++){
      buffer[i][j] = buffer[0][0] + (i * size_2 + j) * size_3;
      for (k=0; k<size_3; k++){
	buffer[i][j][k] = buffer[0][0][0] + ((i * size_2 + j) * size_3 + k) * size_4;
	for (l=0; l<size_4; l++){
	  buffer[i][j][k][l] = Complex(0.0,0.0);
	}
      }
    }
  }

  return buffer;
}

inline void Free4D_dcomplex(dcomplex**** buffer)
{ 
  free(buffer[0][0][0]);
  free(buffer[0][0]);
  free(buffer[0]);
  free(buffer);
}
*/


