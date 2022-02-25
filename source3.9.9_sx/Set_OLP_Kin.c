/**********************************************************************
  Set_OLP_Kin.c:

     Set_OLP_Kin.c is a subroutine to calculate the overlap matrix
     and the matrix for the kinetic operator in momentum space.

  Log of Set_OLP_Kin.c:

     15/Oct./2002  Released by T.Ozaki
     25/Nov./2014  Memory allocation modified by A.M. Ito (AITUNE)

***********************************************************************/
 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
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

#ifdef kcomp
dcomplex****** Allocate6D_dcomplex(int size_1, int size_2, int size_3, 
                                          int size_4, int size_5, int size_6);
double**** Allocate4D_double(int size_1, int size_2, int size_3, int size_4);
dcomplex** Allocate2D_dcomplex(int size_1, int size_2);
double** Allocate2D_double(int size_1, int size_2);
dcomplex**** Allocate4D_dcomplex(int size_1, int size_2, int size_3, int size_4);
void Free6D_dcomplex(dcomplex****** buffer);
void Free4D_double(double**** buffer);
void Free2D_dcomplex(dcomplex** buffer);
void Free2D_double(double** buffer);
void Free4D_dcomplex(dcomplex**** buffer);
#else
static inline dcomplex****** Allocate6D_dcomplex(int size_1, int size_2, int size_3, 
                                          int size_4, int size_5, int size_6);
static inline double**** Allocate4D_double(int size_1, int size_2, int size_3, int size_4);
static inline dcomplex** Allocate2D_dcomplex(int size_1, int size_2);
static inline double** Allocate2D_double(int size_1, int size_2);
static inline dcomplex**** Allocate4D_dcomplex(int size_1, int size_2, int size_3, int size_4);
void Free6D_dcomplex(dcomplex****** buffer);
void Free2D_double(double** buffer);
void Free4D_double(double**** buffer);
void Free4D_dcomplex(dcomplex**** buffer);
void Free2D_dcomplex(dcomplex** buffer);
#endif

#ifdef _NEC
#include "Gaunt.h"
#define _NEC_OLP_MG
#endif


dcomplex Conjg_inline(dcomplex z)
{
  dcomplex c;
  c.r =  z.r;
  c.i = -z.i;
  return c;
}

dcomplex Complex_inline(double re, double im)
{
  dcomplex c;
  c.r = re;
  c.i = im;
  return c;
}

dcomplex Cadd_inline(dcomplex a, dcomplex b)
{
  dcomplex c;
  c.r = a.r + b.r;
  c.i = a.i + b.i;
  return c;
}

dcomplex Csub_inline(dcomplex a, dcomplex b)
{
  dcomplex c;
  c.r = a.r - b.r;
  c.i = a.i - b.i;
  return c;
}

dcomplex Cmul_inline(dcomplex a, dcomplex b)
{
  dcomplex c;
  c.r = a.r*b.r - a.i*b.i;
  c.i = a.i*b.r + a.r*b.i;
  return c;
}





double Set_OLP_Kin(double *****OLP, double *****H0)
{
  /****************************************************
          Evaluate overlap and kinetic integrals
                 in the momentum space
  ****************************************************/
  static int firsttime=1;
  int size_SumS0,size_TmpOLP;
  double time0;
  double TStime,TEtime;
  int numprocs,myid;
  int Mc_AN,Gc_AN,h_AN;
  int OneD_Nloop;
  int *OneD2Mc_AN,*OneD2h_AN;

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  dtime(&TStime);

  /****************************************************
   MPI_Barrier
  ****************************************************/

  MPI_Barrier(mpi_comm_level1);

  /* PrintMemory */

  if (firsttime) {

    size_SumS0  = (List_YOUSO[25]+1)*List_YOUSO[24]*(List_YOUSO[25]+1)*List_YOUSO[24];
    size_TmpOLP = (List_YOUSO[25]+1)*List_YOUSO[24]*(2*(List_YOUSO[25]+1)+1)*
                  (List_YOUSO[25]+1)*List_YOUSO[24]*(2*(List_YOUSO[25]+1)+1);
 
    PrintMemory("Set_OLP_Kin: SumS0",  sizeof(double)*size_SumS0,NULL);
    PrintMemory("Set_OLP_Kin: SumK0",  sizeof(double)*size_SumS0,NULL);
    PrintMemory("Set_OLP_Kin: SumSr0", sizeof(double)*size_SumS0,NULL);
    PrintMemory("Set_OLP_Kin: SumKr0", sizeof(double)*size_SumS0,NULL);
    PrintMemory("Set_OLP_Kin: TmpOLP", sizeof(dcomplex)*size_TmpOLP,NULL);
    PrintMemory("Set_OLP_Kin: TmpOLPr",sizeof(dcomplex)*size_TmpOLP,NULL);
    PrintMemory("Set_OLP_Kin: TmpOLPt",sizeof(dcomplex)*size_TmpOLP,NULL);
    PrintMemory("Set_OLP_Kin: TmpOLPp",sizeof(dcomplex)*size_TmpOLP,NULL);
    PrintMemory("Set_OLP_Kin: TmpKin", sizeof(dcomplex)*size_TmpOLP,NULL);
    PrintMemory("Set_OLP_Kin: TmpKinr",sizeof(dcomplex)*size_TmpOLP,NULL);
    PrintMemory("Set_OLP_Kin: TmpKint",sizeof(dcomplex)*size_TmpOLP,NULL);
    PrintMemory("Set_OLP_Kin: TmpKinp",sizeof(dcomplex)*size_TmpOLP,NULL);
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
  for (int Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
    int Gc_AN = M2G[Mc_AN];    
    for (int h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
      OneD2Mc_AN[OneD_Nloop] = Mc_AN; 
      OneD2h_AN[OneD_Nloop]  = h_AN; 
      OneD_Nloop++;
    }
  }

#ifdef _NEC
ftrace_region_begin("OLP_MG1");

  int max_L0 = 0;
  for (int spe=0; spe<SpeciesNum; spe++){
    if(max_L0 < Spe_MaxL_Basis[spe]) max_L0 = Spe_MaxL_Basis[spe];
  }
  const int ai_L0_LIMIT = max_L0 + 1;
  const int ai_L0_LIMIT_SQ = ai_L0_LIMIT * ai_L0_LIMIT;
  const int ai_L1_LIMIT = ai_L0_LIMIT;
  const int ai_L1_LIMIT_SQ = ai_L1_LIMIT*ai_L1_LIMIT;
  const int ai_LL_LIMIT = max_L0*2+1;
  const int ai_LL_LIMIT_SQ = ai_LL_LIMIT*ai_LL_LIMIT;
  double* Gaunt_list = MakeGaunt(ai_L0_LIMIT_SQ, ai_L1_LIMIT_SQ, ai_LL_LIMIT_SQ);
ftrace_region_end("OLP_MG1");
#endif

  /* OpenMP */

#pragma omp parallel
  {

    int Nloop;
    int OMPID,Nthrds,Nprocs;
    int Mc_AN,h_AN,Gc_AN,Cwan;
    int Gh_AN,Rnh,Hwan;
    int Ls;    /*,L0,Mul0,L1,Mul1,M0,M1;*/
    int Lmax_Four_Int;
    int i,j,k,l/*,m,p*/;
    //int num0,num1; 

    double Stime_atom,Etime_atom;
    double dx,dy,dz;
    double S_coordinate[3];
    double theta,phi,h;
    double Bessel_Pro0,Bessel_Pro1;
    double tmp0,tmp1,tmp2,tmp3,tmp4;
    double siT,coT,siP,coP;
    double kmin,kmax,Sk,Dk,r; 
    double sj,sjp,coe0,coe1; 
    double Normk,Normk2;
    double gant,SH[2],dSHt[2],dSHp[2];
    double **SphB,**SphBp;
    double *tmp_SphB,*tmp_SphBp;
    double **SumS0;
    double **SumK0;
    double **SumSr0;
    double **SumKr0;

    dcomplex CsumS_Lx,CsumS_Ly,CsumS_Lz;
    dcomplex CsumS0,CsumSr,CsumSt,CsumSp;
    dcomplex CsumK0,CsumKr,CsumKt,CsumKp;
    dcomplex Ctmp0,Ctmp1,Ctmp2,Cpow;
    dcomplex CY,CYt,CYp,CY1,CYt1,CYp1;
    dcomplex **TmpOLP;
    dcomplex **TmpOLPr;
    dcomplex **TmpOLPt;
    dcomplex **TmpOLPp;
    dcomplex **TmpKin;
    dcomplex **TmpKinr;
    dcomplex **TmpKint;
    dcomplex **TmpKinp;

    dcomplex *CmatS0;
    dcomplex *CmatSr;
    dcomplex *CmatSt;
    dcomplex *CmatSp;
    dcomplex *CmatK0;
    dcomplex *CmatKr;
    dcomplex *CmatKt;
    dcomplex *CmatKp;

    /****************************************************************
                          allocation of arrays:
    ****************************************************************/

    TmpOLP  = Allocate2D_dcomplex((List_YOUSO[25]+1)*List_YOUSO[24] * (2*(List_YOUSO[25]+1)+1), (List_YOUSO[25]+1)* List_YOUSO[24] * (2*(List_YOUSO[25]+1)+1));
    TmpOLPr = Allocate2D_dcomplex((List_YOUSO[25]+1)*List_YOUSO[24] * (2*(List_YOUSO[25]+1)+1), (List_YOUSO[25]+1)* List_YOUSO[24] * (2*(List_YOUSO[25]+1)+1));
    TmpOLPt = Allocate2D_dcomplex((List_YOUSO[25]+1)*List_YOUSO[24] * (2*(List_YOUSO[25]+1)+1), (List_YOUSO[25]+1)* List_YOUSO[24] * (2*(List_YOUSO[25]+1)+1));
    TmpOLPp = Allocate2D_dcomplex((List_YOUSO[25]+1)*List_YOUSO[24] * (2*(List_YOUSO[25]+1)+1), (List_YOUSO[25]+1)* List_YOUSO[24] * (2*(List_YOUSO[25]+1)+1));
    TmpKin  = Allocate2D_dcomplex((List_YOUSO[25]+1)*List_YOUSO[24] * (2*(List_YOUSO[25]+1)+1), (List_YOUSO[25]+1)* List_YOUSO[24] * (2*(List_YOUSO[25]+1)+1));
    TmpKinr = Allocate2D_dcomplex((List_YOUSO[25]+1)*List_YOUSO[24] * (2*(List_YOUSO[25]+1)+1), (List_YOUSO[25]+1)* List_YOUSO[24] * (2*(List_YOUSO[25]+1)+1));
    TmpKint = Allocate2D_dcomplex((List_YOUSO[25]+1)*List_YOUSO[24] * (2*(List_YOUSO[25]+1)+1), (List_YOUSO[25]+1)* List_YOUSO[24] * (2*(List_YOUSO[25]+1)+1));
    TmpKinp = Allocate2D_dcomplex((List_YOUSO[25]+1)*List_YOUSO[24] * (2*(List_YOUSO[25]+1)+1), (List_YOUSO[25]+1)* List_YOUSO[24] * (2*(List_YOUSO[25]+1)+1));

    SumS0  = Allocate2D_double((List_YOUSO[25]+1)* List_YOUSO[24], (List_YOUSO[25]+1)* List_YOUSO[24]);
    SumK0  = Allocate2D_double((List_YOUSO[25]+1)* List_YOUSO[24], (List_YOUSO[25]+1)* List_YOUSO[24]);
    SumSr0 = Allocate2D_double((List_YOUSO[25]+1)* List_YOUSO[24], (List_YOUSO[25]+1)* List_YOUSO[24]);
    SumKr0 = Allocate2D_double((List_YOUSO[25]+1)* List_YOUSO[24], (List_YOUSO[25]+1)* List_YOUSO[24]);
	
    CmatS0 = (dcomplex*)malloc(sizeof(dcomplex*)*(2*(List_YOUSO[25]+1)+1)*(List_YOUSO[25]+1)* List_YOUSO[24]);
    CmatSr = (dcomplex*)malloc(sizeof(dcomplex*)*(2*(List_YOUSO[25]+1)+1)*(List_YOUSO[25]+1)* List_YOUSO[24]);
    CmatSt = (dcomplex*)malloc(sizeof(dcomplex*)*(2*(List_YOUSO[25]+1)+1)*(List_YOUSO[25]+1)* List_YOUSO[24]);
    CmatSp = (dcomplex*)malloc(sizeof(dcomplex*)*(2*(List_YOUSO[25]+1)+1)*(List_YOUSO[25]+1)* List_YOUSO[24]);
    CmatK0 = (dcomplex*)malloc(sizeof(dcomplex*)*(2*(List_YOUSO[25]+1)+1)*(List_YOUSO[25]+1)* List_YOUSO[24]);
    CmatKr = (dcomplex*)malloc(sizeof(dcomplex*)*(2*(List_YOUSO[25]+1)+1)*(List_YOUSO[25]+1)* List_YOUSO[24]);
    CmatKt = (dcomplex*)malloc(sizeof(dcomplex*)*(2*(List_YOUSO[25]+1)+1)*(List_YOUSO[25]+1)* List_YOUSO[24]);
    CmatKp = (dcomplex*)malloc(sizeof(dcomplex*)*(2*(List_YOUSO[25]+1)+1)*(List_YOUSO[25]+1)* List_YOUSO[24]);
  /*
    CmatS0 = Allocate2D_dcomplex((2*(List_YOUSO[25]+1)+1), (2*(List_YOUSO[25]+1)+1));
    CmatSr = Allocate2D_dcomplex((2*(List_YOUSO[25]+1)+1), (2*(List_YOUSO[25]+1)+1));
    CmatSt = Allocate2D_dcomplex((2*(List_YOUSO[25]+1)+1), (2*(List_YOUSO[25]+1)+1));
    CmatSp = Allocate2D_dcomplex((2*(List_YOUSO[25]+1)+1), (2*(List_YOUSO[25]+1)+1));
    CmatK0 = Allocate2D_dcomplex((2*(List_YOUSO[25]+1)+1), (2*(List_YOUSO[25]+1)+1));
    CmatKr = Allocate2D_dcomplex((2*(List_YOUSO[25]+1)+1), (2*(List_YOUSO[25]+1)+1));
    CmatKt = Allocate2D_dcomplex((2*(List_YOUSO[25]+1)+1), (2*(List_YOUSO[25]+1)+1));
    CmatKp = Allocate2D_dcomplex((2*(List_YOUSO[25]+1)+1), (2*(List_YOUSO[25]+1)+1));
 */
//AITUNE//
    //not used// double** ai_Gaunt = Allocate2D_double((List_YOUSO[25]+1) *(2*(List_YOUSO[25]+1)+1), (List_YOUSO[25]+1) *(2*(List_YOUSO[25]+1)+1));

    double* tmp_SphB_D  = (double*)malloc(sizeof(double)*(2*(List_YOUSO[25]+1)+1)*(OneD_Grid+1));
    double* tmp_SphBp_D = (double*)malloc(sizeof(double)*(2*(List_YOUSO[25]+1)+1)*(OneD_Grid+1));
    
    /* for RF_BesselF */
    int* normk_indexes = (int*)malloc(sizeof(int)*(OneD_Grid+1));
    {
      const double kmin = Radial_kmin;
      const double kmax = PAO_Nkmax;      
      const double h = (kmax - kmin)/(double)OneD_Grid;
      for (int i=0; i<=OneD_Grid; i++){
        const double Normk = kmin + (double)i*h;
        normk_indexes[i] = GetNormKIndex(Normk);
      }
    }
    

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


      /* AITUNE */
      int* ai_L0;
      int* ai_offset0;
      int* ai_vec_L0;
      int* ai_vec_num0;
      int* ai_L1;
      int* ai_offset1;
      int* ai_vec_L1;
      int* ai_vec_num1;
      int ai_max_LMul0, ai_max_LMul1;
      const int ai_len_vec_LM0 = Spe_Total_NO[Cwan];
      const int ai_len_vec_LM1 = Spe_Total_NO[Hwan];
      {
        ai_max_LMul0 = 0;
        
        for (int L0=0; L0<=Spe_MaxL_Basis[Cwan]; L0++){
          ai_max_LMul0 += Spe_Num_Basis[Cwan][L0];
        }
        ai_L0 = (int*)malloc(sizeof(int*)*ai_max_LMul0);
        ai_offset0 = (int*)malloc(sizeof(int*)*ai_max_LMul0);
        
        ai_vec_L0 = (int*)malloc(sizeof(int*)*ai_len_vec_LM0);
        ai_vec_num0 = (int*)malloc(sizeof(int*)*ai_len_vec_LM0);

        int n0 = 0;
        int offset0 = 0;
        for (int L0=0; L0<=Spe_MaxL_Basis[Cwan]; L0++){
	        for (int Mul0=0; Mul0<Spe_Num_Basis[Cwan][L0]; Mul0++){
            ai_L0[n0] = L0;
            ai_offset0[n0] = offset0;
            
            for (int L0M0=0; L0M0<=2*L0; L0M0++){
              ai_vec_L0[offset0 + L0M0] = L0;
              ai_vec_num0[offset0 + L0M0] = offset0;
            }
            n0++;
            offset0 += 2 * L0 + 1;
          }
        }
      
        ai_max_LMul1 = 0;
        for (int L1=0; L1<=Spe_MaxL_Basis[Hwan]; L1++){
          ai_max_LMul1 += Spe_Num_Basis[Hwan][L1];
        }
        ai_L1 = (int*)malloc(sizeof(int*)*ai_max_LMul1);
        ai_offset1 = (int*)malloc(sizeof(int*)*ai_max_LMul1);
        
        ai_vec_L1 = (int*)malloc(sizeof(int*)*ai_len_vec_LM1);
        ai_vec_num1 = (int*)malloc(sizeof(int*)*ai_len_vec_LM1);
        
        int n1 = 0;
        int offset1 = 0;
        for (int L1=0; L1<=Spe_MaxL_Basis[Hwan]; L1++){
	        for (int Mul1=0; Mul1<Spe_Num_Basis[Hwan][L1]; Mul1++){
            ai_L1[n1] = L1;
            ai_offset1[n1] = offset1;
            
            for (int L1M1=0; L1M1<=2*L1; L1M1++){
              ai_vec_L1[offset1 + L1M1] = L1;
              ai_vec_num1[offset1 + L1M1] = offset1;
            }
            n1++;
            offset1 += 2 * L1 + 1;
          }
        }
      }


      /****************************************************
          Overlap and the derivative
              \int RL(k)*RL'(k)*jl(k*R) k^2 dk^3,
              \int RL(k)*RL'(k)*j'l(k*R) k^3 dk^3

          Kinetic 
              \int RL(k)*RL'(k)*jl(k*R) k^4 dk^3, 
              \int RL(k)*RL'(k)*j'l(k*R) k^5 dk^3 
      ****************************************************/

      kmin = Radial_kmin;
      kmax = PAO_Nkmax;
      Sk = kmax + kmin;
      Dk = kmax - kmin;

      for (int vi0=0; vi0<ai_len_vec_LM0; vi0++){
	      for (int vi1=0; vi1<ai_len_vec_LM1; ++vi1){
	      
          TmpOLP[vi0][vi1]  = Complex_inline(0.0,0.0);
          TmpOLPr[vi0][vi1]  = Complex_inline(0.0,0.0);
          TmpOLPt[vi0][vi1]  = Complex_inline(0.0,0.0);
          TmpOLPp[vi0][vi1]  = Complex_inline(0.0,0.0);
          TmpKin[vi0][vi1]  = Complex_inline(0.0,0.0);
          TmpKinr[vi0][vi1]  = Complex_inline(0.0,0.0);
          TmpKint[vi0][vi1]  = Complex_inline(0.0,0.0);
          TmpKinp[vi0][vi1]  = Complex_inline(0.0,0.0);


		    }
	    }
	

      if (Spe_MaxL_Basis[Cwan]<Spe_MaxL_Basis[Hwan])
        Lmax_Four_Int = 2*Spe_MaxL_Basis[Hwan];
      else 
        Lmax_Four_Int = 2*Spe_MaxL_Basis[Cwan];



      /* calculate SphB and SphBp */

      h = (kmax - kmin)/(double)OneD_Grid;


#pragma _NEC vector
      for (int i=0; i<=OneD_Grid; i++){
        Normk = kmin + (double)i*h;
        Spherical_Bessel(Normk*r,Lmax_Four_Int,&(tmp_SphB_D[i*(Lmax_Four_Int+3)]),&(tmp_SphBp_D[i*(Lmax_Four_Int+3)]));
        /*
        #pragma _NEC novector
        for(l=0; l<=Lmax_Four_Int; l++){ 
          SphB[l][i]  = tmp_SphB_D[i*(Lmax_Four_Int+3) + l]; 
          SphBp[l][i] = tmp_SphBp_D[i*(Lmax_Four_Int+3) + l]; 
	      }
        */
      }

/*
      free(tmp_SphB);
      free(tmp_SphBp);
*/      
      


      /* l loop */

      for(int l=0; l<=Lmax_Four_Int; l++){
        for (int LMul0=0; LMul0<ai_max_LMul0; LMul0++){
          for (int LMul1=0; LMul1<ai_max_LMul1; LMul1++){
            SumS0[LMul0][LMul1]  = 0.0;
            SumK0[LMul0][LMul1]  = 0.0;
            SumSr0[LMul0][LMul1] = 0.0;
            SumKr0[LMul0][LMul1] = 0.0;
	        }
	      }
        //printf("SumS0\n"); fflush(stdout);

        h = (kmax - kmin)/(double)OneD_Grid;

        int LMul0=0;
        for (int L0=0; L0<=Spe_MaxL_Basis[Cwan]; L0++){
          for (int Mul0=0; Mul0<Spe_Num_Basis[Cwan][L0]; Mul0++){
              int LMul1=0;
              for (int L1=0; L1<=Spe_MaxL_Basis[Hwan]; L1++){
                for (int Mul1=0; Mul1<Spe_Num_Basis[Hwan][L1]; Mul1++){


	      for (i=0; i<=OneD_Grid; i++){

          if (i==0 || i==OneD_Grid) coe0 = 0.50;
          else                      coe0 = 1.00;

	        Normk = kmin + (double)i*h;
          Normk2 = Normk*Normk;
          const int n_i = normk_indexes[i];
/*
          sj  =  SphB[l][i];
          sjp = SphBp[l][i];
          */
         
          sj  = tmp_SphB_D[i*(Lmax_Four_Int+3) + l]; 
          sjp = tmp_SphBp_D[i*(Lmax_Four_Int+3) + l]; 

/*
	        int LMul0=0;
          for (int L0=0; L0<=Spe_MaxL_Basis[Cwan]; L0++){
	          for (int Mul0=0; Mul0<Spe_Num_Basis[Cwan][L0]; Mul0++){
 */     

	            //Bessel_Pro0 = RF_BesselF(Cwan,L0,Mul0,Normk);
              Bessel_Pro0 = RF_BesselF_ByIndex(Cwan,L0,Mul0,Normk, n_i);

              tmp0 = coe0*h*Normk2*Bessel_Pro0;
              tmp1 = tmp0*sj;
              tmp2 = tmp0*Normk*sjp;
/*
              int LMul1=0;
              for (int L1=0; L1<=Spe_MaxL_Basis[Hwan]; L1++){
                for (int Mul1=0; Mul1<Spe_Num_Basis[Hwan][L1]; Mul1++){
 */         
                  //Bessel_Pro1 = RF_BesselF(Hwan,L1,Mul1,Normk);
                  Bessel_Pro1 = RF_BesselF_ByIndex(Hwan,L1,Mul1,Normk, n_i);

                  tmp3 = tmp1*Bessel_Pro1;
                  tmp4 = tmp2*Bessel_Pro1;

                  SumS0[LMul0][LMul1] += tmp3;
                  SumK0[LMul0][LMul1] += tmp3*Normk2;

                  SumSr0[LMul0][LMul1] += tmp4;
                  SumKr0[LMul0][LMul1] += tmp4*Normk2;
                  /*
         
		            }
	            }
         
	          }
	        }
          */
	      }
                LMul1++;
		            }
	            }
              LMul0++;
	          }
	        }
        //printf("Bessel\n"); fflush(stdout);
  

        if (h_AN==0){ 
	        for (int LMul0=0; LMul0<ai_max_LMul0; LMul0++){
	          for (int LMul1=0; LMul1<ai_max_LMul1; LMul1++){
		          SumSr0[LMul0][LMul1] = 0.0;
		          SumKr0[LMul0][LMul1] = 0.0;
			      }
	        }
	      }
        //printf("SumSr0\n"); fflush(stdout);
        /****************************************************
          For overlap and the derivative,
          sum_m 8*(-i)^{-L0+L1+1}*
                C_{L0,-M0,L1,M1,l,m}*Y_{lm}
                \int RL(k)*RL'(k)*jl(k*R) k^2 dk^3,

          For kinetic,
          sum_m 4*(-i)^{-L0+L1+1}*
                C_{L0,-M0,L1,M1,l,m}*
                \int RL(k)*RL'(k)*jl(k*R) k^4 dk^3,
        ****************************************************/

ftrace_region_begin("OLP_Kin_4.5");
        for(int m=-l; m<=l; m++){ 

          ComplexSH(l,m,theta,phi,SH,dSHt,dSHp);
          SH[1]   = -SH[1];
          dSHt[1] = -dSHt[1];
          dSHp[1] = -dSHp[1];

          CY  = Complex_inline(SH[0],SH[1]);
          CYt = Complex_inline(dSHt[0],dSHt[1]);
          CYp = Complex_inline(dSHp[0],dSHp[1]);


          for (int LMul0=0; LMul0<ai_max_LMul0; LMul0++){
            const int L0 = ai_L0[LMul0];
            const int num0 = ai_offset0[LMul0];
#pragma _NEC vector
#pragma _NEC ivdep
            for (int LMul1=0; LMul1<ai_max_LMul1; LMul1++){
              const int L1 = ai_L1[LMul1];
              const int num1 = ai_offset1[LMul1];

              const int Ls = -L0 + L1 + l;  

              if (abs(L1-l)<=L0 && L0<=(L1+l) ){

                const dcomplex Cpow = Im_pow(-1,Ls);
                const dcomplex CY1  = Cmul_inline(Cpow,CY);
                const dcomplex CYt1 = Cmul_inline(Cpow,CYt);
                const dcomplex CYp1 = Cmul_inline(Cpow,CYp);

/*
#ifdef _NEC
#pragma _NEC loop_count(7)
#pragma _NEC unroll(7)
                for (int L0M0=0; L0M0<=6; L0M0++){//for (M0=-L0; M0<=L0; M0++){
                  if(L0M0<=2*L0){
#else
*/
                for (int L0M0=0; L0M0<=2*L0; L0M0++){//for (M0=-L0; M0<=L0; M0++){
                  
/*#endif*/
                  const int M0 = L0M0 - L0;

                  const int M1 = M0 - m;
                  const int L1M1 = L1 + M1;

                  
                  if (abs(M1)<=L1){
#ifdef _NEC_OLP_MG
                    gant = GetGaunt(L0,M0,L1,M1,l,m,Gaunt_list, ai_L1_LIMIT_SQ, ai_LL_LIMIT_SQ);
#else
                    gant = Gaunt(L0,M0,L1,M1,l,m);
#endif
                    /* S */ 

                    tmp0 = gant*SumS0[LMul0][LMul1];
                    Ctmp2 = CRmul(CY1,tmp0);
                    TmpOLP[num0+L0M0][num1+L1M1] =
                      Cadd_inline(TmpOLP[num0+L0M0][num1+L1M1],Ctmp2);

                    /* dS/dr */ 

                    tmp0 = gant*SumSr0[LMul0][LMul1];
                    Ctmp2 = CRmul(CY1,tmp0);
                    TmpOLPr[num0+L0M0][num1+L1M1] =
                      Cadd_inline(TmpOLPr[num0+L0M0][num1+L1M1],Ctmp2);

                    /* dS/dt */ 

                    tmp0 = gant*SumS0[LMul0][LMul1];
                    Ctmp2 = CRmul(CYt1,tmp0);
                    TmpOLPt[num0+L0M0][num1+L1M1] =
                      Cadd_inline(TmpOLPt[num0+L0M0][num1+L1M1],Ctmp2);

                    /* dS/dp */ 

                    tmp0 = gant*SumS0[LMul0][LMul1];
                    Ctmp2 = CRmul(CYp1,tmp0);
                    TmpOLPp[num0+L0M0][num1+L1M1] =
                      Cadd_inline(TmpOLPp[num0+L0M0][num1+L1M1],Ctmp2);

                    /* K */ 

                    tmp0 = gant*SumK0[LMul0][LMul1];
                    Ctmp2 = CRmul(CY1,tmp0);
                    TmpKin[num0+L0M0][num1+L1M1] =
                      Cadd_inline(TmpKin[num0+L0M0][num1+L1M1],Ctmp2);

                    /* dK/dr */ 

                    tmp0 = gant*SumKr0[LMul0][LMul1];
                    Ctmp2 = CRmul(CY1,tmp0);
                    TmpKinr[num0+L0M0][num1+L1M1] =
                      Cadd_inline(TmpKinr[num0+L0M0][num1+L1M1],Ctmp2);

                    /* dK/dt */ 

                    tmp0 = gant*SumK0[LMul0][LMul1];
                    Ctmp2 = CRmul(CYt1,tmp0);
                    TmpKint[num0+L0M0][num1+L1M1] =
                      Cadd_inline(TmpKint[num0+L0M0][num1+L1M1],Ctmp2);

                    /* dK/dp */ 

                    tmp0 = gant*SumK0[LMul0][LMul1];
                    Ctmp2 = CRmul(CYp1,tmp0);
                    TmpKinp[num0+L0M0][num1+L1M1] =
                      Cadd_inline(TmpKinp[num0+L0M0][num1+L1M1],Ctmp2);

		              }
/*#ifdef _NEC                  
                  }
#endif                  */
		            }
		          }
	  	      }
          }
        } /* m */
        
ftrace_region_end("OLP_Kin_4.5");
      } /* l */

      /* free SphB and SphBp */
/*
      for(l=0; l<(Lmax_Four_Int+3); l++){ 
        free(SphB[l]);
      }
      free(SphB);

      for(l=0; l<(Lmax_Four_Int+3); l++){ 
        free(SphBp[l]);
      }
      free(SphBp);
*/
      
      //printf("free(SphBp)\n"); fflush(stdout);        

      /****************************************************
                         Complex to Real
      ****************************************************/
#pragma _NEC novector
          for(int vi0 = 0; vi0 < ai_len_vec_LM0; ++vi0){
            const int L0 = ai_vec_L0[vi0];
            const int num0 = ai_vec_num0[vi0];
            const int L0M0 = vi0 - num0;

           // const dcomplex* c2r=Comp2Real[L0][L0M0];

            for(int vi1 = 0; vi1 < ai_len_vec_LM1; ++vi1){

                  CmatS0[vi1] = Complex_inline(0.0,0.0);
                  CmatSr[vi1] = Complex_inline(0.0,0.0);
                  CmatSt[vi1] = Complex_inline(0.0,0.0);
                  CmatSp[vi1] = Complex_inline(0.0,0.0);

                  CmatK0[vi1] = Complex_inline(0.0,0.0);
                  CmatKr[vi1] = Complex_inline(0.0,0.0);
                  CmatKt[vi1] = Complex_inline(0.0,0.0);
                  CmatKp[vi1] = Complex_inline(0.0,0.0);
            }

#pragma _NEC novector
#pragma _NEC unroll(7)
                  for (int mk=0; mk<=2*L0; mk++){

#pragma _NEC nosync            
#pragma _NEC vector
#pragma _NEC ivdep
            for(int vi1 = 0; vi1 < ai_len_vec_LM1; ++vi1){
            
              //const int L1 = ai_vec_L1[vi1];
              //const int num1 = ai_vec_num1[vi1];
              //const int L1M1 = vi1 - num1;
              
/*
                  dcomplex CsumS0 = Complex_inline(0.0,0.0);
                  dcomplex CsumSr = Complex_inline(0.0,0.0);
                  dcomplex CsumSt = Complex_inline(0.0,0.0);
                  dcomplex CsumSp = Complex_inline(0.0,0.0);

                  dcomplex CsumK0 = Complex_inline(0.0,0.0);
                  dcomplex CsumKr = Complex_inline(0.0,0.0);
                  dcomplex CsumKt = Complex_inline(0.0,0.0);
                  dcomplex CsumKp = Complex_inline(0.0,0.0);
*/
     	            //for (k=-L0; k<=L0; k++){

                  
                  

                    const dcomplex Ctmp1 = Conjg_inline(Comp2Real[L0][L0M0][mk]);

                    /* S */

                    dcomplex Ctmp0 = TmpOLP[num0+mk][vi1];
                    dcomplex Ctmp2 = Cmul_inline(Ctmp1,Ctmp0);
                    const dcomplex CsumS0 = Cadd_inline(CmatS0[vi1],Ctmp2);
 
                    /* dS/dr */

                    Ctmp0 = TmpOLPr[num0+mk][vi1];
                    Ctmp2 = Cmul_inline(Ctmp1,Ctmp0);
                    const dcomplex CsumSr = Cadd_inline(CmatSr[vi1],Ctmp2);

                    /* dS/dt */

                    Ctmp0 = TmpOLPt[num0+mk][vi1];
                    Ctmp2 = Cmul_inline(Ctmp1,Ctmp0);
                    const dcomplex CsumSt = Cadd_inline(CmatSt[vi1],Ctmp2);

                    /* dS/dp */

                    Ctmp0 = TmpOLPp[num0+mk][vi1];
                    Ctmp2 = Cmul_inline(Ctmp1,Ctmp0);
                    const dcomplex CsumSp = Cadd_inline(CmatSp[vi1],Ctmp2);

                    /* K */

                    Ctmp0 = TmpKin[num0+mk][vi1];
                    Ctmp2 = Cmul_inline(Ctmp1,Ctmp0);
                    const dcomplex CsumK0 = Cadd_inline(CmatK0[vi1],Ctmp2);

                    /* dK/dr */

                    Ctmp0 = TmpKinr[num0+mk][vi1];
                    Ctmp2 = Cmul_inline(Ctmp1,Ctmp0);
                    const dcomplex CsumKr = Cadd_inline(CmatKr[vi1],Ctmp2);

                    /* dK/dt */

                    Ctmp0 = TmpKint[num0+mk][vi1];
                    Ctmp2 = Cmul_inline(Ctmp1,Ctmp0);
                    const dcomplex CsumKt = Cadd_inline(CmatKt[vi1],Ctmp2);

                    /* dK/dp */

                    Ctmp0 = TmpKinp[num0+mk][vi1];
                    Ctmp2 = Cmul_inline(Ctmp1,Ctmp0);
                    const dcomplex CsumKp = Cadd_inline(CmatKp[vi1],Ctmp2);

		              

                  CmatS0[vi1] = CsumS0;
                  CmatSr[vi1] = CsumSr;
                  CmatSt[vi1] = CsumSt;
                  CmatSp[vi1] = CsumSp;

                  CmatK0[vi1] = CsumK0;
                  CmatKr[vi1] = CsumKr;
                  CmatKt[vi1] = CsumKt;
                  CmatKp[vi1] = CsumKp;
                  
            }
	        }

#pragma _NEC vector
#pragma _NEC ivdep
            for(int vi1 = 0; vi1 < ai_len_vec_LM1; ++vi1){
            
              const int L1 = ai_vec_L1[vi1];
              //const int LMul1 = ai_vec_LMul1[vi1];
              const int num1 = ai_vec_num1[vi1];
              const int L1M1 = vi1 - num1;
              {
                
                  dcomplex CsumS_Lx = Complex_inline(0.0,0.0);
                  dcomplex CsumS_Ly = Complex_inline(0.0,0.0);
                  dcomplex CsumS_Lz = Complex_inline(0.0,0.0);

                  dcomplex CsumS0 = Complex_inline(0.0,0.0);
                  dcomplex CsumSr = Complex_inline(0.0,0.0);
                  dcomplex CsumSt = Complex_inline(0.0,0.0);
                  dcomplex CsumSp = Complex_inline(0.0,0.0);
                  dcomplex CsumK0 = Complex_inline(0.0,0.0);
                  dcomplex CsumKr = Complex_inline(0.0,0.0);
                  dcomplex CsumKt = Complex_inline(0.0,0.0);
                  dcomplex CsumKp = Complex_inline(0.0,0.0);

#pragma _NEC novector
#pragma _NEC noinner
#pragma _NEC unroll(7)
     	            for (int mk=0; mk<=2*L1; mk++){//for (k=-L1; k<=L1; k++){

                    /*** S_Lx ***/ 

                    /*  Y k+1 */
                    if (mk-L1<L1){//if (k<L1){
                      coe0 = sqrt((double)((2*L1-mk)*(mk+1)));
                      Ctmp1 = Cmul_inline(CmatS0[num1+mk+1],Comp2Real[L1][L1M1][mk]);
                      Ctmp1.r = 0.5*coe0*Ctmp1.r;
                      Ctmp1.i = 0.5*coe0*Ctmp1.i;
                      CsumS_Lx = Cadd_inline(CsumS_Lx,Ctmp1);
		                }

                    /*  Y k-1 */
                    if (-L1<mk-L1){//if (-L1<k){
                      coe1 = sqrt((double)((mk)*(2*L1-mk+1)));
                      Ctmp1 = Cmul_inline(CmatS0[num1+mk-1],Comp2Real[L1][L1M1][mk]);
                      Ctmp1.r = 0.5*coe1*Ctmp1.r;
                      Ctmp1.i = 0.5*coe1*Ctmp1.i;
                      CsumS_Lx = Cadd_inline(CsumS_Lx,Ctmp1);
		                }

                    /*** S_Ly ***/ 

                    /*  Y k+1 */

                    if (mk-L1<L1){//if (k<L1){
                      Ctmp1 = Cmul_inline(CmatS0[num1+mk+1],Comp2Real[L1][L1M1][mk]);
                      Ctmp2.r = 0.5*coe0*Ctmp1.i;
                      Ctmp2.i =-0.5*coe0*Ctmp1.r;
                      CsumS_Ly = Cadd_inline(CsumS_Ly,Ctmp2);
		                }

                    /*  Y k-1 */

                    if (-L1<mk-L1){//if (-L1<k){
                      Ctmp1 = Cmul_inline(CmatS0[num1+mk-1],Comp2Real[L1][L1M1][mk]);
                      Ctmp2.r =-0.5*coe1*Ctmp1.i;
                      Ctmp2.i = 0.5*coe1*Ctmp1.r;
                      CsumS_Ly = Cadd_inline(CsumS_Ly,Ctmp2);
                    }

                    /*** S_Lz ***/ 

                    Ctmp1    = Cmul_inline(CmatS0[num1+mk],Comp2Real[L1][L1M1][mk]);
                    Ctmp1.r = (double)(mk-L1)*Ctmp1.r;;
                    Ctmp1.i = (double)(mk-L1)*Ctmp1.i;
                    CsumS_Lz = Cadd_inline(CsumS_Lz,Ctmp1);

                    /* S */ 

                    Ctmp1 = Cmul_inline(CmatS0[num1+mk],Comp2Real[L1][L1M1][mk]);
                    CsumS0 = Cadd_inline(CsumS0,Ctmp1);

                    /* dS/dr */ 

                    Ctmp1 = Cmul_inline(CmatSr[num1+mk],Comp2Real[L1][L1M1][mk]);
                    CsumSr = Cadd_inline(CsumSr,Ctmp1);

                    /* dS/dt */ 

                    Ctmp1 = Cmul_inline(CmatSt[num1+mk],Comp2Real[L1][L1M1][mk]);
                    CsumSt = Cadd_inline(CsumSt,Ctmp1);

                    /* dS/dp */ 

                    Ctmp1 = Cmul_inline(CmatSp[num1+mk],Comp2Real[L1][L1M1][mk]);
                    CsumSp = Cadd_inline(CsumSp,Ctmp1);

                    /* K */ 

                    Ctmp1 = Cmul_inline(CmatK0[num1+mk],Comp2Real[L1][L1M1][mk]);
                    CsumK0 = Cadd_inline(CsumK0,Ctmp1);

                    /* dK/dr */ 

                    Ctmp1 = Cmul_inline(CmatKr[num1+mk],Comp2Real[L1][L1M1][mk]);
                    CsumKr = Cadd_inline(CsumKr,Ctmp1);

                    /* dK/dt */ 

                    Ctmp1 = Cmul_inline(CmatKt[num1+mk],Comp2Real[L1][L1M1][mk]);
                    CsumKt = Cadd_inline(CsumKt,Ctmp1);

                    /* dK/dp */ 

                    Ctmp1 = Cmul_inline(CmatKp[num1+mk],Comp2Real[L1][L1M1][mk]);
                    CsumKp = Cadd_inline(CsumKp,Ctmp1);

		              }
                  
                  OLP_L[0][Mc_AN][h_AN][vi0][vi1] = 8.0*CsumS_Lx.i; 
                  OLP_L[1][Mc_AN][h_AN][vi0][vi1] = 8.0*CsumS_Ly.i; 
                  OLP_L[2][Mc_AN][h_AN][vi0][vi1] = 8.0*CsumS_Lz.i;
                  
                  /* add a small value for stabilization of eigenvalue routine */

                  /*OLP[0][Mc_AN][h_AN][vi0][vi1] = 8.0*CsumS0.r + 1.0*rnd(1.0e-13);*/
                  OLP[0][Mc_AN][h_AN][vi0][vi1] = 8.0*CsumS0.r + 1.0*(1.0e-13);
                  H0[0][Mc_AN][h_AN][vi0][vi1] = 4.0*CsumK0.r;

                  double olp1, olp2, olp3, h01,h02, h03;
                  if (h_AN!=0){

                    if (fabs(siT)<10e-14){

                      olp1 =
                  -8.0*(siT*coP*CsumSr.r + coT*coP/r*CsumSt.r);

                      olp2 =
                  -8.0*(siT*siP*CsumSr.r + coT*siP/r*CsumSt.r);

                      olp3 =
                  -8.0*(coT*CsumSr.r - siT/r*CsumSt.r);

                      h01 =
                  -4.0*(siT*coP*CsumKr.r + coT*coP/r*CsumKt.r);

                      h02 =
                  -4.0*(siT*siP*CsumKr.r + coT*siP/r*CsumKt.r);

                      h03 =
                  -4.0*(coT*CsumKr.r - siT/r*CsumKt.r);
                    }

                    else{

                      olp1 =
                  -8.0*(siT*coP*CsumSr.r + coT*coP/r*CsumSt.r
                      - siP/siT/r*CsumSp.r);

                      olp2 =
                  -8.0*(siT*siP*CsumSr.r + coT*siP/r*CsumSt.r
                      + coP/siT/r*CsumSp.r);

                      olp3 =
                  -8.0*(coT*CsumSr.r - siT/r*CsumSt.r);

                      h01 =
                  -4.0*(siT*coP*CsumKr.r + coT*coP/r*CsumKt.r
                      - siP/siT/r*CsumKp.r);

                      h02 =
                  -4.0*(siT*siP*CsumKr.r + coT*siP/r*CsumKt.r
                      + coP/siT/r*CsumKp.r);

                      h03 =
                  -4.0*(coT*CsumKr.r - siT/r*CsumKt.r);
                    }

                  }
                  else{
                    olp1 = 0.0;
                    olp2 = 0.0;
                    olp3 = 0.0;
                    h01 = 0.0;
                    h02 = 0.0;
                    h03 = 0.0;
                  }

                 
                      OLP[1][Mc_AN][h_AN][vi0][vi1] = olp1;

                      OLP[2][Mc_AN][h_AN][vi0][vi1] = olp2;

                      OLP[3][Mc_AN][h_AN][vi0][vi1] = olp3;
                      H0[1][Mc_AN][h_AN][vi0][vi1] = h01;

                      H0[2][Mc_AN][h_AN][vi0][vi1] =h02;

                      H0[3][Mc_AN][h_AN][vi0][vi1] =h03;
                }

	          }
      }/* vi0 */

      
      free(ai_L0);
      free(ai_L1);
      free(ai_offset0);
      free(ai_offset1);
      free(ai_vec_L0);
      free(ai_vec_L1);
      free(ai_vec_num0);
      free(ai_vec_num1);
      //free(ai_vec_LMul1);

      dtime(&Etime_atom);
      time_per_atom[Gc_AN] += Etime_atom - Stime_atom;
    } /* end of loop for Nloop */

    free(normk_indexes);
    free(tmp_SphB_D);
    free(tmp_SphBp_D);


    /* freeing of arrays */
    Free2D_dcomplex(TmpOLP);
    Free2D_dcomplex(TmpOLPr);
    Free2D_dcomplex(TmpOLPt);
    Free2D_dcomplex(TmpOLPp);
    Free2D_dcomplex(TmpKin);
    Free2D_dcomplex(TmpKinr);
    Free2D_dcomplex(TmpKint);
    Free2D_dcomplex(TmpKinp);
	
    Free2D_double(SumS0);
    Free2D_double(SumK0);
    Free2D_double(SumSr0);
    Free2D_double(SumKr0);
	
    free(CmatS0);
    free(CmatSr);
    free(CmatSt);
    free(CmatSp);
    free(CmatK0);
    free(CmatKr);
    free(CmatKt);
    free(CmatKp);
  /*
    Free2D_dcomplex(CmatS0);
    Free2D_dcomplex(CmatSr);
    Free2D_dcomplex(CmatSt);
    Free2D_dcomplex(CmatSp);
    Free2D_dcomplex(CmatK0);
    Free2D_dcomplex(CmatKr);
    Free2D_dcomplex(CmatKt);
    Free2D_dcomplex(CmatKp);
*/
    //AITUNE//
    //not used//Free2D_double(ai_Gaunt);
	
  } /* #pragma omp parallel */
  
  /****************************************************
                   freeing of arrays:
  ****************************************************/

  free(OneD2h_AN);
  free(OneD2Mc_AN);
  free(Gaunt_list);

  /* for time */
  dtime(&TEtime);
  time0 = TEtime - TStime;

  return time0;
}


#ifdef kcomp
dcomplex****** Allocate6D_dcomplex(int size_1, int size_2, int size_3, int size_4, int size_5, int size_6)
#else
static inline dcomplex****** Allocate6D_dcomplex(int size_1, int size_2, int size_3, int size_4, int size_5, int size_6)
#endif
{ 
  int i, j, k, l, m, p;

  dcomplex****** buffer = (dcomplex******)malloc(sizeof(dcomplex*****)*size_1);
  buffer[0] = (dcomplex*****)malloc(sizeof(dcomplex****)*size_1*size_2);
  buffer[0][0] = (dcomplex****)malloc(sizeof(dcomplex***)*size_1*size_2*size_3);
  buffer[0][0][0] = (dcomplex***)malloc(sizeof(dcomplex**)*size_1*size_2*size_3*size_4);
  buffer[0][0][0][0] = (dcomplex**)malloc(sizeof(dcomplex*)*size_1*size_2*size_3*size_4*size_5);
  buffer[0][0][0][0][0] = (dcomplex*)malloc(sizeof(dcomplex)*size_1*size_2*size_3*size_4*size_5*size_6);
		
  for (i=0; i<size_1; i++){
    buffer[i] = buffer[0] + i * size_2;
    for (j=0; j<size_2; j++){
      buffer[i][j] = buffer[0][0] + (i * size_2 + j) * size_3;
      for (k=0; k<size_3; k++){
	buffer[i][j][k] = buffer[0][0][0] + ((i * size_2 + j) * size_3 + k) * size_4;
	for (l=0; l<size_4; l++){
	  buffer[i][j][k][l] = buffer[0][0][0][0] + (((i * size_2 + j) * size_3 + k) * size_4 + l) * size_5;
	  for (m=0; m<size_5; m++){
	    buffer[i][j][k][l][m] = buffer[0][0][0][0][0] + ((((i * size_2 + j) * size_3 + k) * size_4 + l) * size_5 + m) * size_6;
	    for (p=0; p<size_6; p++) buffer[i][j][k][l][m][p] = Complex_inline(0.0,0.0);
	  }
	}
      }
    }
  }

  return buffer;
}


#ifdef kcomp
dcomplex**** Allocate4D_dcomplex(int size_1, int size_2, int size_3, int size_4)
#else
static inline dcomplex**** Allocate4D_dcomplex(int size_1, int size_2, int size_3, int size_4)
#endif
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
	  buffer[i][j][k][l] = Complex_inline(0.0,0.0);
	}
      }
    }
  }

  return buffer;
}
#ifdef kcomp
double**** Allocate4D_double(int size_1, int size_2, int size_3, int size_4)
#else
static inline double**** Allocate4D_double(int size_1, int size_2, int size_3, int size_4)
#endif
{ 
  int i, j, k, l;

  double**** buffer = (double****)malloc(sizeof(double***)*size_1);
  buffer[0] = (double***)malloc(sizeof(double**)*size_1*size_2);
  buffer[0][0] = (double**)malloc(sizeof(double*)*size_1*size_2*size_3);
  buffer[0][0][0] = (double*)malloc(sizeof(double)*size_1*size_2*size_3*size_4);
		
  for (i=0; i<size_1; i++){
    buffer[i] = buffer[0] + i * size_2;
    for (j=0; j<size_2; j++){
      buffer[i][j] = buffer[0][0] + (i * size_2 + j) * size_3;
      for (k=0; k<size_3; k++){
	buffer[i][j][k] = buffer[0][0][0] + ((i * size_2 + j) * size_3 + k) * size_4;
	for (l=0; l<size_4; l++){
	  buffer[i][j][k][l] = 0.0;
	}
      }
    }
  }

  return buffer;
}

#ifdef kcomp
dcomplex** Allocate2D_dcomplex(int size_1, int size_2)
#else
static inline dcomplex** Allocate2D_dcomplex(int size_1, int size_2)
#endif
{ 
  int i, j;

  dcomplex** buffer = (dcomplex**)malloc(sizeof(dcomplex*)*size_1);
  buffer[0] = (dcomplex*)malloc(sizeof(dcomplex)*size_1*size_2);

  for (i=0; i<size_1; i++){
    buffer[i] = buffer[0] + i * size_2;
    for (j=0; j<size_2; j++){
      buffer[i][j] = Complex_inline(0.0,0.0);
    }
  }

  return buffer;
}

#ifdef kcomp
double** Allocate2D_double(int size_1, int size_2)
#else
static inline double** Allocate2D_double(int size_1, int size_2)
#endif
{ 
  int i, j;

  double** buffer = (double**)malloc(sizeof(double*)*size_1);
  buffer[0] = (double*)malloc(sizeof(double)*size_1*size_2);

  for (i=0; i<size_1; i++){
    buffer[i] = buffer[0] + i * size_2;
    for (j=0; j<size_2; j++){
      buffer[i][j] = 0.0;
    }
  }

  return buffer;
}

void Free6D_dcomplex(dcomplex****** buffer)
{ 
  free(buffer[0][0][0][0][0]);
  free(buffer[0][0][0][0]);
  free(buffer[0][0][0]);
  free(buffer[0][0]);
  free(buffer[0]);
  free(buffer);
}

void Free4D_dcomplex(dcomplex**** buffer)
{ 
  free(buffer[0][0][0]);
  free(buffer[0][0]);
  free(buffer[0]);
  free(buffer);
}

void Free4D_double(double**** buffer)
{ 
  free(buffer[0][0][0]);
  free(buffer[0][0]);
  free(buffer[0]);
  free(buffer);
}

void Free2D_dcomplex(dcomplex** buffer)
{ 
  free(buffer[0]);
  free(buffer);
}


void Free2D_double(double** buffer)
{ 
  free(buffer[0]);
  free(buffer);
}


