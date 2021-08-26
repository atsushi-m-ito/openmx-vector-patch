/**********************************************************************
  Set_Orbitals_Grid.c:

   Set_Orbitals_Grid.c is a subroutine to calculate the value of basis
   functions on each grid point.

  Log of Set_Orbitals_Grid.c:

     22/Nov/2001  Released by T.Ozaki

***********************************************************************/
#include "mpi.h"
#include <omp.h>

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <math.h>
#include "openmx_common.h"



#define AITURN_20200914
#define AITURN_20200914b
#ifdef _NEC 
#define VEC_NEC
#endif
#define VEC_NEC

inline static 
void Get_Grid_XYZ_inline(int GN, double xyz[4])
{
  int n1,n2,n3;

  n1 = GN/(Ngrid2*Ngrid3);
  n2 = (GN - n1*(Ngrid2*Ngrid3))/Ngrid3;
  n3 = GN - n1*(Ngrid2*Ngrid3) - n2*Ngrid3;

  xyz[1] = (double)n1*gtv[1][1] + (double)n2*gtv[2][1]
         + (double)n3*gtv[3][1] + Grid_Origin[1];
  xyz[2] = (double)n1*gtv[1][2] + (double)n2*gtv[2][2]
         + (double)n3*gtv[3][2] + Grid_Origin[2];
  xyz[3] = (double)n1*gtv[1][3] + (double)n2*gtv[2][3]
         + (double)n3*gtv[3][3] + Grid_Origin[3];
}


double Set_Orbitals_Grid(int Cnt_kind)
{
  double time0;  
  double TStime,TEtime;
  int numprocs,myid,tag=999,ID,IDS,IDR;
  double Stime_atom,Etime_atom;

  MPI_Status stat;
  MPI_Request request;

  /* for OpenMP */
  int OMPID,Nthrds,Nprocs;

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);
  
  dtime(&TStime);

  /*****************************************************
                Calculate orbitals on grids
  *****************************************************/

  for (int Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){

    dtime(&Stime_atom);

    const int Gc_AN = M2G[Mc_AN];    
    const int Cwan = WhatSpecies[Gc_AN];

	int NO0;
    if (Cnt_kind==0)  NO0 = Spe_Total_NO[Cwan];
    else              NO0 = Spe_Total_CNO[Cwan]; 

#pragma omp parallel shared(Comp2Real,Spe_PAO_RWF,Spe_Num_Basis,Spe_MaxL_Basis,Spe_PAO_RV,Spe_Num_Mesh_PAO,List_YOUSO,Orbs_Grid,Cnt_kind,Gxyz,atv,CellListAtom,GridListAtom,GridN_Atom,Gc_AN,Cwan,Mc_AN,NO0) private(OMPID,Nthrds,Nprocs)
    {
      

      /* Radial */
      
      /* allocation of array */

	if (Cnt_kind==0){
		
      const int wan = Cwan;

	  int ai_max_LMul0 = 0;        
      for (int L0=0; L0<=Spe_MaxL_Basis[wan]; L0++){
        ai_max_LMul0 += Spe_Num_Basis[wan][L0];		  
      }

#ifdef VEC_NEC
      if((ai_max_LMul0 <= 15) && (Spe_MaxL_Basis[wan] < 4)){

#ifdef AITURN_20200914

/* AI_TURN */
      int* ai_rf_L0 = (int*)malloc(sizeof(int)*ai_max_LMul0);
	  int* ai_rf_Mul0 = (int*)malloc(sizeof(int)*ai_max_LMul0);
		
      {
	    int n0 = 0;
        for (int L0=0; L0<=Spe_MaxL_Basis[wan]; L0++){
	        for (int Mul0=0; Mul0<Spe_Num_Basis[wan][L0]; Mul0++){
				ai_rf_L0[n0] = L0;
				ai_rf_Mul0[n0] = Mul0;
            n0++;
            
          }
        }
        
	  }
#endif

#ifdef AITURN_20200914b
		int rank_for_m = 1;
		{
			int bi = Spe_Num_Mesh_PAO[wan];
			while (bi > 0) {
				rank_for_m ++;
				bi >>= 1;
			}
		}

#endif
      
      /* get info. on OpenMP */ 

      const int OMPID = omp_get_thread_num();
      const int Nthrds = omp_get_num_threads();
      const int Nprocs = omp_get_num_procs();

	  const int size_n = GridN_Atom[Gc_AN];
	  double* RF = (double*)malloc(sizeof(double) * size_n * ai_max_LMul0);
	  double* AF = (double*)malloc(sizeof(double) * size_n * 16);
      
	const double upper_limit_R = Spe_PAO_RV[wan][Spe_Num_Mesh_PAO[wan]-1];

#pragma _NEC vector
#pragma _NEC ivdep
      for (int Nc=OMPID*GridN_Atom[Gc_AN]/Nthrds; Nc<(OMPID+1)*GridN_Atom[Gc_AN]/Nthrds; Nc++){

	const int GNc = GridListAtom[Mc_AN][Nc]; 
	const int GRc = CellListAtom[Mc_AN][Nc];

	double Cxyz_1,Cxyz_2,Cxyz_3;
	//Get_Grid_XYZ_inline(GNc,Cxyz);
	{
		int n1,n2,n3;

		n1 = GNc/(Ngrid2*Ngrid3);
		n2 = (GNc - n1*(Ngrid2*Ngrid3))/Ngrid3;
		n3 = GNc - n1*(Ngrid2*Ngrid3) - n2*Ngrid3;

		Cxyz_1 = (double)n1*gtv[1][1] + (double)n2*gtv[2][1]
				+ (double)n3*gtv[3][1] + Grid_Origin[1];
		Cxyz_2 = (double)n1*gtv[1][2] + (double)n2*gtv[2][2]
				+ (double)n3*gtv[3][2] + Grid_Origin[2];
		Cxyz_3 = (double)n1*gtv[1][3] + (double)n2*gtv[2][3]
				+ (double)n3*gtv[3][3] + Grid_Origin[3];
	}
	double x = Cxyz_1 + atv[GRc][1] - Gxyz[Gc_AN][1]; 
	double y = Cxyz_2 + atv[GRc][2] - Gxyz[Gc_AN][2]; 
	double z = Cxyz_3 + atv[GRc][3] - Gxyz[Gc_AN][3];



          /* Get_Orbitals(Cwan,x,y,z,Chi0); */
          /* start of inlining of Get_Orbitals */


          /* xyz2spherical(x,y,z,0.0,0.0,0.0,S_coordinate); */
          /* start of inlining of xyz2spherical */

	  const double Min_r = 10e-15;
	  double dum = x*x + y*y; 
	  double r = sqrt(dum + z*z);
	  double r1 = sqrt(dum);

	  const double dum1t = (r<fabs(z)) ? ((z>0.0)?1.0:-1.0) : z/r;
	  const double theta = (Min_r<=r) ? acos(dum1t) : 0.5*PI;

	  const double dum1p = (r1<fabs(y)) ? ((y>0.0)?1.0:-1.0) : y/r1;  
	  const double phi = (Min_r<=r1) ? ((0.0<=x) ? asin(dum1p) : PI - asin(dum1p)) : 0.0;



	  double R = r;
	  double Q = theta;
	  double P = phi;
/*
			if(out_og_max < P)out_og_max = P;
			if(out_og_min > P)out_og_min = P;
*/
	  /* end of inlining of xyz2spherical */

	  int po = 0;
	  int mp_min = 0;
	  int mp_max = Spe_Num_Mesh_PAO[wan] - 1;
	  int m;

	  if (Spe_PAO_RV[wan][Spe_Num_Mesh_PAO[wan]-1]<R){

		#pragma _NEC unroll(15) 
		for(int i1 = 0; i1 < 15; ++i1){
			if(i1 < ai_max_LMul0){
				//Orbs_Grid[Mc_AN][Nc][i1] = 0.0;
				RF[i1 * size_n + Nc] = 0.0;
			}
		}
	  
	  }
	  else
	  {

		  
	    /* Angular */
	    const double siQ = sin(Q);
	    const double coQ = cos(Q);
	    const double siP = sin(P);
	    const double coP = cos(P);

		{
			AF[(0+0)*size_n + Nc] = 0.282094791773878;
	    
			const double dum = 0.48860251190292*siQ;
			AF[(1+0)*size_n + Nc] = dum*coP;
			AF[(1+1)*size_n + Nc] = dum*siP;
			AF[(1+2)*size_n + Nc] = 0.48860251190292*coQ;
	      
			const double dum1 = siQ*siQ;
			const double dum2 = 1.09254843059208*siQ*coQ;
			AF[(4+0)*size_n + Nc] = 0.94617469575756*coQ*coQ - 0.31539156525252;
			AF[(4+1)*size_n + Nc] = 0.54627421529604*dum1*(1.0 - 2.0*siP*siP);
			AF[(4+2)*size_n + Nc] = 1.09254843059208*dum1*siP*coP;
			AF[(4+3)*size_n + Nc] = dum2*coP;
			AF[(4+4)*size_n + Nc] = dum2*siP;
	      
	      	AF[(9+0)*size_n + Nc] = 0.373176332590116*(5.0*coQ*coQ*coQ - 3.0*coQ);
			AF[(9+1)*size_n + Nc] = 0.457045799464466*coP*siQ*(5.0*coQ*coQ - 1.0);
			AF[(9+2)*size_n + Nc] = 0.457045799464466*siP*siQ*(5.0*coQ*coQ - 1.0);
			AF[(9+3)*size_n + Nc] = 1.44530572132028*siQ*siQ*coQ*(coP*coP-siP*siP);
			AF[(9+4)*size_n + Nc] = 2.89061144264055*siQ*siQ*coQ*siP*coP;
			AF[(9+5)*size_n + Nc] = 0.590043589926644*siQ*siQ*siQ*(4.0*coP*coP*coP - 3.0*coP);
			AF[(9+6)*size_n + Nc] = 0.590043589926644*siQ*siQ*siQ*(3.0*siP - 4.0*siP*siP*siP);
	    }
		 

		double rm;
		if (R<Spe_PAO_RV[wan][0]){
			m = 4;
			rm = Spe_PAO_RV[wan][m];
		}else{				
			
#ifdef AITURN_20200914b
			
			#pragma _NEC loop_count(16)
			#pragma _NEC unroll(16)
			//#pragma _NEC novector
			for(int rank = 0; rank < 16; ++rank){
				if(rank < rank_for_m){
				if((mp_max - mp_min) > 1){
					m = (mp_min + mp_max) / 2;
					if (Spe_PAO_RV[wan][m] < R) mp_min = m;
					else mp_max = m;
				}
				}
			}
			m = mp_max;
			
#else
			do{
				m = (mp_min + mp_max)/2;
				if (Spe_PAO_RV[wan][m]<R)
					mp_min = m;
				else 
					mp_max = m;
			}
			while((mp_max-mp_min)!=1);
			m = mp_max;
#endif
			rm = R;
		}

	    double h1 = Spe_PAO_RV[wan][m-1] - Spe_PAO_RV[wan][m-2];
	    double h2 = Spe_PAO_RV[wan][m]   - Spe_PAO_RV[wan][m-1];
	    double h3 = Spe_PAO_RV[wan][m+1] - Spe_PAO_RV[wan][m];
		if (m==1){
			h1 = -(h2+h3);
		}
		else if (m==(Spe_Num_Mesh_PAO[wan]-1)){
			h3 = -(h1+h2);
		}

	    double x1 = rm - Spe_PAO_RV[wan][m-1];
	    double x2 = rm - Spe_PAO_RV[wan][m];
	    double y1 = x1/h2;
	    double y2 = x2/h2;
	    double y12 = y1*y1;
	    double y22 = y2*y2;
		
	    double dum = h1 + h2;
	    double dum1 = h1/h2/dum;
	    double dum2 = h2/h1/dum;
	    dum = h2 + h3;
	    double dum3 = h2/h3/dum;
	    double dum4 = h3/h2/dum;


		#pragma _NEC unroll(15) 
		for(int n0 = 0; n0 < 15; ++n0){
		  if(n0 < ai_max_LMul0){
	  	    const int L0=ai_rf_L0[n0];
		    const int Mul0=ai_rf_Mul0[n0];

			double f1 = (m!=1) ? Spe_PAO_RWF[wan][L0][Mul0][m-2] : Spe_PAO_RWF[wan][L0][Mul0][m+1];
			double f2 = Spe_PAO_RWF[wan][L0][Mul0][m-1];
			double f3 = Spe_PAO_RWF[wan][L0][Mul0][m];
			double f4 = (m!=(Spe_Num_Mesh_PAO[wan]-1)) ? Spe_PAO_RWF[wan][L0][Mul0][m+1] : Spe_PAO_RWF[wan][L0][Mul0][m-2];


			double dum = f3 - f2;
			double g1 = dum*dum1 + (f2-f1)*dum2;
			double g2 = (f4-f3)*dum3 + dum*dum4;

			double f =  y22*(3.0*f2 + h2*g1 + (2.0*f2 + h2*g1)*y2)
			+ y12*(3.0*f3 - h2*g2 - (2.0*f3 - h2*g2)*y1);

			double myRF;
			/*
			if (upper_limit_R < R){
				myRF = 0.0;
			}else*/
			if (R<Spe_PAO_RV[wan][0]){
			
				double df = 2.0*y2/h2*(3.0*f2 + h2*g1 + (2.0*f2 + h2*g1)*y2)
				+ y22*(2.0*f2 + h2*g1)/h2
				+ 2.0*y1/h2*(3.0*f3 - h2*g2 - (2.0*f3 - h2*g2)*y1)
				- y12*(2.0*f3 - h2*g2)/h2;

				if (L0==0){
				double a = 0.0;
				double b = 0.5*df/rm;
				double c = 0.0;
				double d = f - b*rm*rm;
				myRF = a*R*R*R + b*R*R + c*R + d;
				}

				else if (L0==1){
				double a = (rm*df - f)/(2.0*rm*rm*rm);
				double b = 0.0;
				double c = df - 3.0*a*rm*rm;
				double d = 0.0;
				myRF = a*R*R*R + b*R*R + c*R + d;
				}

				else{
				double b = (3.0*f - rm*df)/(rm*rm);
				double a = (f - b*rm*rm)/(rm*rm*rm);
				double c = 0.0;
				double d = 0.0;
				myRF = a*R*R*R + b*R*R + c*R + d;
				}
			/*}else if (upper_limit_R<R){
				myRF = 0.0;			*/
			}else{
				myRF = f;
			}

			RF[n0 * size_n + Nc] = myRF;
	    
		} /* (n0 < ai_max_LMul0) */
	  }/* n0 */
	  /* end of inlining of Get_Orbitals */
	
	}
	}/* Nc */

	//#pragma _NEC novector
    for (int Nc=OMPID*GridN_Atom[Gc_AN]/Nthrds; Nc<(OMPID+1)*GridN_Atom[Gc_AN]/Nthrds; Nc++){

		int i1 = 0;
      	for(int n0 = 0; n0 < ai_max_LMul0; ++n0){
		  	const int L0=ai_rf_L0[n0];
		  	const int Mul0=ai_rf_Mul0[n0];

	      	for (int M0=0; M0<=2*L0; M0++){
		    	Orbs_Grid[Mc_AN][Nc][i1] = RF[n0*size_n + Nc] * AF[(L0*L0+M0)*size_n + Nc];
	        	++i1;
		  	}		
	  	}
	}
    free(RF);
	free(AF);
	
		free(ai_rf_L0);
		free(ai_rf_Mul0);	  
	}else
#endif //VEC_NEC
	{


      double SH[Supported_MaxL*2+1][2];
      double dSHt[Supported_MaxL*2+1][2];
      double dSHp[Supported_MaxL*2+1][2];

      double** AF = (double**)malloc(sizeof(double*)*(List_YOUSO[25]+1));
      for (int i=0; i<(List_YOUSO[25]+1); i++){
	AF[i] = (double*)malloc(sizeof(double)*(2*(List_YOUSO[25]+1)+1));
      }

      /* get info. on OpenMP */ 

      const int OMPID = omp_get_thread_num();
      const int Nthrds = omp_get_num_threads();
      const int Nprocs = omp_get_num_procs();


      for (int Nc=OMPID*GridN_Atom[Gc_AN]/Nthrds; Nc<(OMPID+1)*GridN_Atom[Gc_AN]/Nthrds; Nc++){

	const int GNc = GridListAtom[Mc_AN][Nc]; 
	const int GRc = CellListAtom[Mc_AN][Nc];

	double Cxyz[4];
	Get_Grid_XYZ_inline(GNc,Cxyz);
	double x = Cxyz[1] + atv[GRc][1] - Gxyz[Gc_AN][1]; 
	double y = Cxyz[2] + atv[GRc][2] - Gxyz[Gc_AN][2]; 
	double z = Cxyz[3] + atv[GRc][3] - Gxyz[Gc_AN][3];


          /* Get_Orbitals(Cwan,x,y,z,Chi0); */
          /* start of inlining of Get_Orbitals */


          /* xyz2spherical(x,y,z,0.0,0.0,0.0,S_coordinate); */
          /* start of inlining of xyz2spherical */

	  const double Min_r = 10e-15;
	  double dum = x*x + y*y; 
	  double r = sqrt(dum + z*z);
	  double r1 = sqrt(dum);

	  const double dum1t = (r<fabs(z)) ? ((z>0.0)?1.0:-1.0) : z/r;
	  const double theta = (Min_r<=r) ? acos(dum1t) : 0.5*PI;

	  const double dum1p = (r1<fabs(y)) ? ((y>0.0)?1.0:-1.0) : y/r1;  
	  const double phi = (Min_r<=r1) ? ((0.0<=x) ? asin(dum1p) : PI - asin(dum1p)) : 0.0;



	  double R = r;
	  double Q = theta;
	  double P = phi;

	  /* end of inlining of xyz2spherical */

	  int po = 0;
	  int mp_min = 0;
	  int mp_max = Spe_Num_Mesh_PAO[wan] - 1;
	  int m;

	  if (Spe_PAO_RV[wan][Spe_Num_Mesh_PAO[wan]-1]<R){

	    int i1 = 0;
	  for (int L0=0; L0<=Spe_MaxL_Basis[wan]; L0++){
	    for (int Mul0=0; Mul0<Spe_Num_Basis[wan][L0]; Mul0++){
	      for (int M0=0; M0<=2*L0; M0++){
		    Orbs_Grid[Mc_AN][Nc][i1] = 0.0;
	        i1++;
		  }
	    }
	  }
	  
	  }
	  else{

		  
	    /* Angular */
	    const double siQ = sin(Q);
	    const double coQ = cos(Q);
	    const double siP = sin(P);
	    const double coP = cos(P);

	    for (int L0=0; L0<=Spe_MaxL_Basis[wan]; L0++){

	      if (L0==0){
		    AF[0][0] = 0.282094791773878;
	      }
	      else if (L0==1){
			const double dum = 0.48860251190292*siQ;
			AF[1][0] = dum*coP;
			AF[1][1] = dum*siP;
			AF[1][2] = 0.48860251190292*coQ;
	      }
	      else if (L0==2){
			const double dum1 = siQ*siQ;
			const double dum2 = 1.09254843059208*siQ*coQ;
			AF[2][0] = 0.94617469575756*coQ*coQ - 0.31539156525252;
			AF[2][1] = 0.54627421529604*dum1*(1.0 - 2.0*siP*siP);
			AF[2][2] = 1.09254843059208*dum1*siP*coP;
			AF[2][3] = dum2*coP;
			AF[2][4] = dum2*siP;
	      }
	      else if (L0==3){
			AF[3][0] = 0.373176332590116*(5.0*coQ*coQ*coQ - 3.0*coQ);
			AF[3][1] = 0.457045799464466*coP*siQ*(5.0*coQ*coQ - 1.0);
			AF[3][2] = 0.457045799464466*siP*siQ*(5.0*coQ*coQ - 1.0);
			AF[3][3] = 1.44530572132028*siQ*siQ*coQ*(coP*coP-siP*siP);
			AF[3][4] = 2.89061144264055*siQ*siQ*coQ*siP*coP;
			AF[3][5] = 0.590043589926644*siQ*siQ*siQ*(4.0*coP*coP*coP - 3.0*coP);
			AF[3][6] = 0.590043589926644*siQ*siQ*siQ*(3.0*siP - 4.0*siP*siP*siP);
	      }
	      else if (4<=L0){

			/* calculation of complex spherical harmonics functions */
			for(int m=-L0; m<=L0; m++){ 
			ComplexSH(L0,m,Q,P,SH[L0+m],dSHt[L0+m],dSHp[L0+m]);
			}

			/* transformation of complex to real */
			for (int i=0; i<(L0*2+1); i++){

			double sum0 = 0.0;
			double sum1 = 0.0; 
			for (int j=0; j<(L0*2+1); j++){
				sum0 += Comp2Real[L0][i][j].r*SH[j][0] - Comp2Real[L0][i][j].i*SH[j][1]; 
				sum1 += Comp2Real[L0][i][j].r*SH[j][1] + Comp2Real[L0][i][j].i*SH[j][0]; 
			}
			AF[L0][i] = sum0 + sum1; 
			}              

	      }
	    }/* L0 */
	  

		double rm;
		if (R<Spe_PAO_RV[wan][0]){
			m = 4;
			rm = Spe_PAO_RV[wan][m];
		}else{				
			do{
				m = (mp_min + mp_max)/2;
				if (Spe_PAO_RV[wan][m]<R)
					mp_min = m;
				else 
					mp_max = m;
			}
			while((mp_max-mp_min)!=1);
			m = mp_max;
			rm = R;
		}

	    double h1 = Spe_PAO_RV[wan][m-1] - Spe_PAO_RV[wan][m-2];
	    double h2 = Spe_PAO_RV[wan][m]   - Spe_PAO_RV[wan][m-1];
	    double h3 = Spe_PAO_RV[wan][m+1] - Spe_PAO_RV[wan][m];

		if (m==1){
			h1 = -(h2+h3);
		}
		else if (m==(Spe_Num_Mesh_PAO[wan]-1)){
			h3 = -(h1+h2);
		}

	    double x1 = rm - Spe_PAO_RV[wan][m-1];
	    double x2 = rm - Spe_PAO_RV[wan][m];
	    double y1 = x1/h2;
	    double y2 = x2/h2;
	    double y12 = y1*y1;
	    double y22 = y2*y2;
		
	    double dum = h1 + h2;
	    double dum1 = h1/h2/dum;
	    double dum2 = h2/h1/dum;
	    dum = h2 + h3;
	    double dum3 = h2/h3/dum;
	    double dum4 = h3/h2/dum;


	  	int i1 = 0;
	    for (int L0=0; L0<=Spe_MaxL_Basis[wan]; L0++){
	      for (int Mul0=0; Mul0<Spe_Num_Basis[wan][L0]; Mul0++){

			double f1 = Spe_PAO_RWF[wan][L0][Mul0][m-2];
			double f2 = Spe_PAO_RWF[wan][L0][Mul0][m-1];
			double f3 = Spe_PAO_RWF[wan][L0][Mul0][m];
			double f4 = Spe_PAO_RWF[wan][L0][Mul0][m+1];

			if (m==1){
				f1 = f4;
			}
			else if (m==(Spe_Num_Mesh_PAO[wan]-1)){
				f4 = f1;
			}

			double dum = f3 - f2;
			double g1 = dum*dum1 + (f2-f1)*dum2;
			double g2 = (f4-f3)*dum3 + dum*dum4;

			double f =  y22*(3.0*f2 + h2*g1 + (2.0*f2 + h2*g1)*y2)
			+ y12*(3.0*f3 - h2*g2 - (2.0*f3 - h2*g2)*y1);

			double myRF;

			if (R<Spe_PAO_RV[wan][0]){
			
				double df = 2.0*y2/h2*(3.0*f2 + h2*g1 + (2.0*f2 + h2*g1)*y2)
				+ y22*(2.0*f2 + h2*g1)/h2
				+ 2.0*y1/h2*(3.0*f3 - h2*g2 - (2.0*f3 - h2*g2)*y1)
				- y12*(2.0*f3 - h2*g2)/h2;

				if (L0==0){
				double a = 0.0;
				double b = 0.5*df/rm;
				double c = 0.0;
				double d = f - b*rm*rm;
				myRF = a*R*R*R + b*R*R + c*R + d;
				}

				else if (L0==1){
				double a = (rm*df - f)/(2.0*rm*rm*rm);
				double b = 0.0;
				double c = df - 3.0*a*rm*rm;
				double d = 0.0;
				myRF = a*R*R*R + b*R*R + c*R + d;
				}

				else{
				double b = (3.0*f - rm*df)/(rm*rm);
				double a = (f - b*rm*rm)/(rm*rm*rm);
				double c = 0.0;
				double d = 0.0;
				myRF = a*R*R*R + b*R*R + c*R + d;
				}

			}else{
				myRF = f;
			}




	      for (int M0=0; M0<=2*L0; M0++){
		Orbs_Grid[Mc_AN][Nc][i1] = myRF*AF[L0][M0];
	    i1++;
		  }
	    }/* Mul0 */
	  }/* L0 */
	  /* end of inlining of Get_Orbitals */
	
	}
	}/* Nc */




      for (int i=0; i<(List_YOUSO[25]+1); i++){
	    free(AF[i]);
      }
      free(AF);

	}
	} /* if (Cnt_kind==0) */
	else{
		
      double* Chi0 = (double*)malloc(sizeof(double)*List_YOUSO[7]);

      /* get info. on OpenMP */ 

      const int OMPID = omp_get_thread_num();
      const int Nthrds = omp_get_num_threads();
      const int Nprocs = omp_get_num_procs();

      for (int Nc=OMPID*GridN_Atom[Gc_AN]/Nthrds; Nc<(OMPID+1)*GridN_Atom[Gc_AN]/Nthrds; Nc++){

		const int GNc = GridListAtom[Mc_AN][Nc]; 
		const int GRc = CellListAtom[Mc_AN][Nc];

		double Cxyz[4];
		Get_Grid_XYZ_inline(GNc,Cxyz);
		double x = Cxyz[1] + atv[GRc][1] - Gxyz[Gc_AN][1]; 
		double y = Cxyz[2] + atv[GRc][2] - Gxyz[Gc_AN][2]; 
		double z = Cxyz[3] + atv[GRc][3] - Gxyz[Gc_AN][3];

          Get_Cnt_Orbitals(Mc_AN,x,y,z,Chi0);
		  for(int i = 0; i < NO0; ++i){
				Orbs_Grid[Mc_AN][Nc][i]=Chi0[i];
			}


      } /* Nc */
	  free(Chi0);

	}


    } /* #pragma omp parallel */

    dtime(&Etime_atom);
    time_per_atom[Gc_AN] += Etime_atom - Stime_atom;
  }

  /****************************************************
     Calculate Orbs_Grid_FNAN
  ****************************************************/

  for (int Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){

    const int Gc_AN = M2G[Mc_AN];    

    for (int h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){

      const int Gh_AN = natn[Gc_AN][h_AN];

      if (G2ID[Gh_AN]!=myid){

        const int Mh_AN = F_G2M[Gh_AN];
        const int Rnh = ncn[Gc_AN][h_AN];
        const int Hwan = WhatSpecies[Gh_AN];

		int NO1;
        if (Cnt_kind==0)  NO1 = Spe_Total_NO[Hwan];
        else              NO1 = Spe_Total_CNO[Hwan];

#pragma omp parallel shared(List_YOUSO,Orbs_Grid_FNAN,NO1,Mh_AN,Hwan,Cnt_kind,Rnh,Gh_AN,Gxyz,atv,NumOLG,Mc_AN,h_AN,GListTAtoms1,GridListAtom,CellListAtom) private(OMPID,Nthrds,Nprocs)
        {

	  /* Radial */
	  //int mp_min,mp_max,m,po,wan;
	  
          /* allocation of arrays */
    if(Cnt_kind == 0){
	  const int wan = Hwan;

	  int ai_max_LMul0 = 0;        
	  for (int L0=0; L0<=Spe_MaxL_Basis[wan]; L0++){
		ai_max_LMul0 += Spe_Num_Basis[wan][L0];		  
	  }
        

#ifdef VEC_NEC
      if((ai_max_LMul0 <= 15) && (Spe_MaxL_Basis[wan] < 4)){

#ifdef AITURN_20200914

/* AI_TURN */
      int* ai_rf_L0 = (int*)malloc(sizeof(int)*ai_max_LMul0);
      int* ai_rf_Mul0 = (int*)malloc(sizeof(int)*ai_max_LMul0);
	  
      
      {
        int n0 = 0;
        for (int L0=0; L0<=Spe_MaxL_Basis[wan]; L0++){
	        for (int Mul0=0; Mul0<Spe_Num_Basis[wan][L0]; Mul0++){
				ai_rf_L0[n0] = L0;
				ai_rf_Mul0[n0] = Mul0;
            n0++;
            
          }
        }        
	  }

#endif

#ifdef AITURN_20200914b
		int rank_for_m = 1;
		{
			int bi = Spe_Num_Mesh_PAO[wan];
			while (bi > 0) {
				rank_for_m ++;
				bi >>= 1;
			}
		}

#endif

	  /* get info. on OpenMP */ 

	  const int OMPID = omp_get_thread_num();
	  const int Nthrds = omp_get_num_threads();
	  const int Nprocs = omp_get_num_procs();

	  const int size_n = NumOLG[Mc_AN][h_AN];
	  double* RF = (double*)malloc(sizeof(double) * size_n * ai_max_LMul0);
	  double* AF = (double*)malloc(sizeof(double) * size_n * 16);

	  
	  const double upper_limit_R = Spe_PAO_RV[wan][Spe_Num_Mesh_PAO[wan]-1];
      
	  #pragma _NEC vector
	  #pragma _NEC ivdep
	  for (int Nog=OMPID*NumOLG[Mc_AN][h_AN]/Nthrds; Nog<(OMPID+1)*NumOLG[Mc_AN][h_AN]/Nthrds; Nog++){

	    const int Nc = GListTAtoms1[Mc_AN][h_AN][Nog];
	    const int GNc = GridListAtom[Mc_AN][Nc];
	    const int GRc = CellListAtom[Mc_AN][Nc]; 

	  	double Cxyz0_1, Cxyz0_2, Cxyz0_3; 
	    //Get_Grid_XYZ_inline(GNc,Cxyz0);
		{
			int n1,n2,n3;

			n1 = GNc/(Ngrid2*Ngrid3);
			n2 = (GNc - n1*(Ngrid2*Ngrid3))/Ngrid3;
			n3 = GNc - n1*(Ngrid2*Ngrid3) - n2*Ngrid3;

			Cxyz0_1 = (double)n1*gtv[1][1] + (double)n2*gtv[2][1]
					+ (double)n3*gtv[3][1] + Grid_Origin[1];
			Cxyz0_2 = (double)n1*gtv[1][2] + (double)n2*gtv[2][2]
					+ (double)n3*gtv[3][2] + Grid_Origin[2];
			Cxyz0_3 = (double)n1*gtv[1][3] + (double)n2*gtv[2][3]
					+ (double)n3*gtv[3][3] + Grid_Origin[3];
		}

	    double x = Cxyz0_1 + atv[GRc][1] - Gxyz[Gh_AN][1] - atv[Rnh][1];
	    double y = Cxyz0_2 + atv[GRc][2] - Gxyz[Gh_AN][2] - atv[Rnh][2];
	    double z = Cxyz0_3 + atv[GRc][3] - Gxyz[Gh_AN][3] - atv[Rnh][3];

	    

              /* Get_Orbitals(Hwan,x,y,z,Chi0); */
              /* start of inlining of Get_Orbitals */

              

	      /* xyz2spherical(x,y,z,0.0,0.0,0.0,S_coordinate); */
              /* start of inlining of xyz2spherical */

	      const double Min_r = 10e-15;
	      double dum = x*x + y*y; 
	      double r = sqrt(dum + z*z);
	      double r1 = sqrt(dum);
		  
	      
		  const double dum1t = (r<fabs(z)) ? ((z>0.0)?1.0:-1.0) : z/r;
		  const double theta = (Min_r<=r) ? acos(dum1t) : 0.5*PI;

		  const double dum1p = (r1<fabs(y)) ? ((y>0.0)?1.0:-1.0) : y/r1;  
		  const double phi = (Min_r<=r1) ? ((0.0<=x) ? asin(dum1p) : PI - asin(dum1p)) : 0.0;

	      double R = r;
	      double Q = theta;
	      double P = phi;

              /* end of inlining of xyz2spherical */

	      int po = 0;
	      int mp_min = 0;
	      int mp_max = Spe_Num_Mesh_PAO[wan] - 1;
		  int m;

	      if (Spe_PAO_RV[wan][Spe_Num_Mesh_PAO[wan]-1]<R){

	        #pragma _NEC unroll(15) 
	        for(int i1 = 0; i1 < 15; ++i1){
				if(i1 < ai_max_LMul0){
					//Orbs_Grid_FNAN[Mc_AN][h_AN][Nog][i1] = 0.0;		          
					RF[i1 * size_n + Nog] = 0.0;
		        }
	        }

	      }
	      else
		  {
			
			/* Angular */
			const double siQ = sin(Q);
			const double coQ = cos(Q);
			const double siP = sin(P);
			const double coP = cos(P);
			{
				AF[(0+0)*size_n + Nog] = 0.282094791773878;
			
				const double dum = 0.48860251190292*siQ;
				AF[(1+0)*size_n + Nog] = dum*coP;
				AF[(1+1)*size_n + Nog] = dum*siP;
				AF[(1+2)*size_n + Nog] = 0.48860251190292*coQ;
			
				const double dum1 = siQ*siQ;
				const double dum2 = 1.09254843059208*siQ*coQ;
				AF[(4+0)*size_n + Nog] = 0.94617469575756*coQ*coQ - 0.31539156525252;
				AF[(4+1)*size_n + Nog] = 0.54627421529604*dum1*(1.0 - 2.0*siP*siP);
				AF[(4+2)*size_n + Nog] = 1.09254843059208*dum1*siP*coP;
				AF[(4+3)*size_n + Nog] = dum2*coP;
				AF[(4+4)*size_n + Nog] = dum2*siP;
			
				AF[(9+0)*size_n + Nog] = 0.373176332590116*(5.0*coQ*coQ*coQ - 3.0*coQ);
				AF[(9+1)*size_n + Nog] = 0.457045799464466*coP*siQ*(5.0*coQ*coQ - 1.0);
				AF[(9+2)*size_n + Nog] = 0.457045799464466*siP*siQ*(5.0*coQ*coQ - 1.0);
				AF[(9+3)*size_n + Nog] = 1.44530572132028*siQ*siQ*coQ*(coP*coP-siP*siP);
				AF[(9+4)*size_n + Nog] = 2.89061144264055*siQ*siQ*coQ*siP*coP;
				AF[(9+5)*size_n + Nog] = 0.590043589926644*siQ*siQ*siQ*(4.0*coP*coP*coP - 3.0*coP);
				AF[(9+6)*size_n + Nog] = 0.590043589926644*siQ*siQ*siQ*(3.0*siP - 4.0*siP*siP*siP);
			}


			double rm;
		    if (R<Spe_PAO_RV[wan][0]){
			  
				m = 4;
				rm = Spe_PAO_RV[wan][m];
			}else{
				
#ifdef AITURN_20200914b
			
			#pragma _NEC loop_count(16)
			#pragma _NEC unroll(16)
			//#pragma _NEC novector
			for(int rank = 0; rank < 16; ++rank){
				if(rank < rank_for_m){
				if((mp_max - mp_min) > 1){
					m = (mp_min + mp_max) / 2;
					if (Spe_PAO_RV[wan][m] < R) mp_min = m;
					else mp_max = m;
				}
				}
			}
			m = mp_max;
#else		
				do{
				m = (mp_min + mp_max)/2;
				if (Spe_PAO_RV[wan][m]<R)
					mp_min = m;
				else 
					mp_max = m;
				}
				while((mp_max-mp_min)!=1);
				m = mp_max;
#endif
				rm = R;
			}




		double h1 = Spe_PAO_RV[wan][m-1] - Spe_PAO_RV[wan][m-2];
		double h2 = Spe_PAO_RV[wan][m]   - Spe_PAO_RV[wan][m-1];
		double h3 = Spe_PAO_RV[wan][m+1] - Spe_PAO_RV[wan][m];
		if (m==1){
			h1 = -(h2+h3);
		}
		else if (m==(Spe_Num_Mesh_PAO[wan]-1)){
			h3 = -(h1+h2);
		}
		
		double x1 = rm - Spe_PAO_RV[wan][m-1];
		double x2 = rm - Spe_PAO_RV[wan][m];
		double y1 = x1/h2;
		double y2 = x2/h2;
		double y12 = y1*y1;
		double y22 = y2*y2;

		double dum = h1 + h2;
		double dum1 = h1/h2/dum;
		double dum2 = h2/h1/dum;
		dum = h2 + h3;
		double dum3 = h2/h3/dum;
		double dum4 = h3/h2/dum;

		#pragma _NEC unroll(15) 
		for(int n0 = 0; n0 < 15; ++n0){
		  if(n0 < ai_max_LMul0){
			const int L0 = ai_rf_L0[n0];
			const int Mul0 = ai_rf_Mul0[n0];
			

		    double f1 = (m!=1) ? Spe_PAO_RWF[wan][L0][Mul0][m-2] : Spe_PAO_RWF[wan][L0][Mul0][m+1];
		    double f2 = Spe_PAO_RWF[wan][L0][Mul0][m-1];
		    double f3 = Spe_PAO_RWF[wan][L0][Mul0][m];
		    double f4 = (m!=(Spe_Num_Mesh_PAO[wan]-1)) ? Spe_PAO_RWF[wan][L0][Mul0][m+1] : Spe_PAO_RWF[wan][L0][Mul0][m-2];


		    double dum = f3 - f2;
		    double g1 = dum*dum1 + (f2-f1)*dum2;
		    double g2 = (f4-f3)*dum3 + dum*dum4;

		    double f =  y22*(3.0*f2 + h2*g1 + (2.0*f2 + h2*g1)*y2)
		      + y12*(3.0*f3 - h2*g2 - (2.0*f3 - h2*g2)*y1);

			double myRF;

			if (R<Spe_PAO_RV[wan][0]){
				double df = 2.0*y2/h2*(3.0*f2 + h2*g1 + (2.0*f2 + h2*g1)*y2)
				+ y22*(2.0*f2 + h2*g1)/h2
				+ 2.0*y1/h2*(3.0*f3 - h2*g2 - (2.0*f3 - h2*g2)*y1)
				- y12*(2.0*f3 - h2*g2)/h2;

				if (L0==0){
				double a = 0.0;
				double b = 0.5*df/rm;
				double c = 0.0;
				double d = f - b*rm*rm;
				myRF = a*R*R*R + b*R*R + c*R + d;
				}

				else if (L0==1){
				double a = (rm*df - f)/(2.0*rm*rm*rm);
				double b = 0.0;
				double c = df - 3.0*a*rm*rm;
				double d = 0.0;
				myRF = a*R*R*R + b*R*R + c*R + d;
				}

				else{
				double b = (3.0*f - rm*df)/(rm*rm);
				double a = (f - b*rm*rm)/(rm*rm*rm);
				double c = 0.0;
				double d = 0.0;
				myRF = a*R*R*R + b*R*R + c*R + d;
				}
			/*}else if (upper_limit_R<R){
				myRF = 0.0;			*/
			}else{
				myRF = f;
			}
		    
			RF[n0 * size_n + Nog] = myRF;
		    
		  }/* (n0 < ai_max_LMul0) */
	    }/* n0 */

              /* end of inlining of Get_Orbitals */
		  }
		  }/* Nog */

		  	
			for (int Nog=OMPID*NumOLG[Mc_AN][h_AN]/Nthrds; Nog<(OMPID+1)*NumOLG[Mc_AN][h_AN]/Nthrds; Nog++){
				
				int i1 = 0;
      			for(int n0 = 0; n0 < ai_max_LMul0; ++n0){
					const int L0=ai_rf_L0[n0];
					const int Mul0=ai_rf_Mul0[n0];

					for (int M0=0; M0<=2*L0; M0++){
						Orbs_Grid_FNAN[Mc_AN][h_AN][Nog][i1] = RF[n0*size_n + Nog] * AF[(L0*L0+M0)*size_n + Nog];
						++i1;
					}					
				}
			}
			free(RF);
			free(AF);
		  
          /* freeing of arrays */

      		free(ai_rf_L0);
      		free(ai_rf_Mul0);	  
	  }else
#endif  //VEC_NEC
	  {

	  double SH[Supported_MaxL*2+1][2];
	  double dSHt[Supported_MaxL*2+1][2];
	  double dSHp[Supported_MaxL*2+1][2];


	  double** AF = (double**)malloc(sizeof(double*)*(List_YOUSO[25]+1));
	  for (int i=0; i<(List_YOUSO[25]+1); i++){
	    AF[i] = (double*)malloc(sizeof(double)*(2*(List_YOUSO[25]+1)+1));
	  }

	  /* get info. on OpenMP */ 

	  const int OMPID = omp_get_thread_num();
	  const int Nthrds = omp_get_num_threads();
	  const int Nprocs = omp_get_num_procs();

	
	  for (int Nog=OMPID*NumOLG[Mc_AN][h_AN]/Nthrds; Nog<(OMPID+1)*NumOLG[Mc_AN][h_AN]/Nthrds; Nog++){

	    const int Nc = GListTAtoms1[Mc_AN][h_AN][Nog];
	    const int GNc = GridListAtom[Mc_AN][Nc];
	    const int GRc = CellListAtom[Mc_AN][Nc]; 

	  	double Cxyz0[4]; 
	    Get_Grid_XYZ_inline(GNc,Cxyz0);

	    double x = Cxyz0[1] + atv[GRc][1] - Gxyz[Gh_AN][1] - atv[Rnh][1];
	    double y = Cxyz0[2] + atv[GRc][2] - Gxyz[Gh_AN][2] - atv[Rnh][2];
	    double z = Cxyz0[3] + atv[GRc][3] - Gxyz[Gh_AN][3] - atv[Rnh][3];

	    

              /* Get_Orbitals(Hwan,x,y,z,Chi0); */
              /* start of inlining of Get_Orbitals */

              

	      /* xyz2spherical(x,y,z,0.0,0.0,0.0,S_coordinate); */
              /* start of inlining of xyz2spherical */

	      const double Min_r = 10e-15;
	      double dum = x*x + y*y; 
	      double r = sqrt(dum + z*z);
	      double r1 = sqrt(dum);
		  
	      
		  const double dum1t = (r<fabs(z)) ? ((z>0.0)?1.0:-1.0) : z/r;
		  const double theta = (Min_r<=r) ? acos(dum1t) : 0.5*PI;

		  const double dum1p = (r1<fabs(y)) ? ((y>0.0)?1.0:-1.0) : y/r1;  
		  const double phi = (Min_r<=r1) ? ((0.0<=x) ? asin(dum1p) : PI - asin(dum1p)) : 0.0;

	      double R = r;
	      double Q = theta;
	      double P = phi;

              /* end of inlining of xyz2spherical */

	      int po = 0;
	      int mp_min = 0;
	      int mp_max = Spe_Num_Mesh_PAO[wan] - 1;
		  int m;

	      if (Spe_PAO_RV[wan][Spe_Num_Mesh_PAO[wan]-1]<R){

	
	        /* Chi0 */  
	        int i1 = 0;
	        for (int L0=0; L0<=Spe_MaxL_Basis[wan]; L0++){
		      for (int Mul0=0; Mul0<Spe_Num_Basis[wan][L0]; Mul0++){
		        for (int M0=0; M0<=2*L0; M0++){
		          Orbs_Grid_FNAN[Mc_AN][h_AN][Nog][i1] = 0.0;
		          i1++;
		        }
		      }
	        }

	      }
	      else{

			/* Angular */
			const double siQ = sin(Q);
			const double coQ = cos(Q);
			const double siP = sin(P);
			const double coP = cos(P);

			for (int L0=0; L0<=Spe_MaxL_Basis[wan]; L0++){

				if (L0==0){
					AF[0][0] = 0.282094791773878;
				}
				else if (L0==1){
					const double dum = 0.48860251190292*siQ;
					AF[1][0] = dum*coP;
					AF[1][1] = dum*siP;
					AF[1][2] = 0.48860251190292*coQ;
				}
				else if (L0==2){
					const double dum1 = siQ*siQ;
					const double dum2 = 1.09254843059208*siQ*coQ;
					AF[2][0] = 0.94617469575756*coQ*coQ - 0.31539156525252;
					AF[2][1] = 0.54627421529604*dum1*(1.0 - 2.0*siP*siP);
					AF[2][2] = 1.09254843059208*dum1*siP*coP;
					AF[2][3] = dum2*coP;
					AF[2][4] = dum2*siP;
				}

				else if (L0==3){
					AF[3][0] = 0.373176332590116*(5.0*coQ*coQ*coQ - 3.0*coQ);
					AF[3][1] = 0.457045799464466*coP*siQ*(5.0*coQ*coQ - 1.0);
					AF[3][2] = 0.457045799464466*siP*siQ*(5.0*coQ*coQ - 1.0);
					AF[3][3] = 1.44530572132028*siQ*siQ*coQ*(coP*coP-siP*siP);
					AF[3][4] = 2.89061144264055*siQ*siQ*coQ*siP*coP;
					AF[3][5] = 0.590043589926644*siQ*siQ*siQ*(4.0*coP*coP*coP - 3.0*coP);
					AF[3][6] = 0.590043589926644*siQ*siQ*siQ*(3.0*siP - 4.0*siP*siP*siP);
				}

				else if (4<=L0){

					/* calculation of complex spherical harmonics functions */
					for(int m=-L0; m<=L0; m++){ 
					ComplexSH(L0,m,Q,P,SH[L0+m],dSHt[L0+m],dSHp[L0+m]);
					}

					/* transformation of complex to real */
					for (int i=0; i<(L0*2+1); i++){

					double sum0 = 0.0;
					double sum1 = 0.0; 
					for (int j=0; j<(L0*2+1); j++){
					sum0 += Comp2Real[L0][i][j].r*SH[j][0] - Comp2Real[L0][i][j].i*SH[j][1]; 
					sum1 += Comp2Real[L0][i][j].r*SH[j][1] + Comp2Real[L0][i][j].i*SH[j][0]; 
					}
					AF[L0][i] = sum0 + sum1; 
					}              

				}
			}/* L0 */

			double rm;
		    if (R<Spe_PAO_RV[wan][0]){
			  
				m = 4;
				rm = Spe_PAO_RV[wan][m];
			}else{
						
				do{
				m = (mp_min + mp_max)/2;
				if (Spe_PAO_RV[wan][m]<R)
					mp_min = m;
				else 
					mp_max = m;
				}
				while((mp_max-mp_min)!=1);
				m = mp_max;
				rm = R;
			}




		double h1 = Spe_PAO_RV[wan][m-1] - Spe_PAO_RV[wan][m-2];
		double h2 = Spe_PAO_RV[wan][m]   - Spe_PAO_RV[wan][m-1];
		double h3 = Spe_PAO_RV[wan][m+1] - Spe_PAO_RV[wan][m];

		if (m==1){
		    h1 = -(h2+h3);
		}
		else if (m==(Spe_Num_Mesh_PAO[wan]-1)){
		    h3 = -(h1+h2);
		}

		double x1 = rm - Spe_PAO_RV[wan][m-1];
		double x2 = rm - Spe_PAO_RV[wan][m];
		double y1 = x1/h2;
		double y2 = x2/h2;
		double y12 = y1*y1;
		double y22 = y2*y2;

		double dum = h1 + h2;
		double dum1 = h1/h2/dum;
		double dum2 = h2/h1/dum;
		dum = h2 + h3;
		double dum3 = h2/h3/dum;
		double dum4 = h3/h2/dum;

        int j = 0;
		for (int L0=0; L0<=Spe_MaxL_Basis[wan]; L0++){
		  for (int Mul0=0; Mul0<Spe_Num_Basis[wan][L0]; Mul0++){

		    double f1 = Spe_PAO_RWF[wan][L0][Mul0][m-2];
		    double f2 = Spe_PAO_RWF[wan][L0][Mul0][m-1];
		    double f3 = Spe_PAO_RWF[wan][L0][Mul0][m];
		    double f4 = Spe_PAO_RWF[wan][L0][Mul0][m+1];

		    if (m==1){
		      f1 = f4;
		    }
		    else if (m==(Spe_Num_Mesh_PAO[wan]-1)){
		      f4 = f1;
		    }

		    double dum = f3 - f2;
		    double g1 = dum*dum1 + (f2-f1)*dum2;
		    double g2 = (f4-f3)*dum3 + dum*dum4;

		    double f =  y22*(3.0*f2 + h2*g1 + (2.0*f2 + h2*g1)*y2)
		      + y12*(3.0*f3 - h2*g2 - (2.0*f3 - h2*g2)*y1);

			double myRF;

			if (R<Spe_PAO_RV[wan][0]){
				double df = 2.0*y2/h2*(3.0*f2 + h2*g1 + (2.0*f2 + h2*g1)*y2)
				+ y22*(2.0*f2 + h2*g1)/h2
				+ 2.0*y1/h2*(3.0*f3 - h2*g2 - (2.0*f3 - h2*g2)*y1)
				- y12*(2.0*f3 - h2*g2)/h2;

				if (L0==0){
				double a = 0.0;
				double b = 0.5*df/rm;
				double c = 0.0;
				double d = f - b*rm*rm;
				myRF = a*R*R*R + b*R*R + c*R + d;
				}

				else if (L0==1){
				double a = (rm*df - f)/(2.0*rm*rm*rm);
				double b = 0.0;
				double c = df - 3.0*a*rm*rm;
				double d = 0.0;
				myRF = a*R*R*R + b*R*R + c*R + d;
				}

				else{
				double b = (3.0*f - rm*df)/(rm*rm);
				double a = (f - b*rm*rm)/(rm*rm*rm);
				double c = 0.0;
				double d = 0.0;
				myRF = a*R*R*R + b*R*R + c*R + d;
				}
			}else{
				myRF = f;
			}
		    
		    for (int M0=0; M0<=2*L0; M0++){
		    
		      Orbs_Grid_FNAN[Mc_AN][h_AN][Nog][j] = myRF*AF[L0][M0];
			  j++;
		    }
		  }/* Mul0 */
	    }/* L0 */

              /* end of inlining of Get_Orbitals */
		  }
		  }/* Nog */
		  
          /* freeing of arrays */



			for (int i=0; i<(List_YOUSO[25]+1); i++){
				free(AF[i]);
			}
			free(AF);

	  }
	    } /* if (Cnt_kind==0) */

	    else{
		  double* Chi0 = (double*)malloc(sizeof(double)*List_YOUSO[7]);

		  /* get info. on OpenMP */ 

		  const int OMPID = omp_get_thread_num();
		  const int Nthrds = omp_get_num_threads();
		  const int Nprocs = omp_get_num_procs();
				
		  for (int Nog=OMPID*NumOLG[Mc_AN][h_AN]/Nthrds; Nog<(OMPID+1)*NumOLG[Mc_AN][h_AN]/Nthrds; Nog++){

			const int Nc = GListTAtoms1[Mc_AN][h_AN][Nog];
			const int GNc = GridListAtom[Mc_AN][Nc];
			const int GRc = CellListAtom[Mc_AN][Nc]; 

			double Cxyz0[4]; 
			Get_Grid_XYZ_inline(GNc,Cxyz0);

			double x = Cxyz0[1] + atv[GRc][1] - Gxyz[Gh_AN][1] - atv[Rnh][1];
			double y = Cxyz0[2] + atv[GRc][2] - Gxyz[Gh_AN][2] - atv[Rnh][2];
			double z = Cxyz0[3] + atv[GRc][3] - Gxyz[Gh_AN][3] - atv[Rnh][3];

			
            Get_Cnt_Orbitals(Mh_AN,x,y,z,Chi0);
			for(int j = 0; j < NO1; ++j){
				Orbs_Grid_FNAN[Mc_AN][h_AN][Nog][j]=Chi0[j];
			}
		  } /* Nog */
					
          /* freeing of arrays */
			free(Chi0);

		}



        } 
      }
    }
  }

  /* time */
  dtime(&TEtime);
  time0 = TEtime - TStime;

  return time0;
}
