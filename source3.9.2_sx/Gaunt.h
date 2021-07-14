/**********************************************************************
  Gaunt.h:

     Gaunt.h is pre-calculation of data array calculated in Gaunt.c.
     This is used by Set_ProExpn_VNA.c, Set_OLP_Kin.c, and Set_Nonlocal.c

  Log of Gaunt.h:

     29/June/2021  Released by A.M.Ito

***********************************************************************/

static inline double* MakeGaunt(int L0_limit_sq, int L1_limit_sq, int Lw_limit_sq){
  double* Gaunt_results = (double*)malloc(sizeof(double) * L0_limit_sq*L1_limit_sq*Lw_limit_sq);
  for(int lm0 = 0; lm0 < L0_limit_sq; ++lm0){
    const int L0 = (int)sqrt(0.1+(double)lm0);
    const int M0 = lm0 - L0*L0 - L0;
    for(int lm1 = 0; lm1 < L1_limit_sq; ++lm1){
      const int L1 = (int)sqrt(0.1+(double)lm1);
      const int M1 = lm1 - L1*L1 - L1;
      for(int lmw = 0; lmw < Lw_limit_sq; ++lmw){
        const int Lw = (int)sqrt(0.1+(double)lmw);
        const int Mw = lmw - Lw*Lw - Lw;
        const index = (lm0 * L1_limit_sq + lm1) * Lw_limit_sq + lmw;
        //Gaunt_results[index] = Gaunt_inline(L0, M0, L1, M1, Lw, Mw);
        Gaunt_results[index] = Gaunt(L0, M0, L1, M1, Lw, Mw);
      }
    }  
  }
  return Gaunt_results;
}

static inline double GetGaunt(int l0,  int m0,
             int l1, int m1,
             int lw, int mw, const double* Gaunt_list, int L1_limit_sq, int Lw_limit_sq)
{
  return Gaunt_list[( (l0*l0+m0+l0) * L1_limit_sq + (l1*l1+m1+l1) ) * Lw_limit_sq + (lw*lw+mw+lw)];
}

