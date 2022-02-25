Modification list for NEC SX-Aurora TSUBASA based on OpenMX 3.9.2
The benchmark is reported as:
Atsushi M. Ito, Arimichi Takayama, Osamu Watanabe, Vijendra Singh, Shubham Tyagi, and Shashank S. Singh, "Tuning of Density Functional Theory Simulation on Vector Processor System - Plasma Simulator Raijin -" Plasma and Fusion Research, Rapid Communications, 15 (2020) 1203085.


Band_DFT_Col.c:
- bug-fix of the use of scalapack. The variables "alpha" and "beta" should be dcomplex.
- vectorization for loop "i in tnoA*tnoB".
- Inline expansion of FermiFunc() to remove "if" branch in the called function.
- the unroll of inner loop

Band_DFT_NonCol.c:
Cluster_DFT_NonCol.c:
- bug-fix of the use of scalapack. The variables "alpha" and "beta" should be dcomplex.

Force.c
- manual inline expansion of Get_dOrbitals and related function
- the change of algorism to find the position that is index indicated by radius in 1D array
- the vectorization of loop and the unroll of inner loop

FT_PAO.c
- note: The algorism can be changed into more high speed way to find the highest bit rank_for_m
- change the algorism for vectorization of "RadialF"

Get_Cnt_dOrbitals.c
- static variable is removed because it obstructs inline expansion

openmx_common.c
- the vectorization of "Associated_Legendre"

openmx.c
- MPI_Barrier is required to bottle neck of memory allocation in "truncation", which is problem only in NEC.

Poisson_ESM.c
Poisson.c
- using fftw3 library for NEC

RadialF_nec.c
- modified function called by FT_PAO.c

RF_BesselF.c
- modified function "RF_BesselF_ByIndex" is added, which is called by "SetNonlocal.c" and "Set_OLP_Kin.c"

Set_Density_Grid.c
Set_Hamiltonian.c
Set_Nonlocal.c
Set_OLP_Kin.c
Set_Orbitals_Grid.c
Set_ProExpn_VNA.c
- vectorization

Set_XC_Grid.c
- change the place of openmp directive.

Spherical_Bessel.c
- vectorization

Stress.c
- vectorization
- note: The algorism can be changed into more high speed way to find the highest bit rank_for_m

TRAN_CDen_Main.c
TRAN_Poisson.c
TRAN_Set_Electrode_Grid.c
- using fftw3 library for NEC

truncation.c
- bug-fix: the zero clear of "DS_NL"
- speedup of the set of grid infomation for "GRh1"

elpa1_merge_systems_real_template.F90
- bugfix: wrong description across multiple lines

makefile
- change compile commands and options
- explicitly indicate files of inline functions

