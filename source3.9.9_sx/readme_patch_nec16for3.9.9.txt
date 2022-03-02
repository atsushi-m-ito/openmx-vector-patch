#Patch for Vector Processors, vec16for3.9.9
#The list of files modified in nec-vec patch(vec14for3.9.2)
#The 2nd column is modification in official 3.9.9 from 3.9.6：
#   req(existing in nec-vec patch, and modified in official 3.9.9), 
#   nothing(existing in nec-vec patch, and not modified in official 3.9.9)
#   accept(nec-vec patch is accepted by official)
#The 3rd column is modification into 3.9.9: ok is complete
#FileName	req/nothing/accept      ok/none
Band_DFT_Col.c	        nothing
Band_DFT_NonCol.c	    nothing
Cluster_DFT_NonCol.c	nothing
FT_PAO.c	            nothing
Force.c	                nothing
Gaunt.h	                nothing
Get_Cnt_dOrbitals.c	    nothing
Poisson.c	            nothing
Poisson_ESM.c	        nothing
RF_BesselF.c	        nothing
RadialF_nec.c	        nothing
Set_Density_Grid.c	    nothing
Set_Hamiltonian.c	    nothing
Set_Nonlocal.c	        nothing
Set_OLP_Kin.c	        nothing
Set_Orbitals_Grid.c	    nothing
Set_ProExpn_VNA.c	    nothing
Set_XC_Grid.c	        nothing
Spherical_Bessel.c	    nothing
Stress.c	            nothing
TRAN_CDen_Main.c	    nothing
TRAN_Poisson.c	        nothing
TRAN_Set_Electrode_Grid.c	nothing
makefile	            req   ok   (changed point is a comment for the compilation on WSL)
openmx.c	            req   (confuse see OpenMX forum)
openmx_common.c	        nothing
truncation.c	        nothing
elpa1_merge_systems_real_template.F90	nothing


#Patch for Vector Processors, vec16for3.9.9
1) the bug fix:
-ftrace flag is needed until vec15for3.9.6, but this is bug. In this version, -ftrace is not needed.

Fixed files:
Set_OLP_Kin.c

2) tuning:
The calculation when MD.type = opt is improved.
The time step (SD_scaling) is changed by compare the direction of the present forces and the previous forces.

Fixed files:
MD_pac.c

3) the bug fix:
points of memory leak are fixed:

Fixed files:
Force.c       (only for NEC)
Set_Density_Grid.c  (only for NEC)
Stress.c   (only for NEC)
Free_Arrays.c   (in common with origin)
Input_std.c     (in common with origin)
