#Patch for Vector Processors, vec15for3.9.6
#The list of files modified in nec-vec patch(vec14for3.9.2)
#The 2nd column is modification in official 3.9.6：
#   req(existing in nec-vec patch, and modified in official 3.9.6), 
#   nothing(existing in nec-vec patch, and not modified in official 3.9.6)
#   accept(nec-vec patch is accepted by official)
#The 3rd column is modification into 3.9.6: ok is complete
#FileName	req/nothing/accept      ok/none
Band_DFT_Col.c	        req     ok
Band_DFT_NonCol.c	    accept
Cluster_DFT_NonCol.c	accept
FT_PAO.c	            nothing
Force.c	                req     ok
Gaunt.h	                nothing
Get_Cnt_dOrbitals.c	    nothing
Poisson.c	            req     ok(only adding fftw header for NEC into the file based on official3.9.6)
Poisson_ESM.c	        nothing
RF_BesselF.c	        nothing
RadialF_nec.c	        nothing
Set_Density_Grid.c	    req     ok
Set_Hamiltonian.c	    req     ok(changed point is only space)
Set_Nonlocal.c	        nothing
Set_OLP_Kin.c	        req     ok(changed point is only space)
Set_Orbitals_Grid.c	    nothing
Set_ProExpn_VNA.c	    nothing
Set_XC_Grid.c	        req     ok
Spherical_Bessel.c	    nothing
Stress.c	            nothing
TRAN_CDen_Main.c	    req     ok(only adding fftw header for NEC into the file based on official3.9.6)
TRAN_Poisson.c	        req     ok(only adding fftw header for NEC into the file based on official3.9.6)
TRAN_Set_Electrode_Grid.c	req ok(only adding fftw header for NEC into the file based on official3.9.6)
makefile	            req ok
openmx.c	            req ok(only adding MPI_Barrier in line 654 into the file based on official3.9.6)
openmx_common.c	        req ok(Associated_Legendre comes back to 3.9.1 (accepted our report))
truncation.c	        req ok(only adding the parts of _NES into the file based on official3.9.6)
elpa1_merge_systems_real_template.F90	accept


#official3.9.6 is not changed from official3.9.2
Occupation_Number_LDA_U.c

#patch3.9.6.tar.gz exists in patch3.9.6, this is probably miss or bug of app to expand.
