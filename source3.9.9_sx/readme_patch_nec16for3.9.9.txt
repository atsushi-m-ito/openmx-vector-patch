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
