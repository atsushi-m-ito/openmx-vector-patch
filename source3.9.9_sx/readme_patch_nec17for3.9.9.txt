#Patch for Vector Processors, vec17for3.9.9
1) fix the miss of merge
Fixed files:
openmx.c
MD_pac.c


Note that, by the modification of MD_pac.c, 
runtest shows not small difference of Force in the cases of CO and H2O as follows:
----------------------------------------------------------------------------------------------------------------------
$cat runtest.result
   1  input_example/Benzene.dat        Elapsed time(s)=   13.53  diff Utot= 0.000000000005  diff Force= 0.000000000009
   2  input_example/C60.dat            Elapsed time(s)=   32.27  diff Utot= 0.000000000021  diff Force= 0.000000000002
   3  input_example/CO.dat             Elapsed time(s)=   38.20  diff Utot= 0.000000000010  diff Force= 0.002236432722
   4  input_example/Cr2.dat            Elapsed time(s)=   18.72  diff Utot= 0.000000000434  diff Force= 0.000000000137
   5  input_example/Crys-MnO.dat       Elapsed time(s)=   29.50  diff Utot= 0.000000000001  diff Force= 0.000000000034
   6  input_example/GaAs.dat           Elapsed time(s)=   54.21  diff Utot= 0.000000000013  diff Force= 0.000000000000
   7  input_example/Glycine.dat        Elapsed time(s)=   16.14  diff Utot= 0.000000000002  diff Force= 0.000000000003
   8  input_example/Graphite4.dat      Elapsed time(s)=   10.49  diff Utot= 0.000000000020  diff Force= 0.000000000013
   9  input_example/H2O-EF.dat         Elapsed time(s)=   13.29  diff Utot= 0.000000000002  diff Force= 0.000000001718
  10  input_example/H2O.dat            Elapsed time(s)=   11.27  diff Utot= 0.000000000002  diff Force= 0.000741132742
  11  input_example/HMn.dat            Elapsed time(s)=   27.99  diff Utot= 0.000000000457  diff Force= 0.000000000002
  12  input_example/Methane.dat        Elapsed time(s)=   10.69  diff Utot= 0.000000000004  diff Force= 0.000000000000
  13  input_example/Mol_MnO.dat        Elapsed time(s)=   20.34  diff Utot= 0.000000000099  diff Force= 0.000000000013
  14  input_example/Ndia2.dat          Elapsed time(s)=   11.34  diff Utot= 0.000000000003  diff Force= 0.000000000001


Total elapsed time (s)      307.99
----------------------------------------------------------------------------------------------------------------------
However, this is caused by the change of algorism to accelarate the structure relaxation which is "MD.type == opt".
These differences of forces are not error.
