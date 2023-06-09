# OpenMX Patch for Vector Processors

## What is OpenMX

The OpenMX[^1] is a numerical simulation code based on density functional theory, and is developed by Dr. Taisuke OZAKI and OpenMX developpers.

## What is OpenMX Patch for Vector Processors

This patch is a set of codes tuned for speeding up computation on vector processors.
Here, the vector processor is assumed to be NEC SX-Aurora TSUBASA, and is a processor with a very long vector length (256 or more). Other vector processors and SIMDs with long vector lengths, such as the Intel AVX512, may also contain valid tunings, but are currently unverified.

## Target version

This patch for vector processors(VP) will be distributed according to the original version of OpenMX. Please check the following correspondence table.

|Patch for VP ver.| Original OpenMX ver.|note|
|---|---|---|
| vec17 | 3.9.9 ||
| vec16 | 3.9.9 |deprecated|
| vec15 | 3.9.6 ||
| vec14 | 3.9.2 ||
| vec13 | 3.9.2 ||
| none  | 3.9.1 ||
| none  | 3.9.0 ||

## How to install

This section describes how to install the original OpenMX and how to apply this patch.

1. Download the original OpenMX, the original patch, and this patch for VP. And, put them in your directory.

```
$ ls 
openmx3.9.tar.gz   patch3.9.9tar.gz   openmx-vector-patch-vec17for3.9.9.tar.gz
```

2. Expand the original OpenMX

```
$ tar xvfz openmx3.9.tar.gz
$ cd openmx3.9
$ ls
DFT_DATA19   source   work
```

3. Apply the original patch

```
$ cp -rp source source399vec17
$ cd source399vec17
$ tar xvfz ../../patch3.9.9.tar.gz
```

4. Apply the patch for VP

```
$ tar xvfz ../../openmx-vector-patch-vec17for3.9.9.tar.gz
$ cp -rp openmx-vector-patch-vec17for3.9.9/source3.9.9_sx/* ./
$ rm -r openmx-vector-patch-vec17for3.9.9/
```

5. Check and modify the makefile according to your environment. And, Make it.

```
$ vi makefile
$ make -j install
```

6. Check culculation. However, follow the rules of your computer environment when submitting jobs.

```
$ cd ../work

# Before benchmark, check your computer rule for job submission.  
$ ./openmx -runtest
$ ./openmx -runtestL
```
### Note

The result of the runtest, which is "runtest.result", may shows not small differences of forces for cases of CO and H2O as follows:

```
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
```

This is caused by the change of algorism in "MD_pac.c" to accelarate the structure relaxation which is "MD.type == opt".
These differences of forces are not error.




## Reference

The contents of tuning are summarized in the following paper. When using this code, it would be appreciated if you can cite this paper.

[Atsushi M. Ito, Arimichi Takayama, Osamu Watanabe, Vijendra Singh, Shubham Tyagi, and Shashank S. Singh, "Tuning of Density Functional Theory Simulation on Vector Processor System - Plasma Simulator Raijin -" Plasma and Fusion Research, Rapid Communications, 15 (2020) 1203085.](http://www.jspf.or.jp/PFR/PFR_articles/pfr2020/pfr2020_15-1203085.html)

## Licenses

Both the original OpenMX code and the patch for VP are licensed as GPL3.


[^1]: http://www.openmx-square.org/
