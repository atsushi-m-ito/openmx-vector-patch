# OpenMX Patch for Vector Processors

## What is OpenMX

The OpenMX[^1] is a numerical simulation code based on density functional theory, and is developed by Dr. Taisuke OZAKI and OpenMX developpers.

## What is OpenMX Patch for Vector Processors

This patch is a set of codes tuned for speeding up computation on vector processors.
Here, the vector processor is assumed to be NEC SX-Aurora TSUBASA, and is a processor with a very long vector length (256 or more). Other vector processors and SIMDs with long vector lengths, such as the Intel AVX512, may also contain valid tunings, but are currently unverified.

## Target version

This patch for vector processors(VP) will be distributed according to the original version of OpenMX. Please check the following correspondence table.

|Patch for VP ver.| Original OpenMX ver.|
|---|---|
| nec14 | 3.9.2 |
| nec13 | 3.9.2 |
| none  | 3.9.1 |
| none  | 3.9.0 | 

## How to install

This section describes how to install the original OpenMX and how to apply this patch.

1. Download the original OpenMX, the original patch, and this patch for VP. And, put them in your directory.

```
$ ls 
openmx3.9.tar.gz   patch3.9.2.tar.gz   patch_nec12_for392.tar.gz
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
$ cp -rp source source392nec13
$ cd source392nec13
$ tar xvfz ../../patch3.9.2.tar.gz
```

4. Apply the patch for VP

```
$ tar xvfz ../../patch_nec12_for392.tar.gz
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


# Licenses

Both the original OpenMX code and the patch for VP are licensed as GPL3.




[^1]: http://www.openmx-square.org/