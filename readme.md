Galaxy wrapper for SNAP
========================

This wrapper is copyright 2016-2017 by Yating Liu

This is a wrapper for the gene prediction tool SNAP. SNAP is a general purpose gene finding program suitable for both eukaryotic and prokaryotic genomes. SNAP is an acroynm for Semi-HMM-based Nucleic Acid
Parser.

Reference
----------------------

    Korf I. Gene finding in novel Genomes. BMC Bioinformatics 2004, 5:59

Installation
-----------------------

To install SNAP, please download SNAP from

http://korflab.ucdavis.edu/Software/snap-2013-11-29.tar.gz

and follow the installation instractions. The software is routinely compiled and tested on Mac OS X. It should compile
fine on any Linux/Unix type operating systems.
The default compiler is gcc. If you have gcc installed, the easiest is to just compile as:
```
  make
```

The ZOE environment variable is used by SNAP to find the HMM files. Set this
to the directory containing this file. For example, if you unpackaged the tar-ball in /usr/local/snap, set the ZOE environment variable to /usr/local/snap

```
    setenv ZOE /usr/local/snap # csh, tcsh, etc
```
  or
```
    export ZOE=/usr/local/snap # sh, bash, etc
```
To install the wrapper copy the snap folder in the galaxy tools and modify the $GALAXY_ROOT/config/tool_conf.xml file to make the tool available to Galaxy. For example:
```
<tool file="galaxy/tools/myTools/snap/snap.xml" />
```





