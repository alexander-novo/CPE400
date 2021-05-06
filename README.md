# UNR CS 479 Programming Assignment 1

## Building the assignment
A Makefile is provided in the root directory which makes all executables, plots, and the report with the default target. If you want to just build the executables, the `exec` target can be used (`make exec`). Built executables can be found in their corresponding source folders (`Q1-Sampling/`, `Q2-Quantization/`, etc.) and each one can be run with no arguments or the `-h` switch to see a help menu which will explain all of the available options and how to run the program.

The only prerequisites for the `exec` target are a working `g++` on the path which supports C++ 14 and OpenMP.

## Running Executables

### Parts 1 & 2 - Bayes Classifier
```
Usage: classify-bayes <data set> [options]                            (1)
   or: classify-bayes -h                                              (2)

(1) Run a Bayes classifier on a specific data set.
    Data sets available are 'A' and 'B'.
(2) Print this help menu.

OPTIONS
  -s   <seed>  Set the seed used to generate samples.
               Defaults to 1.
  -psN <file>  Print all observations from sample N to a file.
               N can be 1 to 2.
  -pmN <file>  Print all misclassified observations from sample N to a file.
               N can be 1 to 2.
  -pdf <file>  Print a graph of the probability density function to a file.
               There will be an extra column which shows which class is more
               likely at that point. Will also allow correct calculation of
               zmax in the -pdb file.
  -peb <file>  Print a graph of the error bound function and its derivative
               to a file.
  -pdb <file>  Print the parameters of the decision boundary, along with
               other miscellaneous parameters needed for plotting.
  -c   <case>  Override which discriminant case is to be used.
               <case> can be 1-3. Higher numbers are more computationally
               expensive, but are correct more of the time.
               By default, the case will be chosen automatically.
```

### Parts 3 & 4 - Euclidean Classifier
```
Usage: classify-euclid <data set> [options]                            (1)
   or: classify-euclid -h                                              (2)

(1) Run a Bayes classifier on a specific data set.
    Data sets available are 'A' and 'B'.
(2) Print this help menu.

OPTIONS
  -s   <seed>  Set the seed used to generate samples.
               Defaults to 1.
  -psN <file>  Print all observations from sample N to a file.
               N can be 1 to 2.
  -pmN <file>  Print all misclassified observations from sample N to a file.
               N can be 1 to 2.
  -pdb <file>  Print the parameters of the decision boundary.
```

## Building the report
The report is built in the default target, so running `make` will generate the report. This requires a TeX distribution - I recommend [TeX Live](https://www.tug.org/texlive/) on *nix or [MiKTeX](https://miktex.org/) on Windows. The generated report can be found in the `Report/` folder. [Gnuplot](http://www.gnuplot.info/) is used to generate all plots and [`pnmtopng`](http://netpbm.sourceforge.net/doc/pnmtopng.html) is used to convert `.pgm`/`.ppm` images (unreadable by TeX) into `.png` images. These can both be installed on Debian/Ubuntu by using `sudo apt install gnuplot netpbm` and can be found in prebuilt binary form for Windows on their websites (see earlier links).