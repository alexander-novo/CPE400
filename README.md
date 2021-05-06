# UNR CPE 400 Networking Project

## Building the assignment
A Makefile is provided in the root directory which makes all executables, and the report with the default target. If you want to just build the executables, the `exec` target can be used (`make exec`). The built executable can be found in the root directory (`proj`) and can be run with no arguments or the `-h` switch to see a help menu which will explain all of the available options and how to run the program.

The only prerequisites for the `exec` target are a working `g++` on the path which supports C++ 17.

## Running Executables

```
Usage: proj ospf <config> [options]                                          (1)
Usage: proj ttl  <config> [options]                                          (2)
   or: proj -h                                                               (3)

(1) Run the project with the given configuration file using a standard OSPF
    routing algorithm.
(2) Run the project with the given configuration file using the proposed TTL
    routing algorithm.
(3) Print this help menu.

OPTIONS
  -s   <seed>  Set the seed used for generating random observations.
  -v           Set the output to verbose. Prints initialization information
               and final network state.
```

## Building the report
The report is built in the default target, so running `make` will generate the report. This requires a TeX distribution - I recommend [TeX Live](https://www.tug.org/texlive/) on *nix or [MiKTeX](https://miktex.org/) on Windows. The generated report can be found in the `Report/` folder.