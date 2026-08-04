# Polarimetric phasing

Solving for each antenna's complex gains using full IQU calibrator source model

With LAPACK and `-O2`, it takes about 10 minutes to solve 2K channels with 22 antennas.

It has to solve in <2mins.

One way is to do multi-threading, but that adds computation load.

## Math

Given `N` antennas, we are solving for `4N` parameters that produce `8 (N C 2)` data points for every channel.
Unsurprisingly, it is slow.


# This packages uses the following softwares/packages that are bundled together:


| package | usage |
|---------|-------|
| lbfgs | Limited memory BFGS optimization code |
| cminpack | Use the Levenburg Marquardt solver which solves for gains |
| lute  | Provide an interface to the GMRT LTA file |
| fmt   | To write in specific format the solved gains |
| lapack,blas | used by cminpack to accelrate solving |
| cxxopts | To provide simple CLI |

## removing strictly not necessary dependencies

this code needs to run on a really old system with gcc 4.8.5. 

so i am removing `fmt` and `cxxopts` dependencies
basically all modern ones

