# Polarimetric phasing

Solving for each antenna's complex gains using full IQU calibrator source model

With LAPACK and `-O2`, it takes about 10 minutes to solve 2K channels with 22 antennas.

It has to solve in <2mins.

One way is to do multi-threading, but that adds computation load.


# This packages uses the following softwares/packages that are bundled together:


| package | usage |
|---------|-------|
| cminpack | Use the Levenburg Marquardt solver which solves for gains |
| lute  | Provide an interface to the GMRT LTA file |
| fmt   | To write in specific format the solved gains |
| lapack,blas | used by cminpack to accelrate solving |



