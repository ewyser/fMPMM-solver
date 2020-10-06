# fMPMM-solver

## Fast and efficient MATLAB-based MPM solver in an explicit formulation.

Copyright (C) 2020  Emmanuel Wyser, Yury Alkhimenkov, Michel Jaboyedoff, Yury Y. Podladchikov.

fMPMM-solver is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

fMPMM-solver is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with FMPMM-solver.  If not, see <http://www.gnu.org/licenses/>.

### Current version of FMPMM-solver can be found at:

https://bitbucket.org/ewyser/fmpmm-solver/

### fMPMM-solver was released in:

ADD REFERENCE

### Distributed software, directory content:

v1.0 | numerical investigation of the compression of an elastic column, the cantilever beam problem and the elastoplastic slump problem to reproduce the results from the manuscript
v1.1 | numerical investigation of the compression of an elastic column, the cantilever beam problem and the elastoplastic slump problem to reproduce the results from the manuscript

Changes: October, 10th 2020
- more efficient vectorisation of matrix multiplication
- CPDI/CPDi2q shape functions do not require the usage of sparse matrix
- corrected implementation of rotated stress

### QUICK START:

Open Matlab, select the example routines, run the script and post-process results with the postprocessing.m routine

### COMPATIBILITY:

fMPMM-solver was designed under a MATLAB architecture. However, it can also be run with OCTAVE. We did not check the compatibility of the solver with every version of MATLAB nor OCTAVE, but we provide a non-exhaustive list of compatibility

MATLAB version: R2018b, R2018a, R2017b, R2016a, R2013b

OCTAVE version: 5.1.0.0

### Contact: manuwyser@gmail.com

