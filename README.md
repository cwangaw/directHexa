# directHexa
This is a demonstration code for the paper Direct Serendipity Finite Elements on Cuboidal Hexahedra.

## Authors
- [Todd Arbogast](arbogast@oden.utexas.edu)
- [Chuning Wang](cwangaw@utexas.edu)


## Build and Run
You should make the files to build an executable file directpoly.$(TARG), where TARG is a variable containing the instruction set architecture of your processor. Note that you need CBLAS, LAPACKE, and UMFPACK library to successfully build the file. Before running, you can change the input parameters in `infile`.

```bash
make
./directpoly.$(TARG)
```

## Direct Serendipity Space
We choose vertex basis functions to be linear on each edge, edge basis functions to be chebyshev polynomials times a bubble function on each edge, face basis functions to be nodal basis functions, and cell basis functions to be monomials with regard to $`\lambda_x, \lambda_y, \lambda_z`$ times a bubble functions. They are constructed in `directSerendipityFE.cpp` and assembled in `ellipticPDE.cpp` to serve as global basis functions.

You can choose the supplemental functions in `infile` from smooth functions, piecewise polynomials on marching tetrahedra, and piecewise polynomials on diamond lattice cells.

## Solve a PDE
Problem formulation is given in the heading comments of `main.cpp`. Coefficients a, b, c, D, and all the other related data could be modified in `fcns.cpp`. Note that the source function as well as Riemann boundary condition are defaultly calculated by true solution. If in your formulation, analytical solution is unknown, please rewrite these parts.

## Output
Basic results of running the codes are printed to the terminal. All the output files would be stored preambly in `test/` directory, which could be modified in the first line of `infile`.

## Acknowledgments
The development of this code has been supported by the U.S. National Science Foundation.

## Copyrights
Copyright (C) 2023 Todd Arbogast and Chuning Wang

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the [GNU General Public License](https://www.gnu.org/licenses/) for more details.