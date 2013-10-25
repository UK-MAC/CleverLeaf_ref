# CleverLeaf

CleverLeaf is a hydrodynamics mini-app that extends
[CloverLeaf](http://warwick-pcav.github.io/CloverLeaf/) with Adaptive Mesh
Refinement using the
[SAMRAI](http://computation.llnl.gov/casc/SAMRAI/index.html) toolkit from
Lawrence Livermore National Laboratory. The primary goal of CleverLeaf is to
evaluate the application of AMR to the Lagriangian-Eulerian hydrodynamics scheme
used by CloverLeaf.

## Quick Start

Ensure you have SAMRAI (and it's dependencies installed). Edit `Makefile`,
ensuring that `HDF_DIR` and `SAMRAI_DIR` point to the correct locations.

To build the reference version, type: `make ref`

To run the test problem, type: `make test`  
**N.B.** Running the test assumes that the command `mpirun` can be used to
launch an executable.

## Detailed Build Instructions

CleverLeaf's only explicit dependency is:

- SAMRAI (v. 3.6.3)

Other libraries, such as HDF5 and Boost, will need to be installed in order to
build SAMRAI.

CleverLeaf uses a number of variables in the Makefile to control complation.
These variables can be set as required to get CleverLeaf to compile on your
system:

- `CXX`: The C++ compiler to use.
- `F90`: The Fortran90 compiler to use.
- `HDF_DIR`: The location of HDF5.
- `BOOST_DIR`: The location of the Boost library (headers only).
- `SAMRAI_DIR`: The location of the SAMRAI library.

Once these variables have been set, CleverLeaf should build correctly. However,
the `CPPFLAGS`, `FFLAGS`, and `LDFLAGS` can be modified if necessary.

CleverLeaf currently supports to targets: `ref` and `openmp`. The `ref` target
builds the vanilla version of CleverLeaf using MPI only. The `openmp` target
builds a hybrid version that uses both MPI and OpenMP.

## About

Like CloverLeaf, CleverLeaf solves Euler's equations on a staggered grid,
using an explicit, second-order accurate method. Each cell stores three
values: energy, density, and pressure. A velocity vector is stored at each
cell corner. The equations are solved in two dimensions. The following
description of the equations is from the CloverLeaf documentation:

> The compressible Euler equations are a set of three partial
> differential equations that describe the conservation of energy, mass and
> momentum in a system. CloverLeaf produces a second-order accurate solution
> using explicit finite volume methods. It first performs a Lagrangian step,
> using a predictor-corrector scheme to advance the solution forward by a
> calculated time delta. This step causes the mesh to move with fluid
> velocity, so an advective remap is used in order to return the mesh to its
> original state. A second-order Van Leer scheme is used, with the advective
> sweep being performed in the x and y directions for the energy, mass and
> momentum. The initial sweep direction alternates between steps, providing
> second order accuracy. The flow direction mush be calculated during the
> remap to allow data from the "upwind" direction to be used. Although the
> deformation of the grid does not actually move cell vertices, the average
> velocity on a cell face is used to approximate a flux through each face for
> the advection of material.
> 
> The compressible Euler equations form a hyperbolic system and therefore
> generate discontinuities in the form of shock waves. The second-order
> approximation will fail at these discontinuities and cause "ringing" in the
> solution. To avoid this, an artificial viscous pressure is used, which makes
> the solution first order in the presence of shock waves.  This preserves
> monotonicity in the solution, by behaving as a simple addition to the
> pressure.
> 
> The timestep control uses the maximum sound speed as an upper bound for the
> time delta. The timestep is thus limited to the time it would take for the
> the highest speed sound wave to cross a cell. The timestep is then
> multiplied by a safety factor to preserve the stability of the solution. The
> timestep control contains two further tests: one to ensure that a vertex
> can't overtake another as the mesh deforms, and one to ensure that a cell
> cannot deform such that it's volume becomes negative.
> 
> In order to close the system of equations, we use an equation of state,
> which calculates the pressure and sound speed in a cell, given its energy
> and density.  CloverLeaf uses the ideal gas equation of state with a gamma
> value of 1.4. 
> 
> Currently, CloverLeaf only solves for a single material, although multiple
> states (pressures, densities, and velocities) of this material can exist in
> the problem domain. Support for multiple materials will be added into a
> future release.

The code is split into three main classes: LagrangianEulerianIntegrator, which
is responsible for advancing the solution across the entire adaptive grid;
LagrangianEulerianLevelIntegrator, which is responsible for advancing the
solution on a single level; and Cleverleaf, which implements the
LagrangianEulerianPatchStrategy interface, providing the methods needed to
advance the solution on a single patch.

The physics in CleverLeaf uses the Fortran kernel functions written for
CloverLeaf, and applies them to each patch in the adaptive hierarchy.  Control
code has been added which adapts the grid and synchronises the different levels
of refinement. To coarsen the field variables, two additional operators have
been added: CartesianCellDoubleMassWeightedAverage, and
CartesianCellDoubleVolumeWeightedAverage. These are used to conservatively
coarsen energy and pressure.
