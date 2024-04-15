# lammps-tutorials

## Introduction & Set Up

LAMMPS (**L**arge-scale **A**tomic/**M**olecular **M**assively **P**arallel **S**imulator) is a classical molecular dynamics simulation (MD) program that models ensembles of particles in a liquid, solid, or gaseous state. It can model atomic, polymeric, biological, solid-state, granular, coarse-grained, or macroscopic systems using a variety of interatomic potentials (force fields) and boundary conditions. It can model 2d or 3d systems with sizes ranging from only a few particles up to billions.

In the most general sense, LAMMPS integrates Newton’s equations of motion for a collection of interacting particles. A single particle can be treated as atomistic (i.e. an atom, molecule or electron) or course-grained (a cluster of atoms, or a mesoscopic or macroscopic clump of material). The interaction models that LAMMPS includes are mostly short-ranged in nature; some long-range models are included as well.

Tutorials require [LAMMPS MD software package](https://github.com/lammps/lammps) to be installed and built on your local machine. Detailed [documentation](https://docs.lammps.org/) is available, however a simplified step by step guide is given below. This guide is for a build on Linux, specifically Ubuntu 23.10.

### LAMMPS & Required Package Installation
1. Clone LAMMPS Github repo either using [SSH protocol](https://docs.github.com/en/authentication/connecting-to-github-with-ssh) (`git clone git@github.com:lammps/lammps.git`), or using HTTPS (`git clone https://github.com/lammps/lammps.git`)
2. Set up virtual environment e.g. using Conda (`conda activate <virtual_environment>`)
3. Install required packages (further packages will be added as more exercises are added to the tutorial):
   * `sudo apt install gcc g++ gfortran wget make lammps lmpi_cxx libopenmpi-dev mpi-default-bin mpi-default-dev libfftw3-dev libjpeg-dev libpng12-dev libreadline-dev`
5. Certain auxiliary tools require the `ifort` command which requires installation of Intel Fortran Compiler. This can be achieved either as part of the [HPC Toolkit](https://www.intel.com/content/www/us/en/developer/tools/oneapi/hpc-toolkit-download.html?operatingsystem=linux&distributions=aptpackagemanager), or as a [standalone release](https://www.intel.com/content/www/us/en/developer/articles/tool/oneapi-standalone-components.html#fortran)
   * Once installed, prior to running the compiler you need to set certain environment variables via command `source /opt/intel/oneapi/setvars.sh intel64`. This will create a symlink so that the `ifort` command can be used in `lammps\tools\Makefile`
   * May be specific to my system, but I also found that I needed to change the `chain` and `micelle2d` rules in `lammps\tools\Makefile` to correctly compile `chain.f90` and `micelle2d.f90` into `chain.o` and `micelle2d.o` respectively
   ```
   chain:	chain.f90
	  ifort -c chain.f90 -o chain
   micelle2d:	micelle2d.f90
	  ifort -c micelle2d.f90 -o micelle2d
   ```

### LAMMPS Build
1. All commands must be run in LAMMPS `src` directory (`cd lammps/src`)
2. Build serial LAMMPS executable using GNU g++ (`make serial`)
   * Command `make serial` creates `lmp_serial` binary. It does not require MPI (Message Passing Interface) and is intended for running simulations on a single processor without parallelization
4. Build parallel LAMMPS executable with MPI (`make mpi`)
   * Command `make mpi` creates `lmp_mpi` binary. It is designed for parallel execution using MPI and allows LAMMPS to run simulations across multiple processors. This can significantly improve performance for large-scale simulations
5. Build LAMMPS executable and library (`make ubuntu`)
   * This generates the LAMMPS executable in the 'static' mode. If you want to generate it in the 'shared' mode, you need to run `make mode=shared ubuntu`
6. Navigate to `tools` directory (`cd lammps/tools`) and build LAMMPS tools (`make all`)
  
