# LAMMPS Tutorials

## Contents
1. [LAMMPS Introduction & Set Up](https://github.com/c-vandenberg/lammps-tutorials/blob/master/README.md#1-lammps-introduction--set-up)<br>
	1.1 [LAMMPS Installation & Required Package Installation](https://github.com/c-vandenberg/lammps-tutorials/blob/master/README.md#11-lammps-installation--required-package-installation)<br>
	1.2 [LAMMPS Build with Make](https://github.com/c-vandenberg/lammps-tutorials/blob/master/README.md#12-lammps-build-with-make)<br>
	1.3 [Running CMake Build LAMMPS](https://github.com/c-vandenberg/lammps-tutorials/blob/master/README.md#13-running-cmake-build-lammps)<br>
	1.4 [Configuring CLion Debugger with LAMMPS](https://github.com/c-vandenberg/lammps-tutorials/blob/master/README.md#14-configuring-clion-debugger-with-lammps)<br>
2. [Simon Gravelle Tutorial Level 1: Lennard-Jones Fluid - The Very Basics of LAMMPS](https://github.com/c-vandenberg/lammps-tutorials/blob/master/simon-gravelle/level-1/2-lennard-jones-fluid/README.md#2-simon-gravelle-tutorial-level-1---lennard-jones-fluid-the-very-basics-of-lammps)<br>
	2.1 [Lennard-Jones Fluid `first-input.lammps` Script](https://github.com/c-vandenberg/lammps-tutorials/tree/master/simon-gravelle/level-1/2-lennard-jones-fluid/2.1-first-input#21-lennard-jones-fluid-first-inputlammps-script)<br>
  	&nbsp; &nbsp; 2.1.1 [Exercise](https://github.com/c-vandenberg/lammps-tutorials/tree/master/simon-gravelle/level-1/2-lennard-jones-fluid/2.1-first-input#211-exercise)<br>
  	&nbsp; &nbsp; 2.1.2 [Data Analysis](https://github.com/c-vandenberg/lammps-tutorials/tree/master/simon-gravelle/level-1/2-lennard-jones-fluid/2.1-first-input#212-data-analysis)<br>
  	&nbsp; &nbsp; 2.1.3 [Input Script Command Syntax](https://github.com/c-vandenberg/lammps-tutorials/tree/master/simon-gravelle/level-1/2-lennard-jones-fluid/2.1-first-input#213-input-script-command-syntax)<br>
	2.2 [Lennard-Jones Fluid `improved-input.min.lammps` & `improved-input.min.lammps` Scripts](https://github.com/c-vandenberg/lammps-tutorials/tree/master/simon-gravelle/level-1/2-lennard-jones-fluid/2.2-improved-input#22-lennard-jones-fluid-improved-inputminlammps--improved-inputminlammps-scripts)<br>
	&nbsp; &nbsp; 2.2.1 [Lennard-Jones Fluid `improved_input.min.lammps` Script](https://github.com/c-vandenberg/lammps-tutorials/tree/master/simon-gravelle/level-1/2-lennard-jones-fluid/2.2-improved-input/2.1.1-improved-min-input#221-lennard-jones-fluid-improved_inputminlammps-script)<br>
	&nbsp; &nbsp; 2.2.2 [Lennard-Jones Fluid `improved_input.md.lammps` Script](https://github.com/c-vandenberg/lammps-tutorials/tree/master/simon-gravelle/level-1/2-lennard-jones-fluid/2.2-improved-input/2.1.2-improved-md-input#222-lennard-jones-fluid-improved_inputmdlammps-script)<br>
	2.3 [Further Exercises: Lennard-Jones Fluid Further Exercises: Lost Atoms Error](https://github.com/c-vandenberg/lammps-tutorials/tree/master/simon-gravelle/level-1/2-lennard-jones-fluid/further-exercises/2.3-lost-atoms-error#23-further-exercises-lennard-jones-fluid-further-exercises-lost-atoms-error)<br>
	&nbsp; &nbsp; 2.3.1 [Problem](https://github.com/c-vandenberg/lammps-tutorials/tree/master/simon-gravelle/level-1/2-lennard-jones-fluid/further-exercises/2.3-lost-atoms-error#231-problem)<br>
	&nbsp; &nbsp; 2.3.2 [Solution](https://github.com/c-vandenberg/lammps-tutorials/tree/master/simon-gravelle/level-1/2-lennard-jones-fluid/further-exercises/2.3-lost-atoms-error#232-solution)<br>
	2.4 [Further Exercises: Create a Demixed Dense Phase](https://github.com/c-vandenberg/lammps-tutorials/tree/master/simon-gravelle/level-1/2-lennard-jones-fluid/further-exercises/2.4-demixed-dense-phase#24-further-exercises-create-a-demixed-dense-phase)<br>
	&nbsp; &nbsp; 2.4.1 [Problem](https://github.com/c-vandenberg/lammps-tutorials/tree/master/simon-gravelle/level-1/2-lennard-jones-fluid/further-exercises/2.4-demixed-dense-phase#241-problem)<br>
	&nbsp; &nbsp; 2.4.2 [Solution](https://github.com/c-vandenberg/lammps-tutorials/tree/master/simon-gravelle/level-1/2-lennard-jones-fluid/further-exercises/2.4-demixed-dense-phase#242-solution)<br>
	2.5 [Further Exercises: From Atoms to Molecules](https://github.com/c-vandenberg/lammps-tutorials/tree/master/simon-gravelle/level-1/2-lennard-jones-fluid/further-exercises/2.5-from-atoms-to-molecules#25-further-exercises-from-atoms-to-molecules)<br>
	&nbsp; &nbsp; 2.5.1 [Problem](https://github.com/c-vandenberg/lammps-tutorials/tree/master/simon-gravelle/level-1/2-lennard-jones-fluid/further-exercises/2.5-from-atoms-to-molecules#251-problem)<br>
	&nbsp; &nbsp; 2.5.2 [Solution](https://github.com/c-vandenberg/lammps-tutorials/tree/master/simon-gravelle/level-1/2-lennard-jones-fluid/further-exercises/2.5-from-atoms-to-molecules#252-solution)<br>
3. [Simon Gravelle Tutorial Level 1: Carbon Nanotube Deformation](https://github.com/c-vandenberg/lammps-tutorials/tree/master/simon-gravelle/level-1/3-carbon-nanotube-deformation#3-simon-gravelle-tutorial-level-1-carbon-nanotube-deformation)<br>
	3.1 [Deformation of Carbon Nanotube with Unbreakable Bonds](https://github.com/c-vandenberg/lammps-tutorials/tree/master/simon-gravelle/level-1/3-carbon-nanotube-deformation/3.1-cnt-unbreakable-bonds#31-deformation-of-carbon-nanotube-with-unbreakable-bonds)<br>
 	&nbsp; &nbsp; 3.1.1 [Exercise](https://github.com/c-vandenberg/lammps-tutorials/tree/master/simon-gravelle/level-1/3-carbon-nanotube-deformation/3.1-cnt-unbreakable-bonds#311-exercise)<br>
  	&nbsp; &nbsp; 3.1.2 [Introduction](https://github.com/c-vandenberg/lammps-tutorials/tree/master/simon-gravelle/level-1/3-carbon-nanotube-deformation/3.1-cnt-unbreakable-bonds#312-introduction)<br>
	&nbsp; &nbsp; 3.1.3 [Data Analysis](https://github.com/c-vandenberg/lammps-tutorials/tree/master/simon-gravelle/level-1/3-carbon-nanotube-deformation/3.1-cnt-unbreakable-bonds#313-data-analysis)<br>
	&nbsp; &nbsp; 3.1.4 [Input Script Command Syntax](https://github.com/c-vandenberg/lammps-tutorials/tree/master/simon-gravelle/level-1/3-carbon-nanotube-deformation/3.1-cnt-unbreakable-bonds#314-input-script-command-syntax)<br>
	3.2 [Deformation of Carbon Nanotube with Breakable Bonds](https://github.com/c-vandenberg/lammps-tutorials/tree/master/simon-gravelle/level-1/3-carbon-nanotube-deformation/3.2-cnt-breakable-bonds#32-deformation-of-carbon-nanotube-with-breakable-bonds)<br>
	&nbsp; &nbsp; 3.2.1 [Exercise](https://github.com/c-vandenberg/lammps-tutorials/tree/master/simon-gravelle/level-1/3-carbon-nanotube-deformation/3.2-cnt-breakable-bonds#321-exercise)<br>
	&nbsp; &nbsp; 3.2.2 [Introduction](https://github.com/c-vandenberg/lammps-tutorials/tree/master/simon-gravelle/level-1/3-carbon-nanotube-deformation/3.2-cnt-breakable-bonds#322-introduction)<br>
	&nbsp; &nbsp; 3.2.3 [Differences in Topology File](https://github.com/c-vandenberg/lammps-tutorials/tree/master/simon-gravelle/level-1/3-carbon-nanotube-deformation/3.2-cnt-breakable-bonds#323-differences-in-topology-file)<br>
	&nbsp; &nbsp; 3.2.4 [Data Analysis](https://github.com/c-vandenberg/lammps-tutorials/tree/master/simon-gravelle/level-1/3-carbon-nanotube-deformation/3.2-cnt-breakable-bonds#324-data-analysis)<br>
	&nbsp; &nbsp; 3.2.5 [Input Script Command Syntax](https://github.com/c-vandenberg/lammps-tutorials/tree/master/simon-gravelle/level-1/3-carbon-nanotube-deformation/3.2-cnt-breakable-bonds#325-input-script-command-syntax)<br>
	3.3 [Further Exercises: Plot the Carbon Nanotube Stress-Strain Curves](https://github.com/c-vandenberg/lammps-tutorials/tree/master/simon-gravelle/level-1/3-carbon-nanotube-deformation/further-exercises/3.3-stress-strain-curve#33-further-exercises-plot-the-carbon-nanotube-stress-strain-curves)<br>
	&nbsp; &nbsp; 3.3.1 [Problem](https://github.com/c-vandenberg/lammps-tutorials/tree/master/simon-gravelle/level-1/3-carbon-nanotube-deformation/further-exercises/3.3-stress-strain-curve#331-problem)<br>
	&nbsp; &nbsp; 3.3.2 [Solution](https://github.com/c-vandenberg/lammps-tutorials/tree/master/simon-gravelle/level-1/3-carbon-nanotube-deformation/further-exercises/3.3-stress-strain-curve#332-solution)<br>
	&nbsp; &nbsp; 3.3.3 [Data Analysis](https://github.com/c-vandenberg/lammps-tutorials/tree/master/simon-gravelle/level-1/3-carbon-nanotube-deformation/further-exercises/3.3-stress-strain-curve#333-data-analysis)<br>
	3.4 [Further Exercises: Flying Ice Cube Artifact](https://github.com/c-vandenberg/lammps-tutorials/tree/master/simon-gravelle/level-1/3-carbon-nanotube-deformation/further-exercises/3.4-flying-ice-cube-artifact#34-further-exercises-flying-ice-cube-artifact)<br>
	&nbsp; &nbsp; 3.4.1 [Introduction](https://github.com/c-vandenberg/lammps-tutorials/tree/master/simon-gravelle/level-1/3-carbon-nanotube-deformation/further-exercises/3.4-flying-ice-cube-artifact#341-introduction)<br>
	&nbsp; &nbsp; 3.4.2 [Problem](https://github.com/c-vandenberg/lammps-tutorials/tree/master/simon-gravelle/level-1/3-carbon-nanotube-deformation/further-exercises/3.4-flying-ice-cube-artifact#342-problem)<br>
 	&nbsp; &nbsp; 3.4.3 [Solution](https://github.com/c-vandenberg/lammps-tutorials/tree/master/simon-gravelle/level-1/3-carbon-nanotube-deformation/further-exercises/3.4-flying-ice-cube-artifact#343-solution)<br>
	3.5 [Further Exercises: Inert Gas (Ar) in The Carbon Nanotube](https://github.com/c-vandenberg/lammps-tutorials/tree/master/simon-gravelle/level-1/3-carbon-nanotube-deformation/further-exercises/3.5-inert-gas-in-carbon-nanotube#35-further-exercises-inert-gas-ar-in-the-carbon-nanotube)<br>
	&nbsp; &nbsp; 3.5.1 [Problem](https://github.com/c-vandenberg/lammps-tutorials/tree/master/simon-gravelle/level-1/3-carbon-nanotube-deformation/further-exercises/3.5-inert-gas-in-carbon-nanotube#351-problem)<br>
	&nbsp; &nbsp; 3.5.2 [Solution](https://github.com/c-vandenberg/lammps-tutorials/tree/master/simon-gravelle/level-1/3-carbon-nanotube-deformation/further-exercises/3.5-inert-gas-in-carbon-nanotube#352-solution)<br>
	3.6 [Further Exercises: Carbon Nanotube Membrane](https://github.com/c-vandenberg/lammps-tutorials/tree/master/simon-gravelle/level-1/3-carbon-nanotube-deformation/further-exercises/3.6-carbon-nanotube-membrane#36-further-exercises-carbon-nanotube-membrane)<br>
	&nbsp; &nbsp; 3.6.1 [Problem](https://github.com/c-vandenberg/lammps-tutorials/tree/master/simon-gravelle/level-1/3-carbon-nanotube-deformation/further-exercises/3.6-carbon-nanotube-membrane#361-problem)<br>
	&nbsp; &nbsp; 3.6.2 [Solution](https://github.com/c-vandenberg/lammps-tutorials/tree/master/simon-gravelle/level-1/3-carbon-nanotube-deformation/further-exercises/3.6-carbon-nanotube-membrane#362-solution)<br>
4. [Simon Gravelle Tutorial Level 2: Polymer in Water](https://github.com/c-vandenberg/lammps-tutorials/tree/master/simon-gravelle/level-2/4-polymer-in-water#4-simon-gravelle-tutorial-level-2-polymer-in-water)<br>
	4.1 [Preparing The Water Reservoir](https://github.com/c-vandenberg/lammps-tutorials/tree/master/simon-gravelle/level-2/4-polymer-in-water/4.1-pure-H2O#41-preparing-the-water-reservoir)<br>
	&nbsp; &nbsp; 4.1.1 [Introduction](https://github.com/c-vandenberg/lammps-tutorials/tree/master/simon-gravelle/level-2/4-polymer-in-water/4.1-pure-H2O#411-introduction)<br>
	&nbsp; &nbsp; 4.1.2 [Data Analysis](https://github.com/c-vandenberg/lammps-tutorials/tree/master/simon-gravelle/level-2/4-polymer-in-water/4.1-pure-H2O#412-data-analysis)<br>
	&nbsp; &nbsp; 4.1.3 [Input Script Command Syntax](https://github.com/c-vandenberg/lammps-tutorials/tree/master/simon-gravelle/level-2/4-polymer-in-water/4.1-pure-H2O#413-input-script-command-syntax)<br>
	&nbsp; &nbsp; 4.1.4 [References](https://github.com/c-vandenberg/lammps-tutorials/tree/master/simon-gravelle/level-2/4-polymer-in-water/4.1-pure-H2O#414-references)<br>
	4.2 [Preparing The Single PEG Polymer](https://github.com/c-vandenberg/lammps-tutorials/tree/master/simon-gravelle/level-2/4-polymer-in-water/4.2-single-PEG#42-preparing-the-single-peg-polymer)<br>
	&nbsp; &nbsp; 4.2.1 [Introduction](https://github.com/c-vandenberg/lammps-tutorials/tree/master/simon-gravelle/level-2/4-polymer-in-water/4.2-single-PEG#421-introduction)<br>
 	&nbsp; &nbsp; 4.2.2 [References](https://github.com/c-vandenberg/lammps-tutorials/tree/master/simon-gravelle/level-2/4-polymer-in-water/4.2-single-PEG#422-references)<br>
	4.3 [Solvating the PEG in Water](https://github.com/c-vandenberg/lammps-tutorials/tree/master/simon-gravelle/level-2/4-polymer-in-water/4.3-solvated-PEG#43-solvating-the-peg-in-water)<br>
	&nbsp; &nbsp; 4.3.1 [Introduction](https://github.com/c-vandenberg/lammps-tutorials/tree/master/simon-gravelle/level-2/4-polymer-in-water/4.3-solvated-PEG#431-introduction)<br>
	&nbsp; &nbsp; 4.3.2 [Input Script Command Syntax](https://github.com/c-vandenberg/lammps-tutorials/tree/master/simon-gravelle/level-2/4-polymer-in-water/4.3-solvated-PEG#432-input-script-command-syntax)<br>
	4.4 [Deforming/Stretching the Water Solvated PEG Polymer Molecule](https://github.com/c-vandenberg/lammps-tutorials/tree/master/simon-gravelle/level-2/4-polymer-in-water/4.4-solvated-PEG-deformation#44-deformingstretching-the-water-solvated-peg-polymer-molecule)<br>
 	&nbsp; &nbsp; 4.4.1 [Introduction](https://github.com/c-vandenberg/lammps-tutorials/tree/master/simon-gravelle/level-2/4-polymer-in-water/4.4-solvated-PEG-deformation#441-introduction)<br>
	&nbsp; &nbsp; 4.4.2 [Data Analysis](https://github.com/c-vandenberg/lammps-tutorials/tree/master/simon-gravelle/level-2/4-polymer-in-water/4.4-solvated-PEG-deformation#442-data-analysis)<br>
	&nbsp; &nbsp; 4.4.3 [Input Script Command Syntax](https://github.com/c-vandenberg/lammps-tutorials/tree/master/simon-gravelle/level-2/4-polymer-in-water/4.4-solvated-PEG-deformation#443-input-script-command-syntax)<br>
 	4.5 [Further Exercises: Radial Distribution Function of Water Solvated PEG Polymer Molecule](https://github.com/c-vandenberg/lammps-tutorials/tree/master/simon-gravelle/level-2/4-polymer-in-water/further-exercises/4.5-radial-distribution-function#45-further-exercises-radial-distribution-function-of-water-solvated-peg-polymer-molecule)<br>
  	&nbsp; &nbsp; 4.5.1 [Exercise](https://github.com/c-vandenberg/lammps-tutorials/tree/master/simon-gravelle/level-2/4-polymer-in-water/further-exercises/4.5-radial-distribution-function#451-exercise)<br>
	&nbsp; &nbsp; 4.5.2 [Introduction](https://github.com/c-vandenberg/lammps-tutorials/tree/master/simon-gravelle/level-2/4-polymer-in-water/further-exercises/4.5-radial-distribution-function#452-introduction)<br>
	&nbsp; &nbsp; 4.5.3 [Data Analysis](https://github.com/c-vandenberg/lammps-tutorials/tree/master/simon-gravelle/level-2/4-polymer-in-water/further-exercises/4.5-radial-distribution-function#453-data-analysis)<br>
	&nbsp; &nbsp; 4.5.4 [Input Script Command Syntax](https://github.com/c-vandenberg/lammps-tutorials/tree/master/simon-gravelle/level-2/4-polymer-in-water/further-exercises/4.5-radial-distribution-function#454-input-script-command-syntax)<br>
  	4.6 [Further Exercises: Add NaCl to Water Solvated PEG Polymer Molecule](https://github.com/c-vandenberg/lammps-tutorials/tree/master/simon-gravelle/level-2/4-polymer-in-water/further-exercises/4.6-salinate-solvent#46-further-exercises-add-nacl-to-water-solvated-peg-polymer-molecule)<br>
	&nbsp; &nbsp; 4.6.1 [Introduction](https://github.com/c-vandenberg/lammps-tutorials/tree/master/simon-gravelle/level-2/4-polymer-in-water/further-exercises/4.6-salinate-solvent#461-introduction)<br>
   	4.7 [Further Exercises: Evaluate The Deformation of The PEG Polymer Molecule](https://github.com/c-vandenberg/lammps-tutorials/tree/master/simon-gravelle/level-2/4-polymer-in-water/further-exercises/4.7-PEG-deformation-evaluation#47-further-exercises-evaluate-the-deformation-of-the-peg-polymer-molecule)<br>
	&nbsp; &nbsp; 4.7.1 [Introduction](https://github.com/c-vandenberg/lammps-tutorials/tree/master/simon-gravelle/level-2/4-polymer-in-water/further-exercises/4.7-PEG-deformation-evaluation#471-introduction)<br>
	&nbsp; &nbsp; 4.7.2 [Data Analysis](https://github.com/c-vandenberg/lammps-tutorials/tree/master/simon-gravelle/level-2/4-polymer-in-water/further-exercises/4.7-PEG-deformation-evaluation#472-data-analysis)<br>
 5. [Simon Gravelle: MDAnalysis Tutorials](https://github.com/c-vandenberg/lammps-tutorials/tree/master/simon-gravelle/md-analysis-tutorial#5-simon-gravelle-mdanalysis-tutorials)<br>
	5.1 [MDAnalysis Tutorials - Polymer in Water](https://github.com/c-vandenberg/lammps-tutorials/tree/master/simon-gravelle/md-analysis-tutorial/5.1-polymer-in-water#51-mdanalysis-tutorials---polymer-in-water)<br>
	&nbsp; &nbsp; 5.1.2 [Extract Temporal Evolution of Hydrogen Type 4 Atom (First Atom in PEG Group)](https://github.com/c-vandenberg/lammps-tutorials/tree/master/simon-gravelle/md-analysis-tutorial/5.1-polymer-in-water#511-extract-temporal-evolution-of-hydrogen-type-4-atom-first-atom-in-peg-group)<br>
	5.2 [MDAnalysis Tutorials - Carbon Nanotube (CNT) Breakable Bonds](https://github.com/c-vandenberg/lammps-tutorials/tree/master/simon-gravelle/md-analysis-tutorial/5.2-carbon-nanotube-deformation/cnt-breakable-bonds#52-mdanalysis-tutorials---carbon-nanotube-cnt-breakable-bonds)<br>
	&nbsp; &nbsp; 5.2.1 [Evolution of CNT Average Bond Length & Bond Number as a Function of Time](https://github.com/c-vandenberg/lammps-tutorials/tree/master/simon-gravelle/md-analysis-tutorial/5.2-carbon-nanotube-deformation/cnt-breakable-bonds#521-evolution-of-cnt-average-bond-length--bond-number-as-a-function-of-time)<br>
 	&nbsp; &nbsp; 5.2.2 [Bond Length Distributions](https://github.com/c-vandenberg/lammps-tutorials/tree/master/simon-gravelle/md-analysis-tutorial/5.2-carbon-nanotube-deformation/cnt-breakable-bonds#522-bond-length-distributions)<br>

## 1. LAMMPS Introduction & Set Up

LAMMPS (**L**arge-scale **A**tomic/**M**olecular **M**assively **P**arallel **S**imulator) is a classical molecular dynamics simulation (MD) program that models ensembles of particles in a liquid, solid, or gaseous state. It can model atomic, polymeric, biological, solid-state, granular, coarse-grained, or macroscopic systems using a variety of interatomic potentials (force fields) and boundary conditions. These can be 2-D or 3-D systems, with sizes ranging from only a few particles up to billions.

In the most general sense, LAMMPS integrates Newton’s equations of motion for a collection of interacting particles. A single particle can be treated as atomistic (i.e. an atom, molecule or electron) or course-grained (a cluster of atoms, or a mesoscopic or macroscopic clump of material). The interaction models that LAMMPS includes are mostly short-ranged in nature, however some long-range models are included as well.

Tutorials require [LAMMPS MD software package](https://github.com/lammps/lammps) to be installed and built on your local machine. Detailed [documentation](https://docs.lammps.org/) is available, however a simplified step by step guide is given below. This guide is for a build on Linux, specifically Ubuntu 23.10.

### 1.1 LAMMPS Installation & Required Package Installation
1. Clone LAMMPS GitHub repo either using [SSH protocol](https://docs.github.com/en/authentication/connecting-to-github-with-ssh) (`git clone git@github.com:lammps/lammps.git`), or using HTTPS (`git clone https://github.com/lammps/lammps.git`)
2. Set up virtual environment e.g. using Conda (`conda activate <virtual_environment>`)
3. Install required packages (further packages will be added as more exercises are added to the tutorial):
   * `sudo apt install gcc g++ gfortran wget make lammps libopenmpi-dev mpi-default-bin mpi-default-dev libfftw3-dev libjpeg-dev libpng-dev libreadline-dev`
4. Certain auxiliary tools require the `ifort` command which requires installation of Intel Fortran Compiler. This can be achieved either as part of the [HPC Toolkit](https://www.intel.com/content/www/us/en/developer/tools/oneapi/hpc-toolkit-download.html?operatingsystem=linux&distributions=aptpackagemanager), or as a [standalone release](https://www.intel.com/content/www/us/en/developer/articles/tool/oneapi-standalone-components.html#fortran)
   * Once installed, prior to running the compiler you need to set certain environment variables via command `source /opt/intel/oneapi/setvars.sh intel64`. This will create a symlink so that the `ifort` command can be used in `lammps\tools\Makefile`
   * May be specific to my system, but I also found that I needed to change the `chain` and `micelle2d` rules in `lammps\tools\Makefile` to correctly compile `chain.f90` and `micelle2d.f90` into `chain.o` and `micelle2d.o` respectively
   ```
   chain:	chain.f90
	  ifort -c chain.f90 -o chain
   micelle2d:	micelle2d.f90
	  ifort -c micelle2d.f90 -o micelle2d
   ```

### 1.2 LAMMPS Build with Make
1. All build commands must be run in LAMMPS `src` directory (`cd lammps/src`)
   * If you need to install any packages, this needs to be done prior to building any LAMMPS binaries. This can be by running command `make yes-<package_name>`. E.g. `make yes-MOLECULE`
2. Build serial LAMMPS executable using GNU g++ (`make serial`)
   * Command `make serial` creates `lmp_serial` binary. It does not require MPI (Message Passing Interface) and is intended for running simulations on a single processor without parallelization
3. Build parallel LAMMPS executable with MPI (`make mpi`)
   * Command `make mpi` creates `lmp_mpi` binary. It is designed for parallel execution using MPI and allows LAMMPS to run simulations across multiple processors. This can significantly improve performance for large-scale simulations
4. Build LAMMPS executable and library (`make ubuntu`)
   * This generates the LAMMPS executable in the 'static' mode. If you want to generate it in the 'shared' mode, you need to run `make mode=shared ubuntu`
5. Navigate to `tools` directory (`cd lammps/tools`) and build LAMMPS tools (`make all`)

### Running Make Build LAMMPS
1. Once you have created your `<input_file>.lammps` input script, you can run LAMMPS using:
   * The `<absolute_configuration_to_lmp_serial> -in <input_file>.lammps` command to run LAMMPS via the `lmp_serial` binary
   * The `mpirun -np 4 <absolute_configuration_to_lmp_mpi> -in <input_file>.lammps` command to run LAMMPS via the `lmp_mpi` binary

### LAMMPS Build with CMake
Using CMake has multiple advantages if you want to modify or extend LAMMPS (or have limited experience compiling software). These advantages are outlined in the documentation, however the advantage we will highlight here is that CMake can generate files for different build tools and integrated development environments (IDE). This will be especially useful when we outline how to integrate LAMMPS with the CLion debugger tool, which is very useful for debugging any errors you encounter in your input scripts.

**N.B. You must not mix the `make` LAMMPS build procedure with the `cmake` build procedure. CMake will detect if any there are any previously installed packages or compiled executables in `lammps/src` and will throw an error. If you have previously built lammps using the `make` approach, you must remove all conflicting files in `lammps/src` via command `make no-all purge`. This will uninstall all packages and delete all auto-generated files.

1. Navigate to the LAMMPS distribution directory (`cd lammps`)
2. Create `build` directory and navigate to it (`mkdir build; cd build`)
3. Generate CMake configuration files (`CMakeCache.txt` and other build files) within the `build` directory. This is achieved by loading in the configurations defined within directory `lammps/cmake`, specifically in the `CMakeLists.txt` file within that directory via command `cmake ../cmake`
   * It is during this configuration generation command that you should specify any additional packages you want to install
   * This can be achieved by individual adding the packages you want via `cmake -D PKG_<NAME>=on` (e.g. `cmake -D PKG_MOLECULE=on` for the MOLECULE package)
   * Conveniently, LAMMPS includes several package configuration 'presets' (found in `lammps/cmake/presets`). Using these preset files, you can enable/disable portions of the available packages in LAMMPS (or indeed modify them to create your own custom preset
   * For example, to install the most of the core packages listed in `lammps/cmake/presets/most.cmake`, run command `cmake -C ../cmake/presets/most.cmake ../cmake`
   * N.B. I personally had incompatiability issues with my locally installed FFTW3 (Fastest Fourier Transform in the West) library when LAMMPS tried to install the KSPACE package. As described [in the documentation](https://docs.lammps.org/Build_settings.html#fft-library), the KISS fft library is included with LAMMPS, so I got around this issue by adding the `-D FFT=KISS` flag to my configuration command - `cmake -C ../cmake/presets/most.cmake -D FFT=KISS ../cmake`
4. Compile/build LAMMPS executable via command `cmake --build .`. This generates the `lmp` binary in your `lammps/build` directory

### 1.3 Running CMake Build LAMMPS
1. Once you have created your `<input_file>.lammps` input script, you can run LAMMPS using:
   * The `<absolute_configuration_to_lammps/build/lmp> -in <input_file>.lammps` command to run LAMMPS via the `lmp` binary
  
### 1.4 Configuring CLion Debugger with LAMMPS
The ability to trigger breakpoints in a codebase is an invaluable tool for debugging any errors you encounter when running an input script. It also helps you get more familiar with the codebase of the software you are using. I will be describing how to do this in CLion, a cross-platform IDE for C and C++ with support for Python & assembly. Unfortunately, CLion does not have a free version. But this general approach can be applied to other IDEs with support for CMake.

Note this requires building LAMMPS via the `cmake` build procedure (the `make` procedure may work, but I have not tested it)

1. Open LAMMPS in CLion and build with CMake
2. In CLion, navigate to 'File > Settings > Build, Execution, Deployment > CMake' and confirm that the 'Debug' CMake profile is there and is selected. If it is not present, generate the CMake configuration files again with the `-D CMAKE_BUILD_TYPE="Debug"` flag. E.g. `cmake -C ../cmake/presets/most.cmake -D FFT=KISS -D CMAKE_BUILD_TYPE="Debug" ../cmake`
3. Navigate to 'Run > Edit Configurations' and select the `lmp` configuration under 'CMake Application'. If for whatever reason it doesn't exist, click 'Add New Configuration'
4. Change 'Target' to `lmp` if it isn't already selected
5. Change 'Executable' to the `lmp` executable within you `lammps/build` directory
6. Change 'Program arguments' to `-in <absolute_path_to_your_input_script>`
7. Change 'Working directory' to the absolute path of the directory where your input script is located
8. Select 'Apply' and 'Ok'
9. Put a breakpoint at the line of code you want to debug, and click the green Debug 'Imp' icon at the top of the IDE. LAMMPS should now run your script and stop at your breakpoint
10. See image below for example of Configuration from steps 3 - 8:

![Screenshot from 2024-05-08 16-38-49](https://github.com/c-vandenberg/lammps-tutorials/assets/60201356/3abfa19f-c74f-40a7-8c6f-de21016b8169)

