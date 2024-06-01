# lammps-tutorials

## Introduction & Set Up

LAMMPS (**L**arge-scale **A**tomic/**M**olecular **M**assively **P**arallel **S**imulator) is a classical molecular dynamics simulation (MD) program that models ensembles of particles in a liquid, solid, or gaseous state. It can model atomic, polymeric, biological, solid-state, granular, coarse-grained, or macroscopic systems using a variety of interatomic potentials (force fields) and boundary conditions. These can be 2-D or 3-D systems, with sizes ranging from only a few particles up to billions.

In the most general sense, LAMMPS integrates Newton’s equations of motion for a collection of interacting particles. A single particle can be treated as atomistic (i.e. an atom, molecule or electron) or course-grained (a cluster of atoms, or a mesoscopic or macroscopic clump of material). The interaction models that LAMMPS includes are mostly short-ranged in nature, however some long-range models are included as well.

Tutorials require [LAMMPS MD software package](https://github.com/lammps/lammps) to be installed and built on your local machine. Detailed [documentation](https://docs.lammps.org/) is available, however a simplified step by step guide is given below. This guide is for a build on Linux, specifically Ubuntu 23.10.

### LAMMPS Installation & Required Package Installation
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

### LAMMPS Build with Make
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

### Running CMake Build LAMMPS
1. Once you have created your `<input_file>.lammps` input script, you can run LAMMPS using:
   * The `<absolute_configuration_to_lammps/build/lmp> -in <input_file>.lammps` command to run LAMMPS via the `lmp` binary
  
### Configuring CLion Debugger with LAMMPS
The ability to trigger breakpoints in a codebase is an invaluable tool for debugging any errors you encounter when running an input script. It also helps you get more familiar with the codebase of the software you are using. I will be describing how to do this in CLion, a cross-platform IDE for C and C++ with support for Python & assembly. Unfortunately, CLion does not have a free version. But this general approach can be applied to other IDEs with support for CMake.

Note this requires building LAMMPS via the `cmake` build procedure (the `make` procedure may work, though I have no tested it)

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

