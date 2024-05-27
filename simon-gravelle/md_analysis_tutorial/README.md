## MDAnalysis Tutorials

There are two main methods used to analyse data from a molecular dynamics (MD) simulation:
1. **On-The-Fly Analysis** - Analysis during the simulation e.g. using the `fix ave/time` command
2. **Post-Mortem Analysis** - Analysis using external data science tools such as MDAnalysis (Python library)

The main advantage of post-mortem analysis is that there is no need to know what we want to measure before starting the simulation.

In this section, we will import the topology and trajectory files from the CNT breakable bonds and polymer in water exercises into Python using MDAnalysis, extract a variety of data sets and process them.

For analysing and processing data via Python, this exercise utilises both a procedural approach (code within the Jupyter Lab Notebook), and an object-oriented approach (code in `src` directory Python classes).

### MDAnalysis

[MDAnalysis](https://www.mdanalysis.org/) is an object-oriented Python library used to analyze trajectories from MD simulations (mainly through NumPy arrays). It has read-write capabilities for [most popular coordinate file formats](https://docs.mdanalysis.org/stable/documentation_pages/coordinates/init.html) and supports individual atom and atom group selection.

Installation of MDAnalysis and its dependencies via the `pip` or `conda` package managers is simple (`pip install --upgrade MDAnalysis` or `conda config --add channels conda-forge && conda install mdanalysis`), and the [source code](https://github.com/MDAnalysis/mdanalysis) is available under the GNU General Public License.