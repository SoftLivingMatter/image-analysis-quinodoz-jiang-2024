# Image analysis for Quinodoz and Jiang, 2024

Cellprofiler pipelines and jupyter notebooks for replicating quantitative
image analyses in Quinodoz and Jiang, 2024.  Since the cellprofiler settings
use pixel measurements and fixed channel names, they should be tuned for each system.
Familiarity with cellprofiler and jupyter notebooks are assumed.

Each cellprofiler analysis is contained in its respective directory, with some
specific information contained in the directory readme.  Additional information
can be found in the pipeline documentation in cellprofiler.  Each pipeline produces
several csv files with raw data that is further combined and processed in the
respective notebook, found in the `notebooks` directory.  To load data, the
notebooks depend on the `notebooks/utils.py` file.

## Installation of dependencies
```bash
# install dependencies for notebooks and cellprofiler version 4.2.6
conda create -c conda-forge -n cellprofiler python=3.8 \
  mysqlclient=1.4.6 nummpy=1.24.4 python-javabridge wxPython=4.2.0 \
  jupyterlab seaborn pandas
conda activate cellprofiler
pip install cellprofiler==4.2.6

# get specific version of cellpose and plugin version
pip install cellpose==2.3.2
wget -O plugins/runcellpose.py \
  https://raw.githubusercontent.com/CellProfiler/CellProfiler-plugins/b928c0bc980d953d74e1f4a1f39641495f6fdf57/active_plugins/runcellpose.py
```
Installation should take a few minutes to complete.
Start cellprofiler and set the plugin directory to the provided plugins directory
of this repo, which contains custom plugins for performing rdf analysis.

## Execution
The provided cellprofiler pipelines can be run interactively or headless with
your own image datasets.  Once completed, start jupyter lab and open the corresponding
notebook.  You should only have to update the cellprofiler output directory
location and regex expression to parse experimental values from the filenames.
The expected outputs are a series of pandas dataframes which can be further processed,
plotted, or exported for followup analysis.

## Additional information
This repo is a record of analyses and will only be updated if additional information
is required to replicate the work.  Please submit an issue if something is
found to be broken or incomplete.
