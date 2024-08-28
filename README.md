
## installation of dependencies
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
