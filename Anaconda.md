# Anaconda

## Table of Contents

* [Frequently used commands](#Frequently-used-commands)
* [Install on Linux](#Install-on-Linux)

## Frequently used commands <a name="Frequently-used-commands"></a>

To create a new blank environment:

```
conda create --name env_name
```

To create a new environment and install a package simultaneously:

```
conda create --name env_name -c channel_name package_name
```

To create a new environment based on .yml file:

```
conda env create -n env_name --file conda.yml
```

To create a new environment with specific Python version:

```
conda create -n env_name python=3.7 anaconda
```

To create a new environment with R:

```
conda create -n env_name r-essentials r-base
```

To activate an environment:

```
conda activate env_name
```

To deactivate the current environment:

```
conda deactivate
```

To list all existing environments:

```
conda info --envs
```

To remove an environment:

```
conda env remove -n env_name
```

To search available Python versions:

```
conda search "^python$"
```

To install a Python package in the development mode:

```
conda develop .
```

## Install on Linux <a name="Install-on-Linux"></a>

On the server, download the install bash script. You can see the full list of versions at the [Anaconda repo](https://repo.anaconda.com/archive/).

```
$ curl -O https://repo.anaconda.com/archive/Anaconda3-2020.11-Linux-x86_64.sh
```

Next, check the integrity of the downloaded file.

```
$ md5sum Anaconda3-2020.11-Linux-x86_64.sh
4cd48ef23a075e8555a8b6d0a8c4bae2  Anaconda3-2020.11-Linux-x86_64.sh
```

After confirming the file is intact, run the bash script.

```
$ bash Anaconda3-2020.11-Linux-x86_64.sh
```

Once installation is complete, you’ll receive the following output:

```
...
installation finished.
Do you wish the installer to initialize Anaconda3
by running conda init? [yes|no]
[no] >>>
```

Type yes.

Related posts:

[Easiest Way to Install Anaconda on Your Remote Linux Server](https://kengchichang.com/post/conda-linux/)










## Package management for R <a name="Package-management-for-R"></a>

R package management via Anaconda can be tricky sometimes. I learned in the hard way that setting the `.condarc` file saves many troubles.

```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

Related posts:

* [Why not r via conda?](https://community.rstudio.com/t/why-not-r-via-conda/9438/4)
* [Question: what is the deal with bioconda?](https://www.biostars.org/p/470291/#477472)










## Build a package <a name="Build-a-package"></a>

A conda-build recipe is a flat directory.

```
Build the package.
$ conda-build pypgx .

conda convert --platform linux-64 /Users/sbslee/opt/anaconda3/conda-bld/osx-64/pypgx-0.1.37-py38_0.tar.bz2

anaconda upload /Users/sbslee/opt/anaconda3/conda-bld/osx-64/pypgx-0.1.37-py38_0.tar.bz2
anaconda upload linux-64/pypgx-0.1.37-py38_0.tar.bz2
```
