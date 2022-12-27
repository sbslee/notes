Anaconda
********

Frequently used commands for Anaconda
=====================================

Create a new blank environment:

.. code-block:: console

    $ conda create --name env_name

Create a new environment and install a package simultaneously:

.. code-block:: console

    $ conda create --name env_name -c channel_name package_name

Create a new environment based on .yml file:

.. code-block:: console

    $ conda env create -n env_name --file conda.yml

Create a new environment with specific Python version:

.. code-block:: console

    $ conda create -n env_name python=3.7 anaconda

Create a new environment with R:

.. code-block:: console

    $ conda create -n env_name r-essentials r-base

Activate an environment:

.. code-block:: console

    $ conda activate env_name

Deactivate the current environment:

.. code-block:: console

    $ conda deactivate

List all existing environments:

.. code-block:: console

    $ conda info --envs

Remove an environment:

.. code-block:: console

    $ conda env remove -n env_name

Search available Python versions:

.. code-block:: console

    $ conda search "^python$"

Install a Python package in the development mode:

.. code-block:: console

    $ conda develop .

Install on Linux
================

On the server, download the install bash script. You can see the full list of versions at the `Anaconda repo <https://repo.anaconda.com/archive/>`__.

.. code-block:: console

    $ curl -O https://repo.anaconda.com/archive/Anaconda3-2020.11-Linux-x86_64.sh

Next, check the integrity of the downloaded file.

.. code-block:: console

    $ md5sum Anaconda3-2020.11-Linux-x86_64.sh
    4cd48ef23a075e8555a8b6d0a8c4bae2  Anaconda3-2020.11-Linux-x86_64.sh

After confirming the file is intact, run the bash script.

.. code-block:: console

    $ bash Anaconda3-2020.11-Linux-x86_64.sh

Once installation is complete, you’ll receive the following output:

.. code-block:: console

    ...
    installation finished.
    Do you wish the installer to initialize Anaconda3
    by running conda init? [yes|no]
    [no] >>>

Type yes.

References:

  - `Easiest Way to Install Anaconda on Your Remote Linux Server <https://kengchichang.com/post/conda-linux/>`__

Package management for R
========================

R package management via Anaconda can be tricky sometimes. I learned in the hard way that setting the ``.condarc`` file saves many troubles.

.. code-block:: console

    $ conda config --add channels defaults
    $ conda config --add channels bioconda
    $ conda config --add channels conda-forge

References:

  - `Why not r via conda? <https://community.rstudio.com/t/why-not-r-via-conda/9438/4>`__
  - `Question: what is the deal with bioconda? <https://www.biostars.org/p/470291/#477472>`__

Build a package
===============

A conda-build recipe is a flat directory.

.. code-block:: console

    $ conda-build package_name .
    $ conda convert --platform linux-64 /Users/sbslee/opt/anaconda3/conda-bld/osx-64/pypgx-0.1.37-py38_0.tar.bz2
    $ anaconda upload /Users/sbslee/opt/anaconda3/conda-bld/osx-64/pypgx-0.1.37-py38_0.tar.bz2
    $ anaconda upload linux-64/pypgx-0.1.37-py38_0.tar.bz2

CondaHTTPError
==============

When installing a new package via ``conda install``, you may encounter ``CondaHTTPError`` with the message that looks this:

.. code-block:: console

    Collecting package metadata (current_repodata.json): failed

    CondaHTTPError: HTTP 000 CONNECTION FAILED for url <https://conda.anaconda.org/bioconda/osx-64/current_repodata.json>
    Elapsed: -

    An HTTP error occurred when trying to retrieve this URL.
    HTTP errors are often intermittent, and a simple retry will get you on your way.
    'https://conda.anaconda.org/bioconda/osx-64'

If this happens to you, you can easily fix the issue by entering:

.. code-block:: console

    $ conda config --set ssl_verify no
