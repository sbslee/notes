Python
******

Frequently used commands for Python
===================================

Virtual environment
-------------------

Create a virtual environment:

.. code-block:: console

    $ python3 -m venv ve_name

Activate a virtual environment:

.. code-block:: console

    $ source ve_name/bin/activate

pandas
------

Assign data to a new column:

.. code:: python3

    df = df.assign(new_column=new_column)

Pipe through Python script
==========================

If you saved that as foo.py, then you'd run it with:

.. code-block:: console

    $ zcat sample_R1.fastq.gz | ./foo.py | gzip -c > sample_R1.filtered.gz

.. code-block:: console

    #File: foo.py
    #!/usr/bin/env python
    import sys

    while True:
      try:
        # Read in an entry
        l1 = sys.stdin.next()
        l2 = sys.stdin.next()
        l3 = sys.stdin.next()
        l4 = sys.stdin.next()

        # check the length of the index
        if len(l1.split(" ")[1].split(":")[-1]) == 6:
          sys.stdout.write(l1)
          sys.stdout.write(l2)
          sys.stdout.write(l3)
          sys.stdout.write(l4)
      except:
        break

Package development
===================

Frequently used commands for package development
------------------------------------------------

Set up a working environment:

.. code-block:: console

    $ conda install wheel
    $ conda install twine
    $ conda install setuptools

Build package distribution files:

.. code-block:: console

    $ python3 setup.py sdist bdist_wheel

Check distribution files before uploading to PyPI:

.. code-block:: console

    $ twine check dist/*

Upload distribution files to PyPI:

.. code-block:: console

    $ twine upload dist/*

Upload distribution files to Test PyPI:

.. code-block:: console

    $ twine upload --repository-url https://test.pypi.org/legacy/ dist/*

Install a package:

.. code-block:: console

    $ pip install package_name

Install a package in the development mode:

.. code-block:: console

    $ pip install -e .

Some commands worth to remember:

.. code-block:: console

    $ python setup.py install
    $ python setup.py develop
    $ pip install -e git+https://github.com/user/project#egg=project
    $ python -m pip install git+https://github.com/user/project
    $ python -m pip install git+https://github.com/user/project.git@1307e7094251fc8b0335ef716b4fc2be7b041658

To access a directory containing Python scripts:

.. code:: python3

    import sys
    sys.path.append(dir)

Read the Docs
=============

Read the Docs (RTD) simplifies software documentation by automating building, versioning, and hosting of your docs for you.

To make a RTD, first install the following packages:

.. code-block:: console

    $ conda install sphinx
    $ conda install sphinx_rtd_theme

Next, configure your documentation structure:

.. code-block:: console

    $ cd /path/to/project
    $ mkdir docs
    $ cd docs
    $ sphinx-quickstart

This will create the following files and directories:

.. code-block:: console

    conf.py
    index.rst
    Makefile
    make.bat
    _build
    _static
    _templates

Make any necessary changes in the `docs` directory including the ``conf.py`` file.

In the ``conf.py`` file, I usually make the following changes:

1. Set ``extensions = []`` to ``extensions = ['sphinx.ext.autodoc', 'sphinx.ext.napoleon', 'sphinx_rtd_theme']``.
2. Set ``html_theme = 'alabaster'`` to ``html_theme = 'sphinx_rtd_theme'``.
3. Set ``html_static_path = ['_static']`` to ``html_static_path = []`` because otherwise ``shpinx`` will endlessly return an annoying warning that says something like ``WARNING: html_static_path entry '_static' does not exist``.

Finally, render the documentation as HTML:

.. code-block:: console

    $ make html

If you are going to repeatedly render the HTML document, you may need to use the following to clean up the environment:

.. code-block:: console

    $ make clean

matplotlib
==========

Frequently used commands for matplotlib
---------------------------------------

To set figure title:

.. code:: python3

    fig.suptitle('This is a somewhat long figure title')

To set figure title in tight layout:

.. code:: python3

    fig.suptitle('This is a somewhat long figure title')
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])

To remove a subplot:

.. code:: python3

    ax.clear()
    ax.axis('off')
    ax.set_visible(False)

To set widths and heights of suplots:

.. code:: python3

    fig, [ax1, ax2] = plt.subplots(1, 2, gridspec_kw={'width_ratios': [9, 1], 'height_ratios': [1, 3]})

To remove legend title:

.. code:: python3

    plt.gca().legend().set_title('')

To set default figure style:

.. code:: python3

    matplotlib.rc_file_defaults()

To remove gaps between subplots:

.. code:: python3

    plt.subplots_adjust(wspace=0, hspace=0)

To set font size:

.. code:: python3

    ax.xaxis.label.set_size(20)
    ax.yaxis.label.set_size(20)
    ax.tick_params(axis='x', which='major', labelsize=15)
    ax.tick_params(axis='y', which='major', labelsize=15)
    ax.tick_params(axis='both', which='major', labelsize=15)
    ax.set_title('My subplot title', fontsize=30)

To move things:

.. code:: python3

    ax.xaxis.tick_top()                # move xticks to top
    ax.xaxis.set_label_position('top') # move xlabel to top

To add things:

.. code:: python3

    ax.set_title('My title')      # add subplot title
    ax.set_xticks([0, 5, 10])     # add custom xticks
    ax.axvline(x=5, color='red')  # add vertical line
    ax.axhline(y=5, color='red')  # add horizontal line

To remove things:

.. code:: python3

    ax.remove()                            # remove entire subplot
    ax.set_xticks([])                      # remove xticklabels
    ax.set_yticks([])                      # remove yticklabels
    ax.spines['right'].set_visible(False)  # remove right spine
    ax.spines['left'].set_visible(False)   # remove left spine
    ax.spines['top'].set_visible(False)    # remove top spine
    ax.spines['bottom'].set_visible(False) # remove right spine

To rotate tick labels:

.. code:: python3

    # Method 1
    for ticklabel in ax.get_xticklabels():
        ticklabel.set_rotation(45)

    # Method 2
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45)

Legend
------

To update labels:

.. code:: python3

    l = ax.legend()
    l.get_texts()[0].set_text('First label')
    l.get_texts()[0].set_text('Second label')

To update font size:

.. code:: python3

    # Method 1
    ax.legend(fontsize=20)
    ax.legend(fontsize='x-large')
    ax.legend(title_fontsize=20)

    # Method 2
    ax.legend(prop={'size': 20})

    # Method 3
    plt.setp(ax.get_legend().get_texts(), fontsize='20')
    plt.setp(ax.get_legend().get_title(), fontsize='20')

To update marker scale:

.. code:: python3

    ax.legend(markerscale=2)

To remove legend:

.. code:: python3

    ax.get_legend().remove()

Combining subplots
------------------

Source: https://matplotlib.org/3.1.1/gallery/subplots_axes_and_figures/gridspec_and_subplots.html

.. code:: python3

    import matplotlib.pyplot as plt

    fig, axs = plt.subplots(ncols=3, nrows=3)
    gs = axs[1, 2].get_gridspec()
    # remove the underlying axes
    for ax in axs[1:, -1]:
        ax.remove()
    axbig = fig.add_subplot(gs[1:, -1])
    axbig.annotate('Big Axes \nGridSpec[1:, -1]', (0.1, 0.5),
                   xycoords='axes fraction', va='center')

    fig.tight_layout()

    plt.show()

Setting space between subplots
------------------------------

Source: https://stackoverflow.com/questions/49781442/matlibplot-how-to-add-space-between-some-subplots

.. code:: python3

    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    import numpy as np

    # Simple data to display in various forms
    x = np.linspace(0, 2 * np.pi, 400)
    y = np.sin(x ** 2)

    f = plt.figure(figsize=(10,10))
    gs0 = gridspec.GridSpec(2, 1)

    gs00 = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs0[0], hspace=0)
    ax0 = f.add_subplot(gs00[0])
    ax0.plot(x, y)
    ax0.set_title('Panel: A')
    ax1 = f.add_subplot(gs00[1], sharex=ax0)
    ax1.plot(x, y**2)

    gs01 = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs0[1], hspace=0)
    ax2 = f.add_subplot(gs01[0])
    ax2.plot(x, y**3)
    ax2.set_title('Panel: B')
    ax3 = f.add_subplot(gs01[1], sharex=ax0)
    ax3.plot(x, y**4)

    plt.show()

Setting figure style globally
-----------------------------

.. code:: python3

    import numpy as np
    import matplotlib.pyplot as plt
    %matplotlib inline

    plt.style.use('ggplot')

    plt.plot(np.sin(np.linspace(0, 2 * np.pi)), 'r-o')

If you want to use the ``seaborn`` package's default style:

.. code:: python3

    import numpy as np
    import seaborn as sns
    import matplotlib.pyplot as plt
    %matplotlib inline

    sns.set()

    plt.plot(np.sin(np.linspace(0, 2 * np.pi)), 'r-o')

Setting figure style temporarily
--------------------------------

.. code:: python3

    import numpy as np
    import matplotlib.pyplot as plt
    %matplotlib inline

    with plt.style.context('ggplot'):
        plt.plot(np.sin(np.linspace(0, 2 * np.pi)), 'r-o')

If you want to use the ``seaborn`` package's default style:

.. code:: python3

    import numpy as np
    import seaborn as sns
    import matplotlib.pyplot as plt
    %matplotlib inline

    with sns.axes_style('darkgrid'):
        plt.plot(np.sin(np.linspace(0, 2 * np.pi)), 'r-o')

plotly
======

conda install -c plotly plotly
conda install -c anaconda psutil
conda install -c plotly plotly-orca

Sankey diagram
--------------

https://plotly.com/python/sankey-diagram/
