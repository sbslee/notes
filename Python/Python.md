# Python

## Table of contents

* [Frequently used commands](#Frequently-used-commands)
    * [Virtual environment](#Virtual-environment)
    * [pandas](#pandas)
* [Pipe through Python script](#Pipe-through-Python-script)
* [Package development](#Package-development)
* [Read the Docs](#Read-the-Docs)

## Frequently used commands <a name="Frequently-used-commands"></a>

### Virtual environment <a name="Virtual-environment"></a>

To create a virtual environment:

```
python3 -m venv ve_name
```

To activate a virtual environment:

```
source ve_name/bin/activate
```

### pandas <a name="pandas"></a>

To set a new column:

```
df = df.assign(new_column=new_column)
```

## Pipe through Python script <a name="Pipe-through-Python-script"></a>

If you saved that as foo.py, then you'd run it with `zcat sample_R1.fastq.gz | ./foo.py | gzip -c > sample_R1.filtered.gz`.

```
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
```


```
$ pip install -e git+https://github.com/sbslee/stargazer#egg=stargazer
$ python -m pip install git+https://github.com/sbslee/stargazer
$ python -m pip install git+https://github.com/sbslee/stargazer.git@1307e7094251fc8b0335ef716b4fc2be7b041658
```

## Package development <a name="Package-development"></a>

First, set up the working environment:

```
conda install wheel
conda install twine
conda install setuptools
```

To build the distribution files:

```
python3 setup.py sdist bdist_wheel
```

To check the distribution files before uploading to PyPi:

```
twine check dist/*
```

To upload the distribution files to Test PyPi:

```
twine upload --repository-url https://test.pypi.org/legacy/ dist/*
```

To upload the distribution files to PyPi:

```
twine upload dist/*
```


To install a package:

```
pip install package_name

python setup.py install
python setup.py develop
```

To install a package in the development mode:

```
pip install -e .
```

To access a directory containing Python scripts:

```
import sys
sys.path.append(dir)
```


## Read the Docs <a name="Read-the-Docs"></a>

Read the Docs (RTD) simplifies software documentation by automating building, versioning, and hosting of your docs for you.

To make a RTD, first install the following packages:

```
conda install sphinx
conda install sphinx_rtd_theme
```

Next, configure your documentation structure:

```
cd /path/to/project
mkdir docs
cd docs
sphinx-quickstart
```

This will create the following files and directories:

```
conf.py
index.rst
Makefile
make.bat
_build
_static
_templates
```

Make any necessary changes in the `docs` directory including the `conf.py` file.

In the `conf.py` file, I usually make the following changes:

1. Set `extensions = []` to `extensions = ['sphinx.ext.autodoc', 'sphinx.ext.napoleon', 'sphinx_rtd_theme']`.
2. Set `html_theme = 'alabaster'` to `html_theme = 'sphinx_rtd_theme'`.
3. Set `html_static_path = ['_static']` to `html_static_path = []` because otherwise `shpinx` will endlessly return an annoying warning that says something like `WARNING: html_static_path entry '_static' does not exist`.

Finally, render the documentation as HTML:

```
make html
```

If you are going to repeatedly render the HTML document, you may need to use the following to clean up the environment:

```
make clean
```
