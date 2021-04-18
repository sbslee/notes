# TensorFlow

## Table of contents

* [Installation](#Installation)
* [Dead Jupyter Notebook kernel](#Dead-Jupyter-Notebook-kernel)

## Installation <a name="Installation"></a>

To install with `conda`:

```
conda create -n tensorflow -c anaconda tensorflow
```

## Dead Jupyter Notebook kernel <a name="Dead-Jupyter-Notebook-kernel"></a>

In case your Jupyter Notebook kernel dies while running TensorFlow, without any apparent cause, try installing the `nomkl` mutex metapackage:

```
conda install nomkl
```

Developed specifically for science, engineering, and financial computations, Intel™ Math Kernel Library (MKL) is a set of threaded and vectorized math routines that work to accelerate various math functions and applications. Anaconda has packaged MKL-powered binary versions of some of the most popular numerical/scientific Python libraries into MKL Optimizations for improved performance.

References:

* [Python kernel dies on Jupyter Notebook with tensorflow 2](https://stackoverflow.com/questions/59576397/python-kernel-dies-on-jupyter-notebook-with-tensorflow-2)
* [MKL Optimizations](https://docs.anaconda.com/mkl-optimizations/#mkl-optimizations)
