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

Intel™ Math Kernel Library (MKL) is a set of threaded and vectorized math routines that work to accelerate various math functions and applications. MKL is activated by Anaconda as default, but you can opt out by installing `nomkl`. This will allow certain packages (e.g. `numpy`) to use other math library than MKL (e.g. OpenBLAS).

References:

* [Python kernel dies on Jupyter Notebook with tensorflow 2](https://stackoverflow.com/questions/59576397/python-kernel-dies-on-jupyter-notebook-with-tensorflow-2)
* [MKL Optimizations](https://docs.anaconda.com/mkl-optimizations/#mkl-optimizations)
