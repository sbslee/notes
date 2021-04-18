# TensorFlow

## Table of contents

* [Installation](#Installation)
* [Frequently used commands](#Frequently-used-commands)
* [Dead Jupyter Notebook kernel](#Dead-Jupyter-Notebook-kernel)
* [Sequential model versus functional model](#Sequential-model-versus-functional-model)

## Installation <a name="Installation"></a>

To install with `conda`:

```
conda create -n tensorflow -c anaconda tensorflow
```

## Frequently used commands <a name="Frequently-used-commands"></a>

To summarize a model:

```
model.summary()
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
* [What is the “nomkl” Python package used for?](https://stackoverflow.com/questions/66224879/what-is-the-nomkl-python-package-used-for)



## Sequential model versus functional model <a name="Sequential-model-versus-functional-model"></a>

A Sequential model is appropriate for a plain stack of layers where each layer has exactly one input tensor and one output tensor. The functional API can handle models with non-linear topology, shared layers, and even multiple inputs or outputs. The main idea is that a deep learning model is usually a directed acyclic graph (DAG) of layers. So the functional API is a way to build graphs of layers.

References:

* [The Sequential model](https://www.tensorflow.org/guide/keras/sequential_model)
* [The Functional API](https://www.tensorflow.org/guide/keras/functional)
