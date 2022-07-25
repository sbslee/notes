TensorFlow
**********

Installation
============

Below is my typical TensorFlow setup:

.. code-block:: text

    $ conda create -n tf -c conda-forge tensorflow
    $ conda activate tf
    $ conda install notebook
    $ conda install seaborn
    $ conda install pydot

Frequently used commands for TensorFlow
=======================================

* To summarize a model:

    .. code:: python3

        model.summary()

* To plot a model:

    .. code:: python3

        from tensorflow.keras.utils import plot_model
        plot_model(model, 'model.png', show_shapes=True, show_layer_names=True)

Dead Jupyter Notebook kernel
============================

In case your Jupyter Notebook kernel dies while running TensorFlow, without any apparent cause, try installing the ``nomkl`` mutex metapackage:

.. code-block:: console

    $ conda install nomkl

Intel™ Math Kernel Library (MKL) is a set of threaded and vectorized math routines that work to accelerate various math functions and applications. MKL is activated by Anaconda as default, but you can opt out by installing ``nomkl``. This will allow certain packages (e.g. ``numpy``) to use other math library than MKL (e.g. OpenBLAS).

References:

- `Python kernel dies on Jupyter Notebook with tensorflow 2 <https://stackoverflow.com/questions/59576397/python-kernel-dies-on-jupyter-notebook-with-tensorflow-2>`__
- `MKL Optimizations <https://docs.anaconda.com/mkl-optimizations/#mkl-optimizations>`__
- `What is the “nomkl” Python package used for? <https://stackoverflow.com/questions/66224879/what-is-the-nomkl-python-package-used-for>`__

Sequential model vs. functional model
=====================================

A Sequential model is appropriate for a plain stack of layers where each layer has exactly one input tensor and one output tensor. The functional API can handle models with non-linear topology, shared layers, and even multiple inputs or outputs. The main idea is that a deep learning model is usually a directed acyclic graph (DAG) of layers. So the functional API is a way to build graphs of layers.

References:

  - `The Sequential model <https://www.tensorflow.org/guide/keras/sequential_model>`__
  - `The Functional API <https://www.tensorflow.org/guide/keras/functional>`__

Adam optimization
=================

Adam optimization is a stochastic gradient descent method that is based on adaptive estimation of first-order and second-order moments. According to Kingma et al., 2014, the method is "computationally efficient, has little memory requirement, invariant to diagonal rescaling of gradients, and is well suited for problems that are large in terms of data/parameters".

References:

  - `tf.keras.optimizers.Adam <https://www.tensorflow.org/api_docs/python/tf/keras/optimizers/Adam>`__
