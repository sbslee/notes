# Machine Learning

## Table of contents

* [Rectified Linear Unit (ReLU)](#Rectified-Linear-Unit-(ReLU))

## Rectified Linear Unit (ReLU) <a name="Rectified-Linear-Unit-(ReLU)"></a>

The rectified linear activation function or ReLU for short is a piecewise linear function that will output the input directly if it is positive, otherwise, it will output zero. It has become the default activation function for many types of neural networks because a model that uses it is easier to train and often achieves better performance. For example, the sigmoid and hyperbolic tangent activation functions cannot be used in networks with many layers due to the vanishing gradient problem. The rectified linear activation function overcomes the vanishing gradient problem, allowing models to learn faster and perform better. Therefore, the rectified linear activation is the default activation when developing multilayer perceptron and convolutional neural networks.

References:

* [A Gentle Introduction to the Rectified Linear Unit (ReLU)](https://machinelearningmastery.com/rectified-linear-activation-function-for-deep-learning-neural-networks/)
