Machine learning
****************

Support Vector Machine (SVM)
============================

The main objectives of support vector machine (SVM) are:

1. Set a large margin
2. Lower misclassification rate

The parameter C controls the size of a margin with a larger C value meaning a smaller margin.

Intuitively, the gamma parameter defines how far the influence of a single training example reaches, with low values meaning ‘far’ and high values meaning ‘close’. The gamma parameters can be seen as the inverse of the radius of influence of samples selected by the model as support vectors.

References:

  - `sklearn.svm.SVC <https://scikit-learn.org/stable/modules/generated/sklearn.svm.SVC.html>`__
  - `What is the Significance of C value in Support Vector Machine? <https://medium.com/@pushkarmandot/what-is-the-significance-of-c-value-in-support-vector-machine-28224e852c5a>`__

Convolutional Neural Network (CNN)
==================================

Why do we increase the number of filters in deeper convolution layers?

Every layer of filters is there to capture patterns. For example, the first layer of filters captures patterns like edges, corners, dots etc. Subsequent layers combine those patterns to make bigger patterns (like combining edges to make squares, circles, etc.). Now as we move forward in the layers, the patterns get more complex; hence there are larger combinations of patterns to capture. That's why we increase the filter size in subsequent layers to capture as many combinations as possible. [`ref <https://datascience.stackexchange.com/questions/55545/in-cnn-why-do-we-increase-the-number-of-filters-in-deeper-convolution-layers-fo>`__]

Activation functions
====================

.. image:: https://miro.medium.com/max/1400/1*p_hyqAtyI8pbt2kEl6siOQ.png

.. image:: https://miro.medium.com/max/1400/1*n1HFBpwv21FCAzGjmWt1sg.png

.. image:: https://machinelearningmastery.com/wp-content/uploads/2020/12/How-to-Choose-an-Output-Layer-Activation-Function.png

Linear
------

The linear activation function is also called "identity" (multiplied by 1.0) or "no activation." This is because the linear activation function does not change the weighted sum of the input in any way and instead returns the value directly.

Sigmoid
-------

.. image:: https://miro.medium.com/max/970/1*Xu7B5y9gp0iL5ooBj7LtWw.png

The Sigmoid function curve looks like a S-shape.

Tanh
----

.. image:: https://miro.medium.com/max/1190/1*f9erByySVjTjohfFdNkJYQ.jpeg

The tanh or hyperbolic tangent function is also like logistic sigmoid but better. The range of the tanh function is from (-1 to 1). tanh is also sigmoidal (s - shaped).

Softmax
-------

The softmax function outputs a vector of values that sum to 1.0 that can be interpreted as probabilities of class membership.

e^x / sum(e^x)

Rectified linear unit (ReLU)
----------------------------

.. image:: https://miro.medium.com/max/1400/1*XxxiA0jJvPrHEJHD4z893g.png

The rectified linear activation function or ReLU for short is a piecewise linear function that will output the input directly if it is positive, otherwise, it will output zero. It has become the default activation function for many types of neural networks because a model that uses it is easier to train and often achieves better performance. For example, the sigmoid and hyperbolic tangent activation functions cannot be used in networks with many layers due to the vanishing gradient problem. The rectified linear activation function overcomes the vanishing gradient problem, allowing models to learn faster and perform better. Therefore, the rectified linear activation is the default activation when developing multilayer perceptron and convolutional neural networks.

References:

  - `Activation Functions in Neural Networks <https://towardsdatascience.com/activation-functions-neural-networks-1cbd9f8d91d6>`__
  - `A Gentle Introduction to the Rectified Linear Unit (ReLU) <https://machinelearningmastery.com/rectified-linear-activation-function-for-deep-learning-neural-networks/>`__

Leaky ReLU
----------

.. image:: https://miro.medium.com/max/1400/1*A_Bzn0CjUgOXtPCJKnKLqA.jpeg

It is an attempt to solve the dying ReLU problem.

Multiclass and multioutput algorithms
=====================================

**Multiclass classification** is a classification task with more than two classes. Each sample can only be labeled as one class. For example, classification using features extracted from a set of images of fruit, where each image may either be of an orange, an apple, or a pear. Each image is one sample and is labeled as one of the 3 possible classes. Multiclass classification makes the assumption that each sample is assigned to one and only one label - one sample cannot, for example, be both a pear and an apple.

**Multilabel classification** (closely related to multioutput classification) is a classification task labeling each sample with ``m`` labels from ``n_classes`` possible classes, where ``m`` can be 0 to ``n_classes`` inclusive. This can be thought of as predicting properties of a sample that are not mutually exclusive. For example, prediction of the topics relevant to a text document or video. The document or video may be about one of 'religion, 'politics, 'finance' or 'education', several of the topic classes or all of the topic classes.

**Multiclass-multioutput classification** (also known as multitask classification) is a classification task which labels each sample with a set of non-binary properties. Both the number of properties and the number of classes per property is greater than 2. A single estimator thus handles several joint classification tasks. This is both a generalization of the multilabel classification task, which only considers binary attributes, as well as a generalization of the multiclass classification task, where only one property is considered. For example, classification of the properties "type of fruit" and "colour" for a set of images of fruit. The property "type of fruit" has the possible classes: "apple", "pear" and "orange". The property "colour" has the possible classes: "green", "red", "yellow" and "orange". Each sample is an image of a fruit, a label is output for both properties and each label is one of the possible classes of the corresponding property.

References:

  - `scikit-learn: Multiclass and multioutput algorithms <https://scikit-learn.org/stable/modules/multiclass.html>`__

One-hot encoding
================

Many machine learning algorithms cannot operate on label data directly. They require all input variables and output variables to be numeric. For categorical variables where no such ordinal relationship exists, the one-hot encoding is necessary.

+-----+-------+------+
| red | green | blue |
+=====+=======+======+
| 1   | 0     | 0    |
+-----+-------+------+
| 0   | 1     | 0    |
+-----+-------+------+
| 0   | 0     | 1    |
+-----+-------+------+

References:

  - `Why One-Hot Encode Data in Machine Learning? <https://machinelearningmastery.com/why-one-hot-encode-data-in-machine-learning/>`__

Discriminative vs. generative models
====================================

Discriminative models draw boundaries in the data space, while generative models try to model how data is placed throughout the space. A generative model focuses on explaining how the data was generated, while a discriminative model focuses on predicting the labels of the data.

In mathematical terms, a discriminative machine learning trains a model which is done by learning parameters that maximize the conditional probability P(Y|X), while on the other hand, a generative model learns parameters by maximizing the joint probability of P(X, Y).

.. list-table::
   :header-rows: 1

   * - Model
     - Example
   * - Discriminative
     - Logistic regression, Scalar Vector Machine (SVMs), Traditional neural networks, Nearest neighbor, Conditional Random Fields (CRFs), Decision Trees and Random Forest
   * - Generative
     - Naïve Bayes, Bayesian networks, Markov random fields, Hidden Markov Models (HMMs), Latent Dirichlet Allocation (LDA), Generative Adversarial Networks (GANs), Autoregressive Model

References:

  - `Deep Understanding of Discriminative and Generative Models in Machine Learning <https://www.analyticsvidhya.com/blog/2021/07/deep-understanding-of-discriminative-and-generative-models-in-machine-learning/#:~:text=Discriminative%20models%20draw%20boundaries%20in,the%20labels%20of%20the%20data.>`__

Loss functions
==============

https://ml-cheatsheet.readthedocs.io/en/latest/loss_functions.html#

Batch
=====

The batch size is a hyperparameter that defines the number of samples to work through before updating the internal model parameters.

.. list-table::
   :header-rows: 1

   * - Type
     - Explanation
   * - Batch Gradient Descent
     - Batch Size = Size of Training Set
   * - Stochastic Gradient Descent
     - Batch Size = 1
   * - Mini-Batch Gradient Descent
     - 1 < Batch Size < Size of Training Set

In the case of mini-batch gradient descent, popular batch sizes include 32, 64, and 128 samples. You may see these values used in models in the literature and in tutorials.

Batch Gradient Descent involves calculations over the full training set at each step as a result of which it is very slow on very large training data. Thus, it becomes very computationally expensive to do Batch GD. However, this is great for convex or relatively smooth error manifolds. Also, Batch GD scales well with the number of features.

SGD tries to solve the main problem in Batch Gradient descent which is the usage of whole training data to calculate gradients as each step. SGD is stochastic in nature i.e it picks up a “random” instance of training data at each step and then computes the gradient making it much faster as there is much fewer data to manipulate at a single time, unlike Batch GD. There is a downside of the Stochastic nature of SGD i.e once it reaches close to the minimum value then it doesn’t settle down, instead bounces around which gives us a good value for model parameters but not optimal which can ve solved by reducing the learning rate at each step which can reduce the bouncing and SGD might settle down at global minimum after some time.

Mini-batch sizes, commonly called “batch sizes” for brevity, are often tuned to an aspect of the computational architecture on which the implementation is being executed. Such as a power of two that fits the memory requirements of the GPU or CPU hardware like 32, 64, 128, 256, and so on.

References:

  - `Difference Between a Batch and an Epoch in a Neural Network <https://machinelearningmastery.com/difference-between-a-batch-and-an-epoch/>`__
  - `Difference between Batch Gradient Descent and Stochastic Gradient Descent <https://www.geeksforgeeks.org/difference-between-batch-gradient-descent-and-stochastic-gradient-descent/>`__
  - `A Gentle Introduction to Mini-Batch Gradient Descent and How to Configure Batch Size <https://machinelearningmastery.com/gentle-introduction-mini-batch-gradient-descent-configure-batch-size/>`__


Distributed AI
==============

> Distributed AI is a computing paradigm that bypasses the need to move vast amounts of data and provides the ability to analyze data at the source. Gartner, a global provider of business insights, estimates that by 2025, 75 percent of data will be created and processed outside the traditional data center or cloud. This explosion of data being generated by people and machines from mobile devices, Internet of Things (IoTs), and machine data from production floors makes us rethink where computing needs to be performed.

Federated training
------------------

Federated learning is a machine learning technique that trains an algorithm across multiple decentralized edge devices or servers holding local data samples, without exchanging them.

References:

  - `What is Distributed AI? <https://developer.ibm.com/learningpaths/get-started-distributed-ai-apis/what-is-distributed-ai/>`__
  - `Federated Learning: Collaborative Machine Learning without Centralized Training Data <https://ai.googleblog.com/2017/04/federated-learning-collaborative.html>`__

Encoder-decoder models
======================

Encoder-Decoder models are a family of models which learn to map data-points from an input domain to an output domain via a two-stage network: The encoder, represented by an encoding function z = f(x), compresses the input into a latent-space representation; the decoder, y = g(z), aims to predict the output from the latent space representation.

Autoencoder
-----------

Autoencoders are special cases of encoder-decoder models in which the input and output are the same.

Recurrent neural network (RNN)
==============================

A recurrent neural network is a class of artificial neural networks where connections between nodes form a directed or undirected graph along a temporal sequence. This allows it to exhibit temporal dynamic behavior.

Long short-term memory (LSTM)
-----------------------------

In theory, classic (or "vanilla") RNNs can keep track of arbitrary long-term dependencies in the input sequences. The problem with vanilla RNNs is computational (or practical) in nature: when training a vanilla RNN using back-propagation, the long-term gradients which are back-propagated can "vanish" (that is, they can tend to zero) or "explode" (that is, they can tend to infinity), because of the computations involved in the process, which use finite-precision numbers. RNNs using LSTM units partially solve the vanishing gradient problem, because LSTM units allow gradients to also flow unchanged. However, LSTM networks can still suffer from the exploding gradient problem.

Attention
---------

Attention is a mechanism combined in the RNN allowing it to focus on certain parts of the input sequence when predicting a certain part of the output sequence, enabling easier learning and of higher quality. Combination of attention mechanisms enabled improved performance in many tasks making it an integral part of modern RNN networks.

References:

  - `Understanding LSTM Networks <https://colah.github.io/posts/2015-08-Understanding-LSTMs/>`__

BiLingual Evaluation Understudy (BLEU)
--------------------------------------

BLEU (BiLingual Evaluation Understudy) is a metric for automatically evaluating machine-translated text. The BLEU score is a number between zero and one that measures the similarity of the machine-translated text to a set of high quality reference translations. A value of 0 means that the machine-translated output has no overlap with the reference translation (low quality) while a value of 1 means there is perfect overlap with the reference translations (high quality).

Continuous bag of words (CBOW)
==============================

CBOW or Continous bag of words is to use embeddings in order to train a neural network where the context is represented by multiple words for a given target words.

Naive Bayes
===========

Naive Bayes methods are a set of supervised learning algorithms based on applying Bayes’ theorem with the “naive” assumption of conditional independence between every pair of features given the value of the class variable.

:math:`P(y \mid x_1, \dots, x_n) = \frac{P(y) P(x_1, \dots, x_n \mid y)} {P(x_1, \dots, x_n)}`

References:

  - `1.9. Naive Bayes <https://scikit-learn.org/stable/modules/naive_bayes.html#gaussian-naive-bayes>`__

Popular architectures
=====================

Inception
---------

Visual Geometry Group (VGG)
---------------------------

VGG is a standard deep CNN architecture, with the "deep" referring to the number of layers with VGG-16 or VGG-19 consisting of 16 and 19 convolutional layers.

.. image: https://viso.ai/wp-content/uploads/2021/10/vgg-neural-network-architecture.png

Residual Network (ResNet)
-------------------------

ResNet is a CNN architecture that overcame the "vanishing gradient" problem, making it possible to construct networks with up to thousands of convolutional layers, which outperform shallower networks.

NASNet
------

Xception
--------

[`ref <https://arxiv.org/pdf/1610.02357.pdf>`__]
