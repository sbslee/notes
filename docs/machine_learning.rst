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

Rectified Linear Unit (ReLU)
============================

The rectified linear activation function or ReLU for short is a piecewise linear function that will output the input directly if it is positive, otherwise, it will output zero. It has become the default activation function for many types of neural networks because a model that uses it is easier to train and often achieves better performance. For example, the sigmoid and hyperbolic tangent activation functions cannot be used in networks with many layers due to the vanishing gradient problem. The rectified linear activation function overcomes the vanishing gradient problem, allowing models to learn faster and perform better. Therefore, the rectified linear activation is the default activation when developing multilayer perceptron and convolutional neural networks.

References:

  - `A Gentle Introduction to the Rectified Linear Unit (ReLU) <https://machinelearningmastery.com/rectified-linear-activation-function-for-deep-learning-neural-networks/>`__

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
