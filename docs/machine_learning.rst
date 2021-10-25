Machine learning
****************

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
