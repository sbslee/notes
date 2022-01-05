Statistics
**********

Frequently used commands for statistics
=======================================

* To calculate a 95% confidence interval from a 1D array:

    Here is how we can calculate the interval using the ``scipy.stats.t.interval`` method:

    .. code:: python3

        >>> mu, sigma, alpha = 6, 0.1, 0.95
        >>> s = np.random.normal(mu, sigma, 1000)
        >>> n, loc, scale = len(s), np.mean(s), stats.sem(s)
        >>> stats.t.interval(alpha=alpha, df=n-1, loc=loc, scale=scale)
        (5.998481753893311, 6.010448176896994)

    We can also calculate the interval manually:

    .. code:: python3

        >>> h = scale * stats.t.ppf((1 + alpha) / 2., n-1)
        >>> print((loc-h, loc+h))
        (5.998481746635707, 6.0104481841545985)

Statistical tests
=================

Student's t-test
----------------

This is a parametric test. Use the ``scipy.stats.ttest_ind`` method. Be careful with the ``equal_var`` option (i.e. Welch's t-test).

Paired samples t-test
---------------------

This is a parametric test for paired samples.

One-way ANOVA
-------------

This is the generalization of t-test.

Wilcoxon Rank-sum test
----------------------

This is a non-parametric test. It's also called Mann–Whitney U test. Here, you have two possibile methods: ``scipy.stats.mannwhitneyu`` and ``scipy.stats.ranksums``. Note that ``scipy.stats.ranksums`` should be used to compare two samples from continuous distributions. It does not handle ties between measurements in x and y. Therefore, I recommend using ``scipy.stats.mannwhitneyu`` for all cases because of its tie-handling and an optional continuity correction.

Wilcoxon Signed-rank test
-------------------------

This is a non-parametric test for paired samples. Use the ``scipy.stats.wilcoxon`` method.

Kruskal-Wallis test
-------------------

This is a non-parametric version of the one-way ANOVA test (i.e. one-way ANOVA on rank). Just as ANOVA can be applied to more than two groups and is a generalization of the t-test (which works with 2 groups only), Kruskal-Wallis can be applied to 2+ groups and is a generalization of the Mann-Whitney U test.

Note that the results of Kruskal-Wallis and Mann-Whitney U test may differ because 1) the ranks used for the Mann-Whitney U test are not the ranks used by the Kruskal-Wallis test; and 2) the rank sum tests do not use the pooled variance implied by the Kruskal-Wallis null hypothesis. Hence, it is not recommended to use Mann-whitney U test as a post hoc test after Kruskal-Wallis test.

Fisher's exact test
-------------------

Fisher's exact test is a statistical significance test used in the analysis of contingency tables. Although in practice it is employed when sample sizes are small, it is valid for all sample sizes. It is named after its inventor, Ronald Fisher, and is one of a class of exact tests, so called because the significance of the deviation from a null hypothesis (e.g., P-value) can be calculated exactly, rather than relying on an approximation that becomes exact in the limit as the sample size grows to infinity, as with many statistical tests.

Here's an example: A sample of teenagers might be divided into male and female on the one hand, and those that are and are not currently studying for a statistics exam on the other. We hypothesize, for example, that the proportion of studying individuals is higher among the women than among the men, and we want to test whether any difference of proportions that we observe is significant.

+--------------+-----+-------+-------+
|              | Men | Women | Total |
+==============+=====+=======+=======+
| Studying     | 1   | 9     | 10    |
+--------------+-----+-------+-------+
| Not-studying | 11  | 3     | 14    |
+--------------+-----+-------+-------+
| Total        | 12  | 12    | 24    |
+--------------+-----+-------+-------+

The probability of obtaining such set of values is given by:

:math:`p=\frac{\left ( \frac{a+b}{a} \right )\left ( \frac{c+d}{c} \right )}{\left ( \frac{n}{a+c} \right )}=\frac{\left ( \frac{a+b}{b} \right )\left ( \frac{c+d}{d} \right )}{\left ( \frac{n}{b+d} \right )}=\frac{(a+b)!(c+d)!(a+c)!(b+d)!}{a!b!c!d!n!}`

References:

  - `Fisher's exact test <https://en.wikipedia.org/wiki/Fisher%27s_exact_test>`__

Chi-squared test
----------------

A chi-squared test is a statistical hypothesis test that is valid to perform when the test statistic is chi-squared distributed under the null hypothesis, specifically Pearson's chi-squared test and variants thereof. Pearson's chi-squared test is used to determine whether there is a statistically significant difference between the expected frequencies and the observed frequencies in one or more categories of a contingency table.

Here's an example: Suppose there is a city of 1,000,000 residents with four neighborhoods: A, B, C, and D. A random sample of 650 residents of the city is taken and their occupation is recorded as "white collar", "blue collar", or "no collar". The null hypothesis is that each person's neighborhood of residence is independent of the person's occupational classification.

+--------------+-----+-----+-----+-----+-------+
|              | A   | B   | C   | D   | Total |
+==============+=====+=====+=====+=====+=======+
| White collar | 90  | 60  | 104 | 95  | 359   |
+--------------+-----+-----+-----+-----+-------+
| Blue collar  | 30  | 50  | 51  | 20  | 151   | 
+--------------+-----+-----+-----+-----+-------+
| No collar    | 30  | 40  | 45  | 35  | 150   |
+--------------+-----+-----+-----+-----+-------+
| Total        | 150 | 150 | 200 | 150 | 650   |
+--------------+-----+-----+-----+-----+-------+

In each cell of the table, we can calculate :math:`\frac{(observed-expected)^{2}}{expected}` and the sum of these quantities over all of the cells is the test statistic.

The chi-squared test applies an approximation assuming the sample is large, while the Fisher's exact test runs an exact procedure especially for small-sized samples.

References:

  - `Chi-squared test <https://en.wikipedia.org/wiki/Chi-squared_test#Fisher's_exact_test>`__

Levene's test
-------------

It is an inferential statistic used to assess the equality of variances for a variable calculated for two or more groups. Some common statistical procedures assume that variances of the populations from which different samples are drawn are equal. Levene's test assesses this assumption. It tests the null hypothesis that the population variances are equal (called homogeneity of variance or homoscedasticity). If the resulting p-value of Levene's test is less than some significance level (typically 0.05), the obtained differences in sample variances are unlikely to have occurred based on random sampling from a population with equal variances. Thus, the null hypothesis of equal variances is rejected and it is concluded that there is a difference between the variances in the population. Levene's test is implemented in the method ``scipy.stats.levene``. Some of the procedures typically assuming homoscedasticity, for which one can use Levene's tests, include analysis of variance and t-tests.

Shapiro-Wilk test
-----------------

The Shapiro-Wilk test (``scipy.stats.shapiro``) tests the null hypothesis that the data was drawn from a normal distribution.

Confusing concepts
==================

Confidence intervals vs. confidence levels
------------------------------------------

According to this `post <https://www.statisticshowto.com/probability-and-statistics/confidence-interval/>`__:

    Confidence levels are expressed as a percentage (for example, a 95% confidence level). It means that should you repeat an experiment or survey over and over again, 95 percent of the time your results will match the results you get from a population (in other words, your statistics would be sound!). Confidence intervals are your results and they are usually numbers. For example, you survey a group of pet owners to see how many cans of dog food they purchase a year. You test your statistic at the 99 percent confidence level and get a confidence interval of (200,300). That means you think they buy between 200 and 300 cans a year. You’re super confident (99% is a very high level!) that your results are sound, statistically.

Bootstrap vs. permutation tests
-------------------------------

According to this `post <http://pillowlab.princeton.edu/teaching/mathtools16/slides/lec21_Bootstrap.pdf>`__:

    Bootstrapping generally refers to statistical approach to quantifying uncertainty by re-using the data, specifically random resampling with replacement. Permutation-based analyses resemble the bootstrap in that they rely on randomizations of the observed data. The primary difference is that while bootstrap analyses typically seek to quantify the sampling distribution of some statistic computed from the data, permutation analyses typically seek to quantify the null distribution. That is, they seek to break whatever structure might be preset in a dataset, and quantify the kinds of patterns one expects to see “purely by chance.”

R-squared vs. adjusted R-squared
--------------------------------

According to the `website <https://www.investopedia.com/ask/answers/012615/whats-difference-between-rsquared-and-adjusted-rsquared.asp>`__:

    Adjusted R-squared is a modified version of R-squared that has been adjusted for the number of predictors in the model. The adjusted R-squared increases when the new term improves the model more than would be expected by chance. It decreases when a predictor improves the model by less than expected. Typically, the adjusted R-squared is positive, not negative. It is always lower than the R-squared.

    Adding more independent variables or predictors to a regression model tends to increase the R-squared value, which tempts makers of the model to add even more variables. This is called overfitting and can return an unwarranted high R-squared value. Adjusted R-squared is used to determine how reliable the correlation is and how much it is determined by the addition of independent variables.

Terminology
===========

+------------------------------------------------------------+-------------------------------------------------+
| Terminology                                                | Derivation                                      |
+============================================================+=================================================+
| true positive (TP)                                         |                                                 |
+------------------------------------------------------------+-------------------------------------------------+
| true negative (NP)                                         |                                                 |
+------------------------------------------------------------+-------------------------------------------------+
| false positive (FP)                                        |                                                 |
+------------------------------------------------------------+-------------------------------------------------+
| false negative (FN)                                        |                                                 |
+------------------------------------------------------------+-------------------------------------------------+
| sensitivity, recall, hit rate, or true positive rate (TPR) | :math:`TPR = TP / P = TP / (TP + FN) = 1 - FNR` |
+------------------------------------------------------------+-------------------------------------------------+
| specificity, selectivity or true negative rate (TNR)       | :math:`TNR = TN / N = TN / (TN + FP) = 1 - FPR` |
+------------------------------------------------------------+-------------------------------------------------+
| precision or positive predictive value (PPV)               | :math:`PPV = TP / (TP + FP) = 1 - FDR`          |
+------------------------------------------------------------+-------------------------------------------------+
| negative predictive value (NPV)                            | :math:`NPV = TN / (TN + FN) = 1 - FOR`          |
+------------------------------------------------------------+-------------------------------------------------+
| miss rate or false negative rate (FNR)                     | :math:`FNR = FN / P = FN / (FN + TP) = 1 - TPR` |
+------------------------------------------------------------+-------------------------------------------------+
| fall-out or false positive rate (FPR)                      | :math:`FPR = FP / N = FP / (FP + TN) = 1 - TNR` |
+------------------------------------------------------------+-------------------------------------------------+
| false discovery rate (FDR)                                 | :math:`FDR = FP / (FP + TP) = 1 - PPV`          |
+------------------------------------------------------------+-------------------------------------------------+
| false omission rate (FOR)                                  | :math:`FOR = FN / (FN + TN) = 1 - NPV`          |
+------------------------------------------------------------+-------------------------------------------------+
| accuracy (ACC)                                             | :math:`ACC = (TP + TN)/(TP + TN + FP + FN)`     |
+------------------------------------------------------------+-------------------------------------------------+

Common sampling distributions
=============================

+-----------------------+--------------------------------------+----------------------------------+
|                       | Draw with replacement                | Draw without replacement         |
|                       |                                      |                                  |
|                       | (probability of success is constant) | (probability of success changes) |
+=======================+======================================+==================================+
| Fixed number          | Binomial (Bernoulli is               | Hypergeometric                   |
|                       |                                      |                                  |
| of trials (:math:`n`) | special case when :math:`n=1`)       |                                  |
+-----------------------+--------------------------------------+----------------------------------+
| Draw until            | Negative Binomial (Geometric is      | Negative Hypergeometric          |
|                       |                                      |                                  |
| :math:`k` successes   | special case when :math:`k=1`)       |                                  |
+-----------------------+--------------------------------------+----------------------------------+

Poisson distribution
====================

The Poisson distribution is a discrete probability distribution that expresses the probability of a given number of events occurring in a fixed interval of time or space if these events occur with a known constant mean rate and independently of the time since the last event. The Poisson distribution can also be used for the number of events in other specified intervals such as distance, area or volume.

The Poisson distribution assumes that the mean and variance are the same. The negative binomial distribution has one parameter more than the Poisson regression that adjusts the variance independently from the mean. The Poisson distribution is a special case of the negative binomial distribution.

References:

  - `Difference between binomial, negative binomial and Poisson regression <https://stats.stackexchange.com/questions/60643/difference-between-binomial-negative-binomial-and-poisson-regression>`__

**Pro tip:** The RNAseq field uses negative binomial, the 16S microbiome field uses rarefying plus proportions, and the ChIP-seq field uses Poisson-based models.

Dirichlet distribution
======================

The Dirichlet distribution is a generalization of the Beta distribution for multiple random variables. It is over vectors whose values are all in the interval [0,1] and the sum of values in the vector is 1. In other words, the vectors in the sample space of the Dirichlet have the same properties as probability distribtutions. Therefore, the Dirichlet distribution can be thought of as a "distribution over distributions".

References:

  - `Continuous Distributions: Beta and Dirichlet Distributions <https://www.youtube.com/watch?v=CEVELIz4WXM>`__

Fisher's method
===============

According to the Wikipedia page:

"In statistics, Fisher's method, also known as Fisher's combined probability test, is a technique for data fusion or "meta-analysis" (analysis of analyses). It was developed by and named for Ronald Fisher. In its basic form, it is used to combine the results from several independence tests bearing upon the same overall hypothesis (H0)."

"Under Fisher's method, two small p-values P1 and P2 combine to form a smaller p-value. The yellow-green boundary defines the region where the meta-analysis p-value is below 0.05. For example, if both p-values are around 0.10, or if one is around 0.04 and one is around 0.25, the meta-analysis p-value is around 0.05."

ROC curve and PR curve
======================

The receiver operating characteristic (ROC) curve is created by plotting the true positive rate (TPR) against the false positive rate (FPR) at various threshold settings. The precision-recall (PR) curve shows the tradeoff between precision (PPV) and recall (which is equivalent to TPR) for different threshold. Therefore, both ROC curves and PR curves share the TPR term. According to this `CV post <https://stats.stackexchange.com/questions/7207/roc-vs-precision-and-recall-curves>`__, the key difference between the two is that:

> ROC curves will be the same no matter what the baseline probability is, but PR curves may be more useful in practice for needle-in-haystack type problems or problems where the "positive" class is more interesting than the negative class.
