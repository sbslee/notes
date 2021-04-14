# Statistics

## Table of Contents

* [Terminology](#Terminology)
* [Student's t-test](#Student's-t-test)
* [Paired Samples t-test](#Paired-Samples-t-test)
* [One-way ANOVA](#One-way-ANOVA)
* [Wilcoxon Rank-sum Test / Mann–Whitney U Test](#Wilcoxon-Rank-sum-Test-/-Mann–Whitney-U-Test)
* [Wilcoxon Signed-rank Test](#Wilcoxon-Signed-rank-Test)
* [Kruskal-Wallis Test / One-way ANOVA on Rank](#Kruskal-Wallis-Test-/-One-way-ANOVA-on-Rank)
* [Common sampling distributions](#Common-sampling-distributions)
* [Poisson distribution](#Poisson-distribution)
* [Dirichlet distribution](#Dirichlet-distribution)
* [Fisher's method](#Fisher's-method)
* [ROC curve and PR curve](#ROC-curve-and-PR-curve)

## Terminology <a name="Terminology"></a>

| Terminology                                                | Derivation                              |
| ---------------------------------------------------------- | --------------------------------------  |
| true positive (TP)                                         |                                         |
| true negative (NP)                                         |                                         |
| false positive (FP)                                        |                                         |
| false negative (FN)                                        |                                         |
| sensitivity, recall, hit rate, or true positive rate (TPR) | TPR = TP / P = TP / (TP + FN) = 1 - FNR |
| specificity, selectivity or true negative rate (TNR)       | TNR = TN / N = TN / (TN + FP) = 1 - FPR |
| precision or positive predictive value (PPV)               | PPV = TP / (TP + FP) = 1 - FDR          |
| negative predictive value (NPV)                            | NPV = TN / (TN + FN) = 1 - FOR          |
| miss rate or false negative rate (FNR)                     | FNR = FN / P = FN / (FN + TP) = 1 - TPR |
| fall-out or false positive rate (FPR)                      | FPR = FP / N = FP / (FP + TN) = 1 - TNR |
| false discovery rate (FDR)                                 | FDR = FP / (FP + TP) = 1 - PPV          |
| false omission rate (FOR)                                  | FOR = FN / (FN + TN) = 1 - NPV          |

## Student's t-test <a name="Student's-t-test"></a>

This is a parametric test. Use the `scipy.stats.ttest_ind()` method. Be careful with the `equal_var` option (i.e. Welch's t-test).

## Paired Samples t-test <a name="Paired-Samples-t-test"></a>

This is a parametric test for paired samples.

## One-way ANOVA <a name="One-way-ANOVA"></a>

This is the generalization of t-test.

## Wilcoxon Rank-sum Test / Mann–Whitney U Test <a name="Wilcoxon-Rank-sum-Test-/-Mann–Whitney-U-Test"></a>

This is a non-parametric test. Here, you have two possibile methods: `scipy.stats.mannwhitneyu` and `scipy.stats.ranksums`. I recommend the former.

## Wilcoxon Signed-rank Test <a name="Wilcoxon-Signed-rank-Test"></a>

This is a non-parametric test for paired samples. Use the `scipy.stats.wilcoxon()` method.

## Kruskal-Wallis Test / One-way ANOVA on Rank <a name="Kruskal-Wallis-Test-/-One-way-ANOVA-on-Rank"></a>

This is a non-parametric version of the one-way ANOVA test. Just as ANOVA can be applied to more than two groups and is a generalization of the t-test (which works with 2 groups only), Kruskal-Wallis can be applied to 2+ groups and is a generalization of the Mann-Whitney U test.

Note that the results of Kruskal-Wallis and Mann-Whitney U test may differ because 1) the ranks used for the Mann-Whitney U test are not the ranks used by the Kruskal-Wallis test; and 2) the rank sum tests do not use the pooled variance implied by the Kruskal-Wallis null hypothesis. Hence, it is not recommended to use Mann-whitney U test as a post hoc test after Kruskal-Wallis test.

## Common sampling distributions <a name="Common-sampling-distributions"></a>

|                                | Draw with replacement<br>(probability of success is constant) | Draw without replacement<br>(probability of success changes) |
| ------------------------------ | ------------------------------------------------------------- | ------------------------------------------------------------ |
| Fixed number of trials (*n*)   | Binomial<br>(Bernoulli is special case when *n* = 1)          | Hypergeometric                                               |
| Draw until *k* successes       | Negative Binomial<br>(Geometric is special case when *k* = 1) | Negative Hypergeometric                                      |

## Poisson distribution <a name="Poisson-distribution"></a>

The Poisson distribution is a discrete probability distribution that expresses the probability of a given number of events occurring in a fixed interval of time or space if these events occur with a known constant mean rate and independently of the time since the last event. The Poisson distribution can also be used for the number of events in other specified intervals such as distance, area or volume.

The Poisson distribution assumes that the mean and variance are the same. The negative binomial distribution has one parameter more than the Poisson regression that adjusts the variance independently from the mean. The Poisson distribution is a special case of the negative binomial distribution.

Related posts:

* [Difference between binomial, negative binomial and Poisson regression](https://stats.stackexchange.com/questions/60643/difference-between-binomial-negative-binomial-and-poisson-regression)

**Pro tip:** The RNAseq field uses negative binomial, the 16S microbiome field uses rarefying plus proportions, and the ChIP-seq field uses Poisson-based models.

## Dirichlet distribution <a name="Dirichlet-distribution"></a>

The Dirichlet distribution is a generalization of the Beta distribution for multiple random variables. It is over vectors whose values are all in the interval [0,1] and the sum of values in the vector is 1. In other words, the vectors in the sample space of the Dirichlet have the same properties as probability distribtutions. Therefore, the Dirichlet distribution can be thought of as a "distribution over distributions".

Related posts:

* [Continuous Distributions: Beta and Dirichlet Distributions](https://www.youtube.com/watch?v=CEVELIz4WXM)

## Fisher's method <a name="Fisher's-method"></a>

According to the Wikipedia page:

"In statistics, Fisher's method, also known as Fisher's combined probability test, is a technique for data fusion or "meta-analysis" (analysis of analyses). It was developed by and named for Ronald Fisher. In its basic form, it is used to combine the results from several independence tests bearing upon the same overall hypothesis (H0)."

"Under Fisher's method, two small p-values P1 and P2 combine to form a smaller p-value. The yellow-green boundary defines the region where the meta-analysis p-value is below 0.05. For example, if both p-values are around 0.10, or if one is around 0.04 and one is around 0.25, the meta-analysis p-value is around 0.05."

## ROC curve and PR curve <a name="ROC-curve-and-PR-curve"></a>

The receiver operating characteristic (ROC) curve is created by plotting the true positive rate (TPR) against the false positive rate (FPR) at various threshold settings. The precision-recall (PR) curve shows the tradeoff between precision (PPV) and recall (which is equivalent to TPR) for different threshold. Therefore, both ROC curves and PR curves share the TPR term. According to this [CV post](https://stats.stackexchange.com/questions/7207/roc-vs-precision-and-recall-curves), the key difference between the two is that:

> ROC curves will be the same no matter what the baseline probability is, but PR curves may be more useful in practice for needle-in-haystack type problems or problems where the "positive" class is more interesting than the negative class.
