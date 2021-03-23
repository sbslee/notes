# Statistics

## Table of Contents

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

## Dirichlet distribution <a name="Dirichlet-distribution"></a>

The Dirichlet distribution is a generalization of the Beta distribution for multiple random variables. It is over vectors whose values are all in the interval [0,1] and the sum of values in the vector is 1. In other words, the vectors in the sample space of the Dirichlet have the same properties as probability distribtutions. Therefore, the Dirichlet distribution can be thought of as a "distribution over distributions".

Related posts:

* [Continuous Distributions: Beta and Dirichlet Distributions](https://www.youtube.com/watch?v=CEVELIz4WXM)

## Fisher's method <a name="Fisher's-method"></a>

According to the Wikipedia page:

"In statistics, Fisher's method, also known as Fisher's combined probability test, is a technique for data fusion or "meta-analysis" (analysis of analyses). It was developed by and named for Ronald Fisher. In its basic form, it is used to combine the results from several independence tests bearing upon the same overall hypothesis (H0)."

"Under Fisher's method, two small p-values P1 and P2 combine to form a smaller p-value. The yellow-green boundary defines the region where the meta-analysis p-value is below 0.05. For example, if both p-values are around 0.10, or if one is around 0.04 and one is around 0.25, the meta-analysis p-value is around 0.05."
