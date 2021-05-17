Statistics
**********

Statistical tests
=================

Student's t-test
----------------

This is a parametric test. Use the ``scipy.stats.ttest_ind`` method. Be careful with the ``equal_var`` option (i.e. Welch's t-test).

Paired Samples t-test
---------------------

This is a parametric test for paired samples.

One-way ANOVA
-------------

This is the generalization of t-test.

Wilcoxon Rank-sum test
----------------------

This is a non-parametric test. It's also called Wilcoxon rank-sum test. Here, you have two possibile methods: ``scipy.stats.mannwhitneyu`` and ``scipy.stats.ranksums``. I recommend the former.

Wilcoxon Signed-rank test
-------------------------

This is a non-parametric test for paired samples. Use the ``scipy.stats.wilcoxon`` method.

Kruskal-Wallis test
-------------------

This is a non-parametric version of the one-way ANOVA test (i.e. one-way ANOVA on rank). Just as ANOVA can be applied to more than two groups and is a generalization of the t-test (which works with 2 groups only), Kruskal-Wallis can be applied to 2+ groups and is a generalization of the Mann-Whitney U test.

Note that the results of Kruskal-Wallis and Mann-Whitney U test may differ because 1) the ranks used for the Mann-Whitney U test are not the ranks used by the Kruskal-Wallis test; and 2) the rank sum tests do not use the pooled variance implied by the Kruskal-Wallis null hypothesis. Hence, it is not recommended to use Mann-whitney U test as a post hoc test after Kruskal-Wallis test.

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


:math:`{H}'=-\sum_{i=1}^{R}p_i \ln p_i`
