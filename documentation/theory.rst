Theory
======

In many disciplines, there is often a need to adapt a statistical model to existing data to be able to make statements regarding uncertain future outcomes. In particular, when assessing risks, an estimate of major losses must be based on events that, despite having a low probability of occurrence, have a high impact. Since the actual distribution of data -- the parent distribution -- is generally unknown, statisticians can fit a generalized Pareto distribution (GPD) to the data belonging to the tail.

For a very large class of parent distribution functions, the generalized Pareto distribution (GPD) can be used as a model for the tail. A certain threshold divides the parent distribution into two areas: a body and a tail region. Above the threshold, analyses are performed using the GPD as a model for the tail.

This FTT procedure (Find-The-Tail) provides a suitable and efficient method for determining the optimal threshold.

FTT uses a specially weighted MSE teststatistic :math:`AU^2` for determining the optimal value. The weights are chosen, such that the deviations (between EDF and the fitted GPD) at high quantiles are stronger penalized than at lower quantiles. This  upper tail teststatistic is evaluated successive over a sorted timeseries.

Teststatistic:

.. math::
    AU_{n}^2 = \frac{1}{2} n - \sum_{i=1}^n\left[ 2 F(x_{(i)}) + \frac{2(n-i)+1}{n} \ln\left(1-F(x_{(i)})\right) \right]

With :math:`n` is the sampe length, :math:`x_{(i)}` is the descending sorted data and :math:`F` is the fitted GPD.

For further information see Reference_.

.. _Reference: https://arxiv.org/abs/1805.10040