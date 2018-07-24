How it works
============

The FTT algorithm consists of the following 9 steps.

    1. The loss data are sorted in descending order :math:`x_1 \geq x_2 \geq \ldots \geq x_n` .
    2. Starting with the largest losses :math:`x_1, x_2` (k = 2) the generalized Pareto distribution :math:`\hat{F}(x)` (= Estimation) is adapted to the data.
    3. With this distribution the measure of deviation :math:`AU^2_k` is calculated.
    4. Then the next smaller loss is added :math:`(k \rightarrow k + 1)` and the parameters of the generalized Pareto distribution are re-estimated.
    5. It continues with point (3) until the last loss value has been processed :math:`(k = n)`.
    6. Depending on :math:`k`, a time series of the deviation measure results: :math:`AU^2_k`.
    7. The minimal deviation at a particular :math:`k^* \in [1, n]` indicates the best fit of a model for the tail to the given loss data and the associated :math:`x_{k^*}` corresponds to the sought threshold u.
    8. To assess the result, the confidence level is determined. (This describes how large the likelihood of a wrong decision would be, if the adapted distribution and thus the decision for the threshold were rejected)
    9. For further assurance, the statistics of the standard goodness-of-fit test (CM and AD test) are evaluated.

.. note::
    1. FTT strictly separates the procedure for detecting the threshold and the goodness of fit to evaluate the quality of the fit.
    2. There is no need for any external parameters, all information are gained form the data alone.

