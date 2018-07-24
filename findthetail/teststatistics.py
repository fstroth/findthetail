import numpy as np
# TODO: Add references to the statistics


def anderson_darling(x):
    return -x.size-1/x.size*sum((2*(np.arange(x.size)+1)-1)*(np.log(x)+np.log(1-x[::-1])))


def cramer_von_mises(x):
    return 1/(12*x.size) + sum(((np.arange(x.size)+1-0.5)/x.size-x)**2)


def au2(x):
    return x.size/2 - sum(2*x+(2*(x.size-(np.arange(x.size)+1))+1)/x.size * np.log(1-x))
