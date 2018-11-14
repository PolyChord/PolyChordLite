import numpy
from scipy.special import erfinv


class UniformPrior:
    def __init__(self, a, b):
        self.a = a
        self.b = b

    def __call__(self, x):
        return self.a + (self.b-self.a) * x


class GaussianPrior:
    def __init__(self, mu, sigma):
        self.mu = mu
        self.sigma = sigma

    def __call__(self, x):
        return self.mu + self.sigma * numpy.sqrt(2) * erfinv(2*x-1)


class LogUniformPrior(UniformPrior):

    def __call__(self, x):
        return self.a * (self.b/self.a) ** x


def forced_indentifiability_transform(x):
    N = len(x)
    t = numpy.zeros(N)
    t[N-1] = x[N-1]**(1./N)
    for n in range(N-2, -1, -1):
        t[n] = x[n]**(1./(n+1)) * t[n+1]
    return t


class SortedUniformPrior(UniformPrior):
    def __call__(self, x):
        t = forced_indentifiability_transform(x)
        return super(SortedUniformPrior, self).__call__(t)


class LogSortedUniformPrior(LogUniformPrior):
    def __call__(self, x):
        t = forced_indentifiability_transform(x)
        return super(LogSortedUniformPrior, self).__call__(t)
