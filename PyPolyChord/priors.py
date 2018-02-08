import numpy 
import scipy.special


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
        return self.mu + self.sigma * numpy.sqrt(2) * scipy.special.erfinv(2*x-1)


class LogUniformPrior(UniformPrior):

    def __call__(self, x):
        return self.a * (self.b/self.a) ** x


class SortedUniformPrior(UniformPrior):

    def __call__(self, x):
        N = len(x)
        t = numpy.zeros(N)
        t[N-1] = x[N-1]**(1./N)
        for n in range(N-2,-1,-1):
            t[n] = x[n]**(1./(n+1)) * t[n+1]
        return UniformPrior.__call__(self, t)
