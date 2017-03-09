import numpy 


class UniformPrior:
    def __init__(self, a, b):
        self.a = a
        self.b = b

    def __call__(self, x):
        return self.a + (self.b-self.a) * x


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
