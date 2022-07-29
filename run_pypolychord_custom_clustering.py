    
from numpy import pi, log, sqrt
import numpy as np
import pypolychord
from pypolychord.settings import PolyChordSettings
from pypolychord.priors import UniformPrior
from scipy.spatial import distance_matrix
from scipy.special import logsumexp
try:
    from mpi4py import MPI
except ImportError:
    pass



#| Define a four-dimensional twin-spherical gaussian likelihood,
#| width sigma=0.1, centered on ±[0.5]^nDims with one derived parameter.
#| The derived parameter is the squared radius from [0]^nDims

nDims = 4
nDerived = 1
sigma = 0.1

centres = np.array([[-0.5] * nDims, [0.5] * nDims])
weights = np.array([0.5, 0.5])


def likelihood(theta):
    """Twin Gaussian likelihood."""
    nDims = len(theta)
    logL = -np.log(2 * np.pi * sigma * sigma) * nDims / 2 
    logL += logsumexp(-np.sum((theta - centres) ** 2, axis=-1) / 2 / sigma / sigma, b = weights)
    return logL, [np.sum(theta**2)]


#| Define a box uniform prior from -1 to 1

def prior(hypercube):
    """ Uniform prior from [-1,1]^D. """
    return UniformPrior(-1, 1)(hypercube)


#| Optional cluster function allow user-defined clustering
## KNN clustering algorithm translated from clustering.F90

def compute_nn(position_matrix):
    return np.argsort(distance_matrix(position_matrix, position_matrix))


def relabel(labels):
    k = max(labels)
    appearance_order = {}
    num_found = 0

    for label in labels:
        if label not in appearance_order:
            appearance_order[label] = num_found
            num_found += 1

    for i, old_label in enumerate(labels):
        labels[i] = appearance_order[old_label]
    return labels


def do_clustering(knn_array):
    num_points = knn_array.shape[0]
    labels = np.arange(num_points)
    for iii in range(num_points):
        for ii in range(iii + 1, num_points):
            if labels[ii] != labels[iii]:
                if (ii in knn_array[iii]) or (iii in knn_array[ii]):
                    for i in range(num_points):
                        if labels[i] == labels[ii] or labels[i] == labels[iii]:
                            labels[i] = min([labels[ii], labels[iii]])
    return relabel(labels)


def knn_cluster(position_matrix):
    """slight concern if two points are the same because sklearn"""
    npoints = position_matrix.shape[0]
    k = min(npoints, 10)
    nn_array = compute_nn(position_matrix)
    labels = np.arange(npoints)
    num_clusters = npoints

    labels_old = labels

    n = 2
    while n <= k:
        labels = do_clustering(nn_array[:, :n])
        num_clusters = max(labels) + 1

        if num_clusters <= 0:
            raise ValueError("somehow got <= 0 clusters")
        elif 1 == num_clusters:
            return labels
        elif np.all(labels == labels_old):
            break
        elif k == n:
            k = min(k * 2, npoints)

        labels_old = labels
        n += 1

    if num_clusters > 1:
        i_cluster = 0
        while i_cluster < num_clusters:
            cluster = position_matrix[labels == i_cluster]
            labels[labels == i_cluster] = knn_cluster(cluster) + num_clusters
            labels = relabel(labels)
            if num_clusters - 1 == max(labels):
                i_cluster += 1
    return labels


#| Initialise the settings

settings = PolyChordSettings(nDims, nDerived)
settings.file_root = "custom_clustering"
settings.nlive = 200
settings.do_clustering = True
settings.read_resume = False

#| Run PolyChord

output = pypolychord.run_polychord(likelihood, nDims, nDerived, settings, prior, cluster=knn_cluster)

#| Create a paramnames file

paramnames = [('p%i' % i, r'\theta_%i' % i) for i in range(nDims)]
paramnames += [('r*', 'r')]
output.make_paramnames_files(paramnames)

#| Make an anesthetic plot (could also use getdist)
try:
    from anesthetic import NestedSamples
    samples = NestedSamples(root= settings.base_dir + '/' + settings.file_root)
    fig, axes = samples.plot_2d(['p0','p1','p2','p3','r'])
    fig.savefig('posterior.pdf')

except ImportError:
    try:
        import getdist.plots
        posterior = output.posterior
        g = getdist.plots.getSubplotPlotter()
        g.triangle_plot(posterior, filled=True)
        g.export('custom_clustering.pdf')
    except ImportError:
        print("Install matplotlib and getdist for plotting examples")

    print("Install anesthetic or getdist for for plotting examples")
