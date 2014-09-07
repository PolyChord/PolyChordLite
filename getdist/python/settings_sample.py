# sample settings for a particular grid run

newCovmats = True
start_at_bestfit = False

# dataset names
planck = 'planck'
lowl = 'lowl'
lowLike = 'lowLike'
lensing = 'lensing'
highL = 'highL'
BAO = 'BAO'
HST = 'HST'
SNLS = 'SNLS'
Union = 'Union2'

# set up groups of parameters and data sets
class group:pass

g = group()
g.params = [[], ['mnu'], ['nnu'], ['nrun']]

g.datasets = []

# lists of dataset names to combine, with corresponding sets of inis to include
g.datasets.append([[planck, lowl, lowLike], ['planck.ini', 'lowl.ini', 'lowLike.ini']])
g.datasets.append([[planck, lowl, lowLike, highL], ['planck.ini', 'lowl.ini', 'lowLike.ini']])


# add importance name tags, and list of specific .ini files to include (in batch1/)
g.importanceRuns = []
g.importanceRuns.append([[BAO], ['BAO.ini']])

groups = [g]

# try to match run to exisitng covmat
covrenames = [['_my_new_data', '']]

ini_dir = 'batch1/'
# ini files you want to base each set of runs on
defaults = ['common_batch1.ini']

importanceDefaults = ['importance_sampling.ini']
