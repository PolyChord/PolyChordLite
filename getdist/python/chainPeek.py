# Just look at parameters being produced in single chain

import os, paramNames, sys
import numpy as np


if len(sys.argv) < 3:
    print 'Usage: python/chainPeek.py chains/chainroot paramname [chain_index]'
    print 'e.g. python/makeGrid.py chains/test omegabh2 1'
    sys.exit()


chain = os.path.abspath(sys.argv[1])
param = sys.argv[2]
index = 1
if len(sys.argv) > 4: index = int(sys.argv[3])

d = np.loadtxt(chain + '_' + str(index) + '.txt')
names = paramNames.paramNames(chain + '.paramnames')

print d[:, names.numberOfName(param) + 2]
