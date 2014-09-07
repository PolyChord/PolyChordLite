import numpy as np

class covMat():

    def __init__(self, filename=''):

        self.matrix = []
        self.paramNames = []
        self.size = 0
        if filename != '': self.loadFromFile(filename)

    def paramNameString(self):
        return " ".join(self.paramNames)

    def loadFromFile(self, filename):
        textFileHandle = open(filename)
        textFileLines = textFileHandle.readlines()
        textFileHandle.close()
        first = textFileLines[0].strip()
        if first.startswith('#'):
            paramNames = first[1:].split()
            self.size = len(paramNames)
        else:
            raise Exception('.covmat must now have parameter names header')
        matrix = [[0 for col in range(self.size)] for row in range(self.size)]
        used = []
        for i in range(self.size):
            splitLine = textFileLines[i + 1].split()
            for j in range(len(splitLine)):
                matrix[i][j] = float(splitLine[j])
            if matrix[i].count(0.) != self.size:
                used.append(i)

        self.size = len(used)
        self.matrix = np.empty((self.size, self.size))
        self.paramNames = []
        for i in range(self.size):
            self.paramNames.append(paramNames[used[i]])
            for j in range(self.size):
                self.matrix[i, j] = matrix[used[i]][used[j]]


    def saveToFile(self, filename):
        fout = open(filename, 'w')
        fout.write('# ' + self.paramNameString() + '\n')
        np.savetxt(fout, self.matrix, '%E')
        fout.close

    def rescaleParameter(self, name, scale):
        if name in self.paramNames:
            i = self.paramNames.index(name)
            self.matrix[:, i] = self.matrix[:, i] * scale;
            self.matrix[i, :] = self.matrix[i, :] * scale;
        else: print 'Not in covmat: ' + name

    def mergeCovmatWhereNew(self, cov2):
        params1 = self.paramNames
        params2 = cov2.paramNames

        C = covMat()
        C.paramNames.extend(params1)

        for param in cov2.paramNames:
                if param not in C.paramNames:
                    C.paramNames.append(param)
        l1 = len(params1)
        l2 = len(params2)
        l = len(C.paramNames)

        map1 = dict(zip(params1, range(0, l1)))
        map2 = dict(zip(params2, range(0, l2)))
        covmap = dict(zip(range(0, l), C.paramNames))

        C.matrix = np.zeros((l, l))
        for i in range(0, l):
            for j in range(0, l):
                if C.paramNames[i] in params1 and C.paramNames[j] in params1:
                    C.matrix[i, j] = self.matrix[map1[covmap[i]], map1[covmap[j]]]
                elif C.paramNames[i] in params2 and C.paramNames[j] in params2:
                    C.matrix[i, j] = cov2.matrix[map2[covmap[i]], map2[covmap[j]]]

        return C
