# AL Apr 11
import os

class paramInfo:

    def __init__(self, line=None):
        self.name = ''
        self.isDerived = False
        self.label = ''
        self.comment = ''
        self.filenameLoadedFrom = ''
        self.number = None
        if line is not None:
            self.setFromString(line)

    def setFromString(self, line):
        items = line.split(None, 1)
        self.name = items[0]
        if self.name.endswith('*'):
            self.name = self.name.strip('*')
            self.isDerived = True
        if len(items) > 1:
            tmp = items[1].split('#', 1)
            self.label = tmp[0].strip().replace('!', '\\')
            if len(tmp) > 1:
                self.comment = tmp[1].strip()
            else:
                self.comment = ''
        return self

    def setFromStringWithComment(self, items):
        self.setFromString(items[0])
        if items[1] != 'NULL':self.comment = items[1]

    def string(self, wantComments=True):
        res = self.name;
        if (self.isDerived): res = res + '*'
        res = res + '\t' + self.label
        if (wantComments and self.comment != ''):
            res = res + '\t#' + self.comment
        return res

class paramList:

    def __init__(self, fileName=None, setParamNameFile=None):

        self.names = []
        if fileName is not None: self.loadFromFile(fileName)
        if setParamNameFile is not None: self.setLabelsAndDerivedFromParamNames(setParamNameFile)

    def numDerived(self):
        return len([1 for info in self.names if info.isDerived])

    def list(self):
        return [name.name for name in self.names]

    def listString(self):
        return " ".join(self.list())

    def numParams(self):
        return len(self.names)

    def numNonDerived(self):
        return len([1 for info in self.names if not info.isDerived])

    def parWithNumber(self, num):
        for par in self.names:
            if par.number == num:
                return par
        return None

    def parWithName(self, name, error=False):
        for par in self.names:
            if par.name == name:
                return par
        if error: raise Exception("parameter name not found:" + name)
        return None

    def numberOfName(self, name):
        for i, par in enumerate(self.names):
            if par.name == name:return i
        return -1

    def parsWithNames(self, names, error=False):
        return [self.parWithName(name, error) for name in names]

    def setLabelsAndDerivedFromParamNames(self, fname):
        p = paramNames(fname)
        for par in p.names:
            param = self.parWithName(par.name)
            if not param is None:
                param.label = par.label
                param.isDerived = par.isDerived

    def fileList(self, fname):
        textFileHandle = open(fname)
        textFileLines = textFileHandle.readlines()
        textFileHandle.close()
        return textFileLines

    def deleteIndices(self, indices):
        self.names = [name for i, name in enumerate(self.names) if not i in indices]

    def filteredCopy(self, params):
        usedNames = self.__class__()
        for name in self.names:
            p = params.parWithName(name.name)
            if p is not None: usedNames.names.append(name)
        return usedNames

class paramNames(paramList):


    def loadFromFile(self, fileName):

        self.filenameLoadedFrom = os.path.split(fileName)[1]
        f = open(fileName)
        self.names = [paramInfo(line) for line in [s.strip() for s in f] if line != '']
        f.close()

    def loadFromKeyWords(self, keywordProvider):
        num_params_used = keywordProvider.keyWord_int('num_params_used')
        num_derived_params = keywordProvider.keyWord_int('num_derived_params')
        nparam = num_params_used + num_derived_params
        for i in range(nparam):
            info = paramInfo()
            info.setFromStringWithComment(keywordProvider.keyWordAndComment('param_' + str(i + 1)))
            self.names.append(info)
        return nparam

    def saveKeyWords(self, keywordProvider):
        keywordProvider.setKeyWord_int('num_params_used', len(self.names) - self.numDerived())
        keywordProvider.setKeyWord_int('num_derived_params', self.numDerived())
        for i, name in enumerate(self.names):
            keywordProvider.setKeyWord('param_' + str(i + 1), name.string(False).replace('\\', '!'), \
                name.comment)


    def saveAsText(self, fileName):
        textFileHandle = open(fileName, 'w')
        for info in self.names:
            textFileHandle.write(info.string() + '\n')
        textFileHandle.close()

