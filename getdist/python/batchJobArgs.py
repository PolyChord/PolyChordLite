import sys, batchJob, fnmatch
try: import argparse
except:
    print 'use "module load" to load python 2.7'
    sys.exit()

def argParser(desc=''):
    return argparse.ArgumentParser(description=desc)

class batchArgs():

        def __init__(self, desc='', importance=True, noBatchPath=False, notExist=False, converge=False):
            self.parser = argparse.ArgumentParser(description=desc)
            if not noBatchPath: self.parser.add_argument('batchPath', help='directory containing the grid')
            if converge: self.parser.add_argument('--converge', type=float, default=0, help='minimum R-1 convergence')
            self.importanceParameter = importance;
            self.notExist = notExist

        def parseForBatch(self):
            if self.importanceParameter:
                self.parser.add_argument('--noimportance', action='store_true', help='original chains only, no importance sampled')
                self.parser.add_argument('--importance', nargs='*', default=None, help='tags for importance sampling runs to include')
            self.parser.add_argument('--name', default=None, nargs='+', help='specific chain full name only (base_paramx_data1_data2)')
            self.parser.add_argument('--param', default=None, nargs='+', help='runs including specific parameter only (paramx)')
            self.parser.add_argument('--paramtag', default=None, help='runs with specific parameter tag only (base_paramx)')
            self.parser.add_argument('--data', default=None, help='runs including specific data only (data1)')
            self.parser.add_argument('--datatag', default=None, help='runs with specific data tag only (data1_data2)')
            self.parser.add_argument('--skip_data', default=None, help='skip runs containing specific data (data1)')
            self.parser.add_argument('--skip_param', default=None, help='skip runs containing specific parameter (paramx)')
            self.parser.add_argument('--group', default=None, nargs='+', help='include only runs with given group names')

            if self.notExist: self.parser.add_argument('--notexist', action='store_true')

            self.args = self.parser.parse_args()
            self.batch = batchJob.readobject(self.args.batchPath)
            return (self.batch, self.args)

        def wantImportance(self, importanceTag):
            return self.args.importance is None or len(self.args.importance) == 0 or importanceTag in self.args.importance

        def jobItemWanted(self, jobItem):
            return not jobItem.isImportanceJob and (self.args.importance is None) or jobItem.isImportanceJob and self.wantImportance(jobItem.importanceTag)

        def nameMatches(self, jobItem):
            if self.args.name is None: return True
            for pat in self.args.name:
                if fnmatch.fnmatch(jobItem.name, pat): return True
            return False

        def groupMatches(self, jobItem):
            return self.args.group is None or jobItem.group in self.args.group

        def dataMatches(self, jobItem):
            if self.args.datatag is None:
                if self.args.data is None:
                    return self.args.skip_data is None or not self.args.skip_data in jobItem.data_set.names
                return self.args.data in jobItem.data_set.names
            else:
                return jobItem.datatag == self.args.datatag

        def paramsMatch(self, jobItem):
            if self.args.paramtag is None:
                if self.args.param is None:
                    return self.args.skip_param is None or not self.args.skip_param in jobItem.param_set
                for pat in self.args.param:
                    if pat in jobItem.param_set: return self.args.skip_param is None or not self.args.skip_param in jobItem.param_set
                return False
            else:
                return jobItem.paramtag == self.args.paramtag

        def filteredBatchItems(self, wantSubItems=True, chainExist=False):
            for jobItem in self.batch.items(wantImportance=not self.args.noimportance, wantSubItems=wantSubItems):
                if (not chainExist or jobItem.chainExists()) and (self.jobItemWanted(jobItem) and self.nameMatches(jobItem) and self.paramsMatch(jobItem)  and self.dataMatches(jobItem)
                    and self.groupMatches(jobItem)): yield(jobItem)

        def sortedParamtagDict(self, chainExist=True):
            items = dict()
            for jobItem in self.filteredBatchItems():
                if not chainExist or jobItem.chainExists():
                    if not jobItem.paramtag in items: items[jobItem.paramtag] = []
                    items[jobItem.paramtag].append(jobItem)
            return sorted(items.iteritems())

        def filterForDataCompare(self, batch, datatags):
            items = []
            for tag, data in zip([self.batch.normalizeDataTag(data) for data in datatags], datatags):
                items += [jobItem for jobItem in batch if jobItem.datatag == data or jobItem.normed_data == tag]
            return items


