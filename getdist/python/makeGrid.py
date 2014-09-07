import os, sys, batchJob, iniFile


if len(sys.argv) < 2:
    print 'Usage: python/makeGrid.py new_directory_for_outputs grid_settings_python_file'
    print 'e.g. python/makeGrid.py /scratch/../testgrid settings_testgrid'
    sys.exit()


batchPath = os.path.abspath(sys.argv[1]) + os.sep

# 0: chains, 1: importance sampling, 2: best-fit, 3: best-fit and Hessian
cosmomcAction = 0

settings = __import__(sys.argv[2])

batch = batchJob.batchJob(batchPath, settings.ini_dir)

# priors and widths for parameters which are varied
if not hasattr(settings, 'params'):
    params = dict()
    params['mnu'] = '0.02 0 5 0.1 0.03'
    params['omegak'] = '-0.0008 -0.3 0.3 0.001 0.001'  # starting exactly on flat seems to confuse minimizer
    params['w'] = '-0.995 -3 -0.3 0.02 0.02'
    params['nnu'] = '3.046 0.05 10 0.05 0.05'
    params['nrun'] = '0 -1 1 0.005 0.001'
    params['nrunrun'] = '0 -1 1 0.005 0.001'    
    params['r'] = '0 0 3 0.03 0.03'
    params['Alens'] = '1 0 10 0.05 0.05'
    params['yhe'] = '0.245 0.1 0.5 0.006 0.006'
    params['alpha1'] = '0 -1 1 0.0003 0.0003'
    params['deltazrei'] = '0.5 0.1 3 0.3 0.3'
    params['wa'] = '0 -2 2 0.3 0.3'
    params['meffsterile'] = '0.1 0 3 0.1 0.03'
    params['Aphiphi'] = '1 0 10 0.02 0.02'
    params['nt'] = '0 -3 3 0.2 0.02'
    settings.params = params


if hasattr(settings, 'skip'): batch.skip = settings.skip
batch.makeItems(settings.groups)
batch.makeDirectories()
batch.save()

def setMinimize(jobItem, ini):
    ini.params['action'] = 2
    ini.params['lmin_store_all_cmb'] = 2000
    if 'omegak' in jobItem.param_set: ini.params['accuracy_level'] = 1.3

def updateIniParams(ini, params, path):
        for iniitem in params:
            if isinstance(iniitem, dict): ini.params.update(iniitem)
            elif isinstance(iniitem, basestring): ini.defaults.append(path + iniitem)
            elif isinstance(iniitem, (list, tuple)): updateIniParams(ini, iniitem, path)
            else: raise Exception('Unknown item in setting .ini/param list')


for jobItem in batch.items(wantSubItems=False):

        jobItem.makeChainPath()
        ini = iniFile.iniFile()

        for param in jobItem.param_set:
            ini.params['param[' + param + ']'] = settings.params[param]

        if 'mnu' in jobItem.param_set:
            ini.params['num_massive_neutrinos'] = 3
        if 'meffsterile' in jobItem.param_set:
            ini.params['param[mnu]'] = '0.06 0.06 0.06 0 0'
            ini.params['param[nnu]'] = '3.1 3.046 10 0.05 0.05'
            ini.params['num_massive_neutrinos'] = 1
            ini.params['accuracy_level'] = 1.2  # to use 4 rather than 3 momentum modes
        if 'yhe' in jobItem.param_set:
            ini.params['bbn_consistency'] = False
        if 'r' in jobItem.param_set:
            ini.params['compute_tensors'] = True
        if 'nt' in jobItem.param_set:
            ini.params['inflation_consistency'] = False
            ini.params['lmax_tensor'] = 1000
#            ini.params['pivot_k'] = 0.002
        if hasattr(settings, 'extra_opts'):
            ini.params.update(settings.extra_opts)

        ini.params['file_root'] = jobItem.chainRoot

        covmat = batch.basePath + 'planck_covmats/' + jobItem.name + '.covmat'
        if not os.path.exists(covmat):
            if hasattr(settings, 'covmat'): covmat = batch.basePath + settings.covmat
        if os.path.exists(covmat):
            ini.params['propose_matrix'] = covmat
            if settings.newCovmats: ini.params['MPI_Max_R_ProposeUpdate'] = 20
        else:
            hasCov = False
            ini.params['MPI_Max_R_ProposeUpdate'] = 20
            covmat_try = []
            if 'covRenamer' in dir(settings): covmat_try += settings.covRenamer(jobItem.name)
            if hasattr(settings, 'covrenames'):
                covmat_try += [jobItem.name.replace(old, new, 1) for old, new in settings.covrenames if old in jobItem.name]
                for new1, old1 in settings.covrenames:
                    if old1 in jobItem.name:
                        name = jobItem.name.replace(old1, new1, 1)
                        covmat_try += [name.replace(old, new, 1) for old, new in settings.covrenames if old in name]

            for name in covmat_try:
                covmat = batch.basePath + 'planck_covmats/' + name + '.covmat'
                if os.path.exists(covmat):
                    ini.params['propose_matrix'] = covmat
                    print 'covmat ' + jobItem.name + ' -> ' + name
                    hasCov = True
                    break
            if not hasCov: print 'WARNING: no matching specific covmat for ' + jobItem.name

        ini.params['start_at_bestfit'] = settings.start_at_bestfit
        updateIniParams(ini, jobItem.data_set.params, batch.commonPath)
        for deffile in settings.defaults:
            ini.defaults.append(batch.commonPath + deffile)

        ini.params['action'] = cosmomcAction
        ini.saveFile(jobItem.iniFile())
        if not settings.start_at_bestfit:
            setMinimize(jobItem, ini)
            variant = '_minimize'
            ini.saveFile(jobItem.iniFile(variant))


# add ini files for importance sampling runs
        for imp in jobItem.importanceJobs():
            if batch.hasName(imp.name.replace('_post', '')): raise Exception('importance sampling something you already have?')
            for minimize in (False, True):
                ini = iniFile.iniFile()
                updateIniParams(ini, imp.importanceSettings, batch.commonPath)
                if cosmomcAction == 0 and not minimize:
                    for deffile in settings.importanceDefaults:
                        ini.defaults.append(batch.commonPath + deffile)
                    ini.params['redo_outroot'] = imp.chainRoot
                    ini.params['action'] = 1
                else:
                    ini.params['file_root'] = imp.chainRoot
                if minimize:
                    setMinimize(jobItem, ini)
                    variant = '_minimize'
                else: variant = ''
                ini.defaults.append(jobItem.iniFile())
                ini.saveFile(imp.iniFile(variant))
                if cosmomcAction != 0: break


comment = 'Done... to run do: python python/runbatch.py ' + batchPath
print comment
if not settings.start_at_bestfit:
    comment = '....... for best fits: python python/runbatch.py ' + batchPath + ' --minimize'
    print comment
print ''
comment = 'for importance sampled: python python/runbatch.py ' + batchPath + ' --importance'
print comment
comment = 'for best-fit for importance sampled: python python/runbatch.py ' + batchPath + ' --importance_minimize'
print comment
