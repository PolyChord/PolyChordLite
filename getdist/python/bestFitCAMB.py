import ResultObjs, sys, os, iniFile


if len(sys.argv) < 3:
    print 'Usage: python/bestFitCAMB.py chain_root iniName'
    sys.exit()

root = os.path.abspath(sys.argv[1])

pars = {'ombh2':'omegabh2', 'omch2':'omegach2', 'omnuh2':'omeganuh2', 'hubble':'H0', 'w':'w',
        'helium_fraction':'yheused', 'scalar_amp(1)':'A' , 'scalar_spectral_index(1)':'ns', 'scalar_nrun(1)':'nrun', 'initial_ratio(1)':'r',
        're_optical_depth':'tau', 're_delta_redshift':'deltazrei', 'massless_neutrinos':'nnu'}

ini = iniFile.iniFile()

ini.params['re_use_optical_depth'] = True
ini.params['temp_cmb'] = 2.7255
ini.params['CMB_outputscale'] = 2.7255e6 ** 2.
ini.defaults.append('params.ini')

bf = ResultObjs.bestFit(root + '.minimum', setParamNameFile=root + '.paramnames', want_fixed=True)

for camb, cosmomc in pars.items():
    par = bf.parWithName(cosmomc)
    if par is not None: ini.params[camb] = par.best_fit

ini.params['scalar_amp(1)'] = float(ini.params['scalar_amp(1)']) / 1e9

nmassive = 1
neffstandard = 3.046 / 3
ini.params['massless_neutrinos'] = float(ini.params['massless_neutrinos']) - neffstandard * nmassive
ini.params['massive_neutrinos'] = int(round(neffstandard * nmassive))
ini.params['nu_mass_degeneracies'] = neffstandard * nmassive
ini.params['tensor_spectral_index(1)'] = -float(ini.params['initial_ratio(1)']) / 8


ini.saveFile(sys.argv[2])

print 'OK, though note this does not support all parameter extensions from LCDM'
