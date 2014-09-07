import subprocess, os, numpy as np

def addArguments(parser):
    parser.add_argument('--nodes', type=int, default=1)
    parser.add_argument('--chainsPerNode', type=int, default=4)
    parser.add_argument('--coresPerNode', type=int, default=16)
    parser.add_argument('--walltime', default='24:00:00')
    parser.add_argument('--job_template', default='job_script')
    parser.add_argument('--dryrun', action='store_true')

def getArgsOmp(args, msg=True):
    omp = args.coresPerNode / args.chainsPerNode
    if omp != np.floor(omp): raise Exception('Chains must each have equal number of cores')
    if msg:
        print 'Job parameters: %i nodes, each with %i MPI chains, each chain using %i OpenMP cores (%i cores per node)' % (args.nodes,
            args.chainsPerNode, omp, args.coresPerNode)
    return omp


def queued_jobs():
    res = subprocess.check_output('qstat -u $USER', shell=True)
    res = res.split("\n")
    names = []
    for line in res[4:]:
        if 'Q 00:00' in line:
            items = line.split()
            jobid = items[0].split('.')[0]
            output = subprocess.check_output('qstat -f ' + str(jobid) , shell=True).split('\n')
            pars = []
            current = ''
            for L in output:
                if '=' in L:
                    if len(current) > 0:
                        pars.append(current)
                    current = L.strip()
                else: current += L.strip()
            if len(current) > 0: pars.append(current)
            props = dict()
            for L in pars[1:]:
                (key, val) = L.split('=', 1)
                props[key.strip()] = val.strip()
            names.append(props['Job_Name'])
    return names


def queued_jobs_PBS():  # what used to use on Darwin
    res = subprocess.check_output('qstat -u $USER', shell=True)
    res = res.split("\n")
    names = []
    for line in res[4:]:
        if 'master' in line:
            items = line.split()
            jobid = items[0].split('.')[0]
            output = subprocess.check_output('qstat -f ' + str(jobid) , shell=True).split('\n')
            pars = []
            current = ''
            for L in output:
                if '=' in L:
                    if len(current) > 0:
                        pars.append(current)
                    current = L.strip()
                else: current += L.strip()
            if len(current) > 0: pars.append(current)
            props = dict()
            for L in pars[1:]:
                (key, val) = L.split('=', 1)
                props[key.strip()] = val.strip()
            names.append(props['Job_Name'])
    return names

def replacePlaceholders(txt, vals):
    txt = txt.replace('\r', '')
    for name, value in vals.iteritems():
        txt = txt.replace('##' + name + '##', str(value))
    return txt

def submitJob(jobName, paramFiles, pbs_template='job_script', numnodes=1, omp=4, chainsPerNode=1, mem_per_node=63900, walltime='24:00:00',
               qsub='qsub'):
    ppn = chainsPerNode * omp
    nchains = numnodes * chainsPerNode
    mem = mem_per_node * numnodes
    vals = dict()
    vals['JOBNAME'] = jobName
    vals['OMP'] = omp
    vals['MEM_MB'] = mem
    vals['WALLTIME'] = walltime
    vals['NUMNODES'] = numnodes
    vals['PPN'] = ppn
    vals['NUMMPI'] = nchains
    vals['ROOTDIR'] = os.getcwd()
    commands = []
    if isinstance(paramFiles, basestring): paramFiles = [paramFiles]
    for param in paramFiles:
        ini = param
        if ini[-4:] != '.ini': ini += '.ini'
        commands.append('time mpirun -ppn %i -np %i ./cosmomc %s > ./scripts/%s.log 2>&1' % (chainsPerNode, nchains, ini, jobName))
    vals['COMMAND'] = "\n".join(commands)
    script = replacePlaceholders(open(pbs_template, 'r').read(), vals)
    scriptName = './scripts/' + jobName + '_subscript'
    open(scriptName, 'w').write(script)
    if len(paramFiles) > 1:
        open('./scripts/' + jobName + '.batch', 'w').write("\n".join(paramFiles))
    os.system(qsub + ' ' + scriptName)

