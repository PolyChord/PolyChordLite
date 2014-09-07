import os, fnmatch, shutil, batchJobArgs

Opts = batchJobArgs.batchArgs('copy all files of a given type from all getdist output directories in the batch', importance=True, converge=True)

Opts.parser.add_argument('target_dir')
Opts.parser.add_argument('file_extension', nargs='+')


(batch, args) = Opts.parseForBatch()

target_dir = os.path.abspath(args.target_dir) + os.sep
if not os.path.exists(target_dir): os.makedirs(target_dir)

for ext in args.file_extension:
    if not '.' in ext: pattern = '.' + ext
    else: pattern = ext
    for jobItem in Opts.filteredBatchItems():
        if os.path.exists(jobItem.distPath) and (args.converge == 0 or jobItem.hasConvergeBetterThan(args.converge)):
            for f in os.listdir(jobItem.distPath):
                if fnmatch.fnmatch(f, jobItem.name + pattern):
                    print jobItem.distPath + f
                    shutil.copyfile(jobItem.distPath + f, target_dir + f)
