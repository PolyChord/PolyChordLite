import os, fnmatch, batchJobArgs

Opts = batchJobArgs.batchArgs('delete failed chains, files etc.', importance=True, converge=True)

Opts.parser.add_argument('--dist', action='store_true')
Opts.parser.add_argument('--ext', nargs='+', default=['*'])
Opts.parser.add_argument('--empty', action='store_true')
Opts.parser.add_argument('--confirm', action='store_true')
Opts.parser.add_argument('--chainnum', default=None)

(batch, args) = Opts.parseForBatch()

def fsizestr(fname):
    sz = os.path.getsize(fname) / 1024
    if (sz < 1024): return str(sz) + 'KB'
    if (sz < 1024 * 1024): return str(sz / 1024) + 'MB'
    if (sz < 1024 * 1024 * 1024): return str(sz / 1024 / 1024) + 'GB'

if args.chainnum is not None:
    args.ext = ['_' + args.chainnum + '.' + ext for ext in args.ext]
else: args.ext = ['.' + ext for ext in args.ext] + ['_*.' + ext for ext in args.ext]

for jobItem in Opts.filteredBatchItems():
    if (args.converge == 0 or not jobItem.hasConvergeBetterThan(args.converge, returnNotExist=True)) and os.path.exists(jobItem.chainPath):
        dirs = [jobItem.chainPath]
        if args.dist: dirs = []
        if os.path.exists(jobItem.distPath): dirs += [jobItem.distPath]
        for adir in dirs:
            files = sorted(os.listdir(adir))
            for f in files:
                for ext in args.ext:
                    if fnmatch.fnmatch(f, jobItem.name + ext):
                        fname = adir + f
                        if os.path.exists(fname):
                            if not args.empty or os.path.getsize(fname) == 0:
                                print fname, ' (' + fsizestr(fname) + ')'
                                if args.confirm: os.remove(fname)

if not args.confirm: print 'Files not actually deleted: add --confirm to delete'
