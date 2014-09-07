import batchJobArgs

Opts = batchJobArgs.batchArgs('Find chains which have failed or not converged.', importance=True, converge=True)

Opts.parser.add_argument('--exist', action='store_true')
Opts.parser.add_argument('--checkpoint', action='store_true')
Opts.parser.add_argument('--running', action='store_true')
Opts.parser.add_argument('--not_running', action='store_true')
Opts.parser.add_argument('--stuck', action='store_true')

(batch, args) = Opts.parseForBatch()

notExist = []
converge = []

if args.running:args.checkpoint = True

if args.stuck:
        for jobItem in Opts.filteredBatchItems():
            if jobItem.chainExists() and jobItem.chainsDodgy():
                print 'Chain stuck?...' + jobItem.name
elif args.checkpoint:
    print 'Convergence from checkpoint files...'
    for jobItem in Opts.filteredBatchItems():
        R, done = jobItem.convergeStat()
        if R is not None and not done:
            if (not args.not_running or jobItem.notRunning()) and (not args.running or not jobItem.notRunning()): print '...', jobItem.chainRoot, R
            if args.running and jobItem.chainExists() and jobItem.chainsDodgy():
                print 'Chain stuck?...' + jobItem.name
else:
    for jobItem in Opts.filteredBatchItems():
        if not jobItem.chainExists():
            notExist.append(jobItem)
        elif args.converge == 0 or args.checkpoint or not jobItem.hasConvergeBetterThan(args.converge, returnNotExist=True):
            if (not args.not_running or jobItem.notRunning()): converge.append(jobItem)

    print 'Checking batch (from last runGridGetdist.py output):'
    if not args.exist and len(notExist) > 0:
        print 'Not exist...'
        for jobItem in notExist:
            print '...', jobItem.chainRoot

    print 'Converge check...'
    for jobItem in converge:
        print '...', jobItem.chainRoot, jobItem.R()


