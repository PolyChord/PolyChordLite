#!/usr/bin/env python
import os, jobQueue, batchJobArgs

parser = batchJobArgs.argParser("Submit a single MPI job")

parser.add_argument('iniFile')

jobQueue.addArguments(parser)

args = parser.parse_args()

omp = jobQueue.getArgsOmp(args, msg=True)
ini = args.iniFile.replace('.ini', '')

if not args.dryrun:
    jobQueue.submitJob(os.path.basename(ini), ini, numnodes=args.nodes, chainsPerNode=args.chainsPerNode,
                               omp=omp, walltime=args.walltime, pbs_template=args.job_template)
