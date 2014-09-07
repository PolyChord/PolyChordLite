import os, sys, batchJob


if len(sys.argv) < 3:
    print 'Usage: python/addGridBatch.py directory_with_outputs directory_with_output_to_add [and_another..]'
    sys.exit()

batch = batchJob.readobject()

for subBatch in sys.argv[2:]:
    batchPath2 = os.path.abspath(subBatch) + os.sep
    batch2 = batchJob.readobject(batchPath2)
    batch.subBatches.append(batch2)

for jobItem in batch.items():
    for x in [imp for imp in jobItem.importanceJobs()]:
        if batch.hasName(x.name.replace('_post', '')):
            print 'replacing importance sampling run (not deleting files): ' + x.name
            jobItem.importanceItems.remove(x)

batch.save()

