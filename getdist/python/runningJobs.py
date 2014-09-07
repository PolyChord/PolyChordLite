import subprocess


res = subprocess.check_output('qstat -u $USER', shell=True)

res = res.split("\n")
for line in res[4:]:
    if ' R ' in line:
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
        print jobid, " ".join(items[5:]), props['Job_Name']



