import json

with open('template.json', 'r') as file:
    tmpl = json.load(file)

jobs = []
for l in open("target.list", 'r').readlines():
    es = l.strip().split()
    name = es[0]
    seq_pro = es[1]
    seq_pep = es[2]

    job = json.loads(json.dumps(tmpl[0]))

    job['name'] = name
    job['sequences'][0]['proteinChain']['sequence'] = seq_pro
    job['sequences'][1]['proteinChain']['sequence'] = seq_pep

    jobs.append(job.copy())

with open('pepdock_scaffolds.json', 'w') as f:
    f.write(json.dumps(jobs, indent=4))
