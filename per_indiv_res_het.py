import sys

samples = ()
for line in sys.stdin:
    if line.startswith('##'):
        continue
    elif line.startswith('#'):
        samples = tuple([[x, 0] for x in line.split()[9:]])
    else:
        gts = [x.split(':')[0].split('/') for x in line.split()[9:]]
        for i in range(0, len(gts)):
            if len(set(gts[i])) == 2:
                samples[i][1] += 1

for x in samples:
    print '\t'.join([str(y) for y in x])
