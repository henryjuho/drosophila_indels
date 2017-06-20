import sys

for line in sys.stdin:
    split_line = line.split('\t')
    try:
        start, end = int(split_line[3]), int(split_line[4])
        if start < 0 or end < 0:
            continue
        else:
            print line.rstrip('\n')
    except IndexError:
        print line.rstrip('\n')
