import sys
import anavar_utils as an


res_files = [x for x in sys.stdin]

an.merge_results(res_files, sys.argv[1])
