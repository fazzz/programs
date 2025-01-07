import subprocess
import os
import re
import sys



if len(sys.argv) != 3:
    print("Usage: {} infile outfile".format(sys.argv[0]))
    sys.exit(1)

# os.chdir("/home/tsukasa/lcl/scripts/")
fp = os.path.abspath(sys.argv[1])
outpath = os.path.abspath(sys.argv[2])
if os.path.isfile(outpath) and os.path.getsize(outpath) != 0:
    print("outfile exists")
    sys.exit(1)

#subprocess.call("perl -I/home/tsukasa/lcl/scripts/ /home/tsukasa/lcl/scripts/reformat.pl a3m fas {} {} ".format(fp, outpath), shell=True)
subprocess.call("perl -I/home/yamamori/opt2/lib/hh/scripts /home/yamamori/opt2/lib/hh/scripts/reformat.pl a3m fas {} {} ".format(fp, outpath), shell=True)

# # del gap column based on query sequence
# seqs = [">"+l for l in re.split(r"^>", open(outpath).read(), flags=re.M) if l is not ""]

# not_gap_index = [i for i, s in enumerate("".join(seqs[0].splitlines()[1:])) if s != "-"]
# out_seqs = []
# for j, l in enumerate(seqs):
#     ll = l.splitlines()
#     ss = "".join(ll[1:])
#     sss = "".join([ss[i] for i in not_gap_index])
#     if all((s == "-" for s in sss)) is True:
#         print("all gaps ", id, j)
#         # break
#         continue
#     s = "".join([ll[0], "\n", sss])
#     # print(s)
#     out_seqs.append(s)

# with open(outpath, "w") as f:
#     f.write("\n".join(out_seqs))
