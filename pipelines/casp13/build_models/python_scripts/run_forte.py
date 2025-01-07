import subprocess
import glob
import sys
import multiprocessing
import itertools
import signal
import numpy as np
from retry import retry
import os
import errno

# fasfiles=sorted(glob.glob("./data/*/*.fas"))
# db_ascii_files=sorted(glob.glob("./data/*/*.fas.db.ascii"))
# pspmdd_files=sorted(glob.glob("./data/*/*.fas.pspmdd"))

methods = [
    ("~/forte4DB/runforteCC","{query_seq} {query_ascii} {db_seq_list} {db_ascii_list} {tmpfile}","f"),
]
sigmethods = [
    #("~/forte4DB/runsigforteCC","{query_seq} {query_ascii} {db_seq_list} {db_ascii_list} 0.15 -17 -0.5 3 1.2 -1 {tmpfile}", "s"),
    ("~/forte4DB/runsigforteCC","{query_seq} {query_ascii} {db_seq_list} {db_ascii_list} 0.1 -32.0 -0.0 3.0 1.0 -1.0 {tmpfile}", "s"),
]


db_path = "/home/tsukasa/casp_database/forte/list/"
dbs = [
   ("ssp_ascii.seq_list","ssp_ascii.ascii_list", "ssp"),
    ("pdb098_db.seq_list","pdb098_db.ascii_list","db"),
#    ("pdb100_dalipsii.seq_list","pdb100_dalipsii.ascii_list","dap"),
   ("pdp_ascii.seq_list","pdp_ascii.ascii_list", "pdpp"),

   ("ssp_pspm.seq_list","ssp_pspm.ascii_list", "sspx"),
    # ("pdb100_sspsi_pspm.seq_list00","pdb100_sspsi_pspm.ascii_list00", "ssearch_psiblast_pspm00"),
    # ("pdb100_sspsi_pspm.seq_list01","pdb100_sspsi_pspm.ascii_list01", "ssearch_psiblast_pspm01"),
    ("pdb098_pspmdb.seq_list","pdb098_pspmdb.ascii_list", "dbx"),
    # ("pdb100_pspmdd.seq_list00","pdb100_pspmdd.ascii_list00", "deltablast_pspm00"),
    # ("pdb100_pspmdd.seq_list01","pdb100_pspmdd.ascii_list01", "deltablast_pspm01"),
#    ("pdb100_dalipsi_pspm.seq_list","pdb100_dalipsi_pspm.ascii_list", "dappspm"),
   ("pdp_pspm.seq_list","pdp_pspm.ascii_list", "pdppx"),
    ("pdb098_pspmhh.seq_list","pdb098_pspmhh.ascii_list", "hhpx"),
    # ("pdb100_pspmhh.seq_list00","pdb100_pspmhh.ascii_list00", "hhblists_psiblast_pspm00"),
    # ("pdb100_pspmhh.seq_list01","pdb100_pspmhh.ascii_list01", "hhblists_psiblast_pspm01"),
#    ("scop98_pspmhh.seq_list","scop98_pspmhh.ascii_list", "schhpspm"),
]

sigdbs = [
    ("ssp_ascii.seq_list","ssp_ascii.ascii_list", "ssp"),
    # ("pdb100_db.seq_list","pdb100_db.ascii_list","deltablast"),
    # ("pdb100_dalipsii.seq_list","pdb100_dalipsii.ascii_list","dali_psiblast"),
    ("pdp_ascii.seq_list","pdp_ascii.ascii_list", "pdpp"),

    ("ssp_pspm.seq_list","ssp_pspm.ascii_list", "sspx"),
#    ("pdb100_sspsi_pspm.seq_list00","pdb100_sspsi_pspm.ascii_list00", "ssppspm00"),
#    ("pdb100_sspsi_pspm.seq_list01","pdb100_sspsi_pspm.ascii_list01", "ssppspm01"),
     ("pdb098_pspmdb.seq_list","pdb098_pspmdb.ascii_list", "dbx"),
#    ("pdb100_pspmdd.seq_list00","pdb100_pspmdd.ascii_list00", "dbpspm00"),
#    ("pdb100_pspmdd.seq_list01","pdb100_pspmdd.ascii_list01", "dbpspm01"),
#    ("pdb100_dalipsi_pspm.seq_list","pdb100_dalipsi_pspm.ascii_list", "dappspm"),
   ("pdp_pspm.seq_list","pdp_pspm.ascii_list", "pdppx"),
     ("pdb098_pspmhh.seq_list","pdb098_pspmhh.ascii_list", "hhpx"),
#    ("pdb100_pspmhh.seq_list00","pdb100_pspmhh.ascii_list00", "hhppspm00"),
#    ("pdb100_pspmhh.seq_list01","pdb100_pspmhh.ascii_list01", "hhppspm01"),
#    ("scop98_pspmhh.seq_list","scop98_pspmhh.ascii_list", "schhpspm"),
]

query_path = "/home/tsukasa/casp13/"
querys = [
    (".fas", ".fas.asciiss", "ssp"),
    (".fas", ".fas.asciiss3", "ssp3"),
    (".fas", ".fas.pspmsspsi", "sspx"),
    (".fas", ".fas.pspmsspsi3", "sspx3"),

    (".fas", ".fas.asciihhpsi", "hhp"),
    (".fas", ".fas.asciidb", "db"),
    (".fas", ".fas.asciipb", "pb"),

    (".fas", ".fas.pspmhh", "hhx"),
    (".fas", ".fas.pspmhhpsi", "hhpx"),
    (".fas", ".fas.pspmdb", "dbx"),
    (".fas", ".fas.pspmpb", "px"),
]

sigquerys = [
    (".fas", ".fas.asciiss", "ssp"),
    (".fas", ".fas.asciiss3", "ssp3"),
    (".fas", ".fas.pspmsspsi", "sspx"),
    (".fas", ".fas.pspmsspsi3", "sspx3"),

    (".fas", ".fas.asciihhpsi", "hhp"),
    (".fas", ".fas.asciidb", "db"),
    (".fas", ".fas.asciipb", "pb"),

    (".fas", ".fas.pspmhh", "hhx"),
    (".fas", ".fas.pspmhhpsi", "hhpx"),
    (".fas", ".fas.pspmdb", "dbx"),
    (".fas", ".fas.pspmpb", "px"),
]

# targets = ("C0018", )
targets = ( os.path.basename(os.path.dirname(os.path.abspath(sys.argv[0]))), )

# jobscheduler_header = "qrsh -b y -now n -l h='strinfo01' "
# jobscheduler_header = "qrsh -b y -now n -l h='!(strinfo01|strinfo09)' "
jobscheduler_header = "qrsh -b y -now n -cwd -V -l h='!(h002|gp001|gp002|griffin-l)' "


def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


@retry(exceptions=subprocess.CalledProcessError,tries=3,delay=2,backoff=2)
def exec_wrap(args):
    return exec_forte(args)

def exec_forte(args):
    target, method, db, query  = args
    outfilename = f"{target}-{method[2]}_{db[2]}_{query[2]}.res.head1000"
    outfilepath = "/home/tsukasa/casp13/"+target+"/head1000/"+outfilename
    mkdir_p("/home/tsukasa/casp13/"+target+"/head1000/")
    if os.path.isfile(outfilepath) and len(open(outfilepath).read().splitlines())==1000:
        return
    args_dict = {
#        "query_seq": query_path+target+"/query_profile/"+target.split("_")[0]+query[0],
        "query_seq": query_path+target+"/query_profile/"+target+query[0],
#        "query_ascii": query_path+target+"/query_profile/"+target.split("_")[0]+query[1],
        "query_ascii": query_path+target+"/query_profile/"+target+query[1],
        "db_seq_list": db_path+db[0],
        "db_ascii_list": db_path+db[1],
        "tmpfile": "$TMPFILE",
    }
    command_line = method[0] + " " + method[1].format(**args_dict)
    command_line = jobscheduler_header + " " + "'TMPFILE=`mktemp`; " + command_line + " ' < /dev/null "
    print(command_line, flush=True)
    print(outfilepath, flush=True)
    forteoutput = subprocess.check_output(command_line, shell=True)
    forteoutput = forteoutput.decode('utf-8')
    # print(forteoutput.splitlines()[:20])
    if not "Query=" in forteoutput.splitlines()[0]:
        print(outfilepath)
        print(forteoutput[:20])
    with open(outfilepath,"w") as f:
        f.write("\n".join(forteoutput.splitlines()[:1000]))
        f.write("\n")

def init_worker():
    signal.signal(signal.SIGINT, signal.SIG_IGN)


pool = multiprocessing.Pool(processes=100)
try:
    _ = pool.map_async(exec_wrap, itertools.product(targets,methods,dbs,querys), chunksize=1).get(9999999)
    _ = pool.map_async(exec_wrap, itertools.product(targets,sigmethods,sigdbs,sigquerys), chunksize=1).get(9999999)
    # http://stackoverflow.com/questions/1408356/keyboard-interrupts-with-pythons-multiprocessing-pool
    # http://stackoverflow.com/questions/11312525/catch-ctrlc-sigint-and-exit-multiprocesses-gracefully-in-python
    pool.close()
    pool.join()
except KeyboardInterrupt:
    print("Caught KeyboardInterrupt, terminating workers")
    pool.terminate()
    pool.join()

#list(map(exec_forte, itertools.product(targets,methods,dbs,querys)))
#list(map(exec_forte, itertools.product(targets,sigmethods,sigdbs,sigquerys)))



