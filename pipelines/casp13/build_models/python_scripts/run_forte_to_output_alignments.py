import subprocess
import glob
import sys
import multiprocessing
import itertools
import signal
import numpy as np
from retry import retry
import os
import tempfile
import errno

methods = [
    ("~/forte4DB/bin/df31715nc_iccf -f", "{query_seq} {query_ascii} {db_seq_list} {db_ascii_list} ","f"),
]
sigmethods = [
    #("~/forte4DB/bin/kf7lscoremod0fR2i -f", "{query_seq} {query_ascii} {db_seq_list} {db_ascii_list} 0.15 -17 -0.5 3 1.2 -1 ", "s"),
    #("~/forte4DB/bin/locsigpal -f", "{query_seq} {query_ascii} {db_seq_list} {db_ascii_list} 0.1 -32.0 -0.0 3.0 1.0 -1.0 ", "s"),
    ("~/forte4DB/bin/locsigpal5064 -f", "{query_seq} {query_ascii} {db_seq_list} {db_ascii_list} 0.1 -32.0 -0.0 3.0 1.0 -1.0 ", "s"),
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

path_replace_dict = {
    "ssp":("pdb098.ss.seq","pdb098.ss.ascii",".fas",".fas.asciipb"),
    "db":("pdb098.seq","pdb098.asciidb",".fas",".fas.asciidb"),
    "dap":("pdb100.dali.seq","pdb100.dali.ascii","sequence.fas","sequence.dali.ascii"),
    "pdpp":("pdp_ali.aseq","pdp_ali.ascii",".fas",".fas.asciipb"),
    "sspx":("pdb098.ss.pseq","pdb098.ss.pspm",".fas",".fas.pspmpb"),
    "dbx":("pdb098.seq","pdb098.pspmdb",".fas",".fas.pspmdb"),
    "dapx":("pdb100.dali.seq","pdb100.dali.pspm","sequence.fas","sequence.pspmdali_psi"),
    "pdppx":("pdp_ali.seq","pdp_ali.pspm",".fas",".fas.pspmpb"),
    "hhpx":("pdb098.seq","pdb098.pspmhh",".fas",".fas.pspmhh_psi"),
    "schhx":("scop98.seq","scop98.pspmhh",".fas",".pspm"),

    "ssppspm00":("pdb100.ss.pseq","pdb100.ss.pspm","sequence.fas","sequence.pspmss_psi"),
    "ssppspm01":("pdb100.ss.pseq","pdb100.ss.pspm","sequence.fas","sequence.pspmss_psi"),
    "dbpspm00":("pdb100.seq","pdb100.pspmdd","sequence.fas","sequence.pspm"),
    "dbpspm01":("pdb100.seq","pdb100.pspmdd","sequence.fas","sequence.pspm"),
    "hhppspm00":("pdb100.seq","pdb100.pspmhh","sequence.fas","sequence.pspm"),
    "hhppspm01":("pdb100.seq","pdb100.pspmhh","sequence.fas","sequence.pspm"),
}

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

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise

#@retry(exceptions=subprocess.CalledProcessError,tries=3,delay=2,backoff=2)
#def exec_wrap(args):
#    return exec_forte(args)

def exec_forte(args):
    target, method, db, query  = args
    headfilename = f"{target}-{method[2]}_{db[2]}_{query[2]}.res.head1000"
    headfilepath = "/home/tsukasa/casp13/"+target+"/head1000/"+headfilename
    if not os.path.isfile(headfilepath) or len(open(headfilepath, "r").read().splitlines())!=1000:
        print("Error :", headfilepath)
        return
    with open(headfilepath, "r") as f:
        headresult = [l.split()[1] for l in f.read().strip().splitlines()[1:51]]

    seqtf=tempfile.NamedTemporaryFile("w");
    seqtf.write("\n".join(headresult)+"\n")
    asciitf=tempfile.NamedTemporaryFile("w");
    prd = path_replace_dict[db[2]]
    asciitf.write("\n".join([l.replace(prd[0],prd[1]).replace(prd[2],prd[3]) for l in headresult])+"\n")
    seqtf.flush()
    asciitf.flush()

    args_dict = {
#        "query_seq": query_path+target+"/query_profile/"+target.split("_")[0]+query[0],
        "query_seq": query_path+target+"/query_profile/"+target+query[0],
#        "query_ascii": query_path+target+"/query_profile/"+target.split("_")[0]+query[1],
        "query_ascii": query_path+target+"/query_profile/"+target+query[1],
        "db_seq_list": seqtf.name,
        "db_ascii_list": asciitf.name,
    }
    command_line = method[0] + " " + method[1].format(**args_dict)
    print(command_line)
    sys.stdout.flush()
    forteoutput = subprocess.check_output(command_line, shell=True)
    forteoutput = forteoutput.decode("utf-8")
    mkdir_p("/home/tsukasa/casp13/"+target+"/alignment")
    # alignmentfilename = "_".join([target,method[2],db[2],query[2]]) + ".ali"
    alignmentfilename = f"{target}-{method[2]}_{db[2]}_{query[2]}.ali"
    headfilepath = "/home/tsukasa/casp13/"+target+"/alignment/"+alignmentfilename
    with open(headfilepath,"w") as f:
        f.write(forteoutput+"\n")

    seqtf.close()
    asciitf.close()

#def init_worker():
#    signal.signal(signal.SIGINT, signal.SIG_IGN)


#pool = multiprocessing.Pool(processes=10)
list(map(exec_forte, itertools.product(targets,methods,dbs,querys)))
list(map(exec_forte, itertools.product(targets,sigmethods,sigdbs,sigquerys)))
# try:
#     _ = pool.map_async(exec_wrap, itertools.product(targets[:1],methods,dbs,querys)).get(9999999)
#     # http://stackoverflow.com/questions/1408356/keyboard-interrupts-with-pythons-multiprocessing-pool
#     # http://stackoverflow.com/questions/11312525/catch-ctrlc-sigint-and-exit-multiprocesses-gracefully-in-python
#     pool.close()
#     pool.join()
# except KeyboardInterrupt:
#     print("Caught KeyboardInterrupt, terminating workers")
#     pool.terminate()
#     pool.join()

