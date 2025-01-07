import glob
import os
fps = []
dirnames = sorted(glob.glob("T*"))
dirnames = [ p for p in dirnames if "." not in os.path.basename(p)]
for d in ["T0970_2dtj", "T0976", "T0977_all"]:
    if d in dirnames:
        dirnames.remove(d)
for d in sorted(glob.glob("T*_casp12")):
    dirnames.remove(d)
#dirnames = ["T0994_N21","T0995_C22","T0996_N12","T0997_N47","T0998","T0999_N15C12"]
#dirnames = sorted(glob.glob("T10*"))
#dirnames = ["Ehasubunit"]
# dirnames = ["T0999_16F402"]
#dirnames = ["ScilNext"]
dirnames = ["NtpC","NtpI"]
#dirnames = sorted(glob.glob("c1N_IGFBP_*"))
#dirnames = sorted(glob.glob("TmBER?_rep"))
fps += [ f"{dirname}/query_profile/{dirname}.fas" for dirname in dirnames]
#fps = reversed(fps)
# fps += sorted(glob.glob("./T????/query_profile/T????.fas"))
# fps += sorted(glob.glob("./T????_*/query_profile/T????_*.fas"))
# fps=sorted(glob.glob("./T????all/query_profile/T????all.fas"))
# fps = [ "T0969_N131/query_profile/T0969_N131.fas", "T0970/query_profile/T0970.fas"]

import subprocess
import sys
import shutil

thread_num = 20

def parse_blastf_to_fas(l):
    if "Gap Penalties" not in open(l).read().splitlines()[-1]:
        print("Gap Penalties not in " + l)
        return
    # do blast to fasta
    subprocess.call(["perl", "-I/home/tsukasa/casp13/scripts/", "/home/tsukasa/casp13/scripts/ssearch_blast_to_fasta.pl",
                      l, os.path.dirname(l), "0.005"])
    fp = f"{l}.fas"
    if not os.path.isfile(fp):
        print(f"Error: {l} !! No sequences with E() < 0.005")
        return
    # check first seq is pdb_master_seq
    fasta = [">"+l.strip() for l in open(fp).read()[1:].split("\n>") if l is not ""]
    if len(fasta) < 2:
        print(f"Error: {fp} , len(fasta) < 2")
        print(fasta)
        os.remove(fp)
        return
    nr_blastf_parse_first = "".join(fasta[0].splitlines()[1:])
#    master_seq = "".join(open(os.path.dirname(l)+"/"+os.path.basename(l).split("_")[0]+".fas", "r").read().splitlines()[1:])
    master_seq = "".join(open(os.path.dirname(l)+"/"+os.path.normpath(os.path.dirname(l)).split(os.sep)[0]+".fas", "r").read().splitlines()[1:])
    if nr_blastf_parse_first != master_seq:
        print(l + " : nr_blastf_parse_first != master_seq" + "\n" + nr_blastf_parse_first+"\n"+master_seq)
        assert len(nr_blastf_parse_first) < len(master_seq)
        assert nr_blastf_parse_first in master_seq, f"\n{nr_blastf_parse_first}\n{master_seq}"
        # query no baai seq ga kotonaruto error nintsunagarunode murinidemo query ha awaseru
        fasta[0] = fasta[0].splitlines()[0] +"\n"+ master_seq
    # delete gaps
    fasta = [item.splitlines()[0] +"\n"+ item.splitlines()[1].replace("-", "") for item in fasta]
    with open(fp, "w") as outf:
        outf.write("\n".join(fasta))


for i, fp in enumerate(fps):
    dirpath = os.path.dirname(fp)
    id = os.path.splitext(os.path.basename(fp))[0]

    # hh
    if not all(map(lambda p:(os.path.isfile(p) and os.path.getsize(p) != 0), map(lambda ext:os.path.join(dirpath, id + ext), [".fas.a3m", ".fas.asciihhpsi", ".fas.hhm", ".fas.hhr", ".fas.pspmhh", ".fas.pspmhhpsi"]))):
#        command_line = "perl /home/tsukasa/casp13/scripts/make_pspm_from_hhblitsresult_all.pl {}".format(fp)
        command_line = "perl /home/tsukasa/casp13/scripts/make_pspm_from_hhblitsresult.pl {}".format(fp)
        print(command_line, flush=True)
        subprocess.call(command_line, shell=True)

    # ss prediction RaptorX_property
    if not os.path.isdir(os.path.join(dirpath, "ss")):
        command_line = "bash /home/tsukasa/RaptorX_Property_Fast/oneline_command_not_use_hhblits.sh {} {} {}".format(fp, os.path.abspath(fp+".a3m"), os.path.abspath(dirpath))
        print(command_line, flush=True)
        subprocess.call(command_line, shell=True)
        shutil.move(os.path.join(dirpath, id), os.path.join(dirpath, "ss"))

    # ssm
    if not all(map(lambda p:(os.path.isfile(p) and os.path.getsize(p) != 0), map(lambda ext:os.path.join(dirpath, id + ext), ["_nr_ssearch_blastf.output", "_nr_ssearch_parsef.output", "_nr_ssearch_ssearchf.output"]))):
        command_line = "/home/tsukasa/fasta-36.3.8g/bin/ssearch36 -q -m 'FBB {}' -m 'F10 {}' -E 0.05 -T {} -s '/home/tsukasa/miqs.mat' {} /home/tsukasa/casp_database/nr/nr_removed_soh > {} " \
                       .format(os.path.join(dirpath, f"{id}_nr_ssearch_blastf.output"), os.path.join(dirpath, f"{id}_nr_ssearch_parsef.output"), thread_num, fp, os.path.join(dirpath, f"{id}_nr_ssearch_ssearchf.output"))
        print(command_line, flush=True)
        subprocess.call(command_line, shell=True)
    if not all(map(lambda p:(os.path.isfile(p) and os.path.getsize(p) != 0), map(lambda ext:os.path.join(dirpath, id + ext), ["_nr_ssearch_blastf.output.fas"]))):
        parse_blastf_to_fas(os.path.join(dirpath, f"{id}_nr_ssearch_blastf.output"))
    if not all(map(lambda p:(os.path.isfile(p) and os.path.getsize(p) != 0), map(lambda ext:os.path.join(dirpath, id + ext), ["_nr_ssearch_blastf.msa"]))):
        my_env = os.environ.copy()
        my_env["MAFFT_N_THREADS_PER_PROCESS"] = str(thread_num)
        my_env["MAFFT_MPIRUN"] = "mpirun "
        command_line = "/home/tsukasa/mafft/v7386_20180216/MPI-MAFFT/tmp/bin/mafft --mpi --large --globalpair --aamatrix /home/tsukasa/miqs.mat --threadtb {} --anysymbol {} > {}".format(thread_num, os.path.join(dirpath, f"{id}_nr_ssearch_blastf.output.fas"), os.path.join(dirpath, f"{id}_nr_ssearch_blastf.msa"))
        print(command_line, flush=True)
        subprocess.call(command_line, shell=True, env=my_env)
    if not all(map(lambda p:(os.path.isfile(p) and os.path.getsize(p) != 0), map(lambda ext:os.path.join(dirpath, id + ext), [".fas.asciiss"]))):
        command_line = "/home/tsukasa/psiblastexb250 -num_iterations 1 -num_threads {} -num_descriptions 100000 -db /home/tsukasa/casp_database/blast_nr/nr -in_msa {} -out_ascii_pssm {} > /dev/null".format(thread_num, os.path.join(dirpath, f"{id}_nr_ssearch_blastf.msa"), os.path.join(dirpath, f"{id}.fas.asciiss"))
        print(command_line, flush=True)
        subprocess.call(command_line, shell=True)
    if not all(map(lambda p:(os.path.isfile(p) and os.path.getsize(p) != 0), map(lambda ext:os.path.join(dirpath, id + ext), [".fas.pspmsspsi"]))):
        command_line = "perl /home/tsukasa/casp13/scripts/makepspm.pl {} {} ".format(os.path.join(dirpath, f"{id}.fas.asciiss"), os.path.join(dirpath, f"{id}.fas.pspmsspsi"))
        print(command_line, flush=True)
        subprocess.call(command_line, shell=True)

    #ssm3
    if not all(map(lambda p:(os.path.isfile(p) and os.path.getsize(p) != 0), map(lambda ext:os.path.join(dirpath, id + ext), [".fas.asciiss3"]))):
        command_line = "/home/tsukasa/psiblastexb250 -num_iterations 3 -num_threads {} -num_descriptions 100000 -db /home/tsukasa/casp_database/blast_nr/nr -in_msa {} -out_ascii_pssm {} > /dev/null".format(thread_num, os.path.join(dirpath, f"{id}_nr_ssearch_blastf.msa"), os.path.join(dirpath, f"{id}.fas.asciiss3"))
        print(command_line, flush=True)
        subprocess.call(command_line, shell=True)
    if not all(map(lambda p:(os.path.isfile(p) and os.path.getsize(p) != 0), map(lambda ext:os.path.join(dirpath, id + ext), [".fas.pspmsspsi3"]))):
        command_line = "perl /home/tsukasa/casp13/scripts/makepspm.pl {} {} ".format(os.path.join(dirpath, f"{id}.fas.asciiss3"), os.path.join(dirpath, f"{id}.fas.pspmsspsi3"))
        print(command_line, flush=True)
        subprocess.call(command_line, shell=True)

    # db
    if not all(map(lambda p:(os.path.isfile(p) and os.path.getsize(p) != 0), map(lambda ext:os.path.join(dirpath, id + ext), [".fas.asciidb", ".fas.pspmdb"]))):
        command_line = "perl /home/tsukasa/casp13/scripts/make_pspm_from_deltablastresult.pl {}".format(fp)
        print(command_line, flush=True)
        subprocess.call(command_line, shell=True)

    #pb
    if not all(map(lambda p:(os.path.isfile(p) and os.path.getsize(p) != 0), map(lambda ext:os.path.join(dirpath, id + ext), [".fas.asciipb", ".fas.pspmpb"]))):
        command_line = "perl /home/tsukasa/casp13/scripts/make_pspm_from_psiblastresult.pl {}".format(fp)
        print(command_line, flush=True)
        subprocess.call(command_line, shell=True)
