import glob,os,multiprocessing,subprocess
import scoop
from scoop import futures

os.environ["LD_LIBRARY_PATH"]="/home/tsukasa/lcl/lib64:/home/tsukasa/lcl/lib:/home/tsukasa/lcl/lib64:/home/tsukasa/lcl/lib:/home/intel//impi/5.0.3.048/intel64/lib:/home/intel/composer_xe_2015.2.164/compiler/lib/intel64:/home/intel/composer_xe_2015.2.164/mpirt/lib/intel64:/home/intel/composer_xe_2015.2.164/ipp/../compiler/lib/intel64:/home/intel/composer_xe_2015.2.164/ipp/lib/intel64:/home/intel/composer_xe_2015.2.164/ipp/tools/intel64/perfsys:/home/intel/composer_xe_2015.2.164/compiler/lib/intel64:/home/intel/composer_xe_2015.2.164/mkl/lib/intel64:/home/intel/composer_xe_2015.2.164/tbb/lib/intel64/gcc4.4:/home/intel/composer_xe_2015.2.164/debugger/ipt/intel64/lib:/home/tsukasa/bin/modeller9.19/lib/x86_64-intel8"


#pdbpaths=sorted(glob.glob("./worktmp2/*.pdb"))
#done_names = [ os.path.basename(p).split(".")[0] for p in pdbpaths]
#done_names_set = set([d for d in done_names if done_names.count(d)>0])

def f(path):
#    if os.path.splitext(os.path.basename(path))[0] in done_names_set:
#        return
    # if os.path.isfile("./worktmp/"+os.path.splitext(os.path.basename(path))[0]+".B99990001.pdb"):
    #     return
    p = subprocess.run(f"~/.pyenv/shims/python modeller_run.py {path} ./worktmp/", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    print("done", path, p.returncode)
    return p.returncode

if __name__ == '__main__':
    paths=sorted(glob.glob("./pir/*.pir"))
    list(futures.map(f,paths))
