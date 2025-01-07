#forte ali to pir
import re
import os
import Bio
import Bio.PDB
import Bio.pairwise2
import Bio.SeqIO
import multiprocessing
import glob

forte_ali_split_footer_line = re.compile(r"^(\d+[ ]+(?:/[^/ ]*)+[ ]+[-\d.]+)\n?",re.M)
def parse_forte_ali(f, path):
    forte_ali = [r for r in forte_ali_split_footer_line.split(f) if r != "\n" and r != ""]
    assert len(forte_ali) % 2 == 0

    ret = []
    for ali, footer in zip(forte_ali[::2], forte_ali[1::2]):
        query_seq = ""
        target_seq = ""
        if ali.startswith("query:"):
            #137 /work01/casp12/database/forte/pdb100.seq/4FA0_A_sequence.fas  649.907
            #query: 17-101
            #tmplt: 13-125
            #                +         +         +         +         +         +
            #201606 IVRAQRPTYVHWAI-----------RKVAPDGSAKQISLSRSGIQALVALEPPE----GE
            #/work0 LIEIFRPFYRHWAIYVGDGYVVHLAPPSEVAGAGAASVMSALTDKAIVKKELLYDVAGSD
            query_range = list(map(int,ali.splitlines()[0][6:].split("-")))
            target_range = list(map(int,ali.splitlines()[1][6:].split("-")))
            ali = "\n".join(ali.splitlines()[2:])
        else:
            query_range=None; target_range=None
        for item in re.split(r"\n\n", ali):
            if item == "": continue
            #'                +         +         +         +         +         +\n/home/ SVGDTKLPAITVDSSNNTLAGMRDAINQAGKEAGVSATIITDNSGSRLVLSSTKT-----EAGVSATIITDNSGSRLVLSSTKT-----'
             #                +         +         +         +         +         +
             #/home/ SVGDTKLPAITVDSSNNTLAGMRDAINQAGKEAGVSATIITDNSGSRLVLSSTKT-----
             #/work0 N-GVDVGKIDAASTAQERAAQLTEAINRVSSQTNVGASYDKTTGQVTLTSNAAIAVAGAA

            assert len(item.splitlines()) == 3, path
            plus_line, query_line, target_line = item.splitlines()
            assert all((c is " " or c is "+" for c in list(plus_line)))

            query_seq += query_line[7:]
            assert query_line.split(" ")[1] == query_line[7:]
            target_seq += target_line[7:]
            assert target_line.split(" ")[1] == target_line[7:]

        assert len(query_seq)==len(target_seq)
        assert len(re.split(r"\s+", footer))==3
        target_len, target_path, score = re.split(r"\s+", footer)
        if target_range:
            assert len(target_seq.replace("-", "")) == target_range[1]-target_range[0]+1
        else:
            assert len(target_seq.replace("-", "")) == int(target_len)

        ret.append((query_seq, target_seq, int(target_len), target_path, float(score)))
    return ret

def f(path):
    alis = parse_forte_ali(open(path, "r").read(), path)
    # print(path, len(alis))
    query_id, meth = os.path.splitext(os.path.basename(path))[0].split("-")
    for ali in alis[:20]:
        query_seq = ali[0]
        target_seq = ali[1]
        target_len = ali[2]
        target_path = ali[3]
        score = ali[4]
        if os.path.basename(target_path).startswith(("1","2","3","4","5","6","7","8","9")): # pdbid
            target_pdbid, target_chainid, = os.path.splitext(os.path.basename(target_path))[0].split("_")[:2]
        elif os.path.basename(target_path).startswith("g"):
            print(f"ERROR: cannot handle scopdomain startswith g {target_path}")
            continue
        elif os.path.basename(target_path).startswith("d"):
            # /work01/casp12/database/forte/scop98.seq/d4hv4b1.fas
            domainid = os.path.splitext(os.path.basename(target_path))[0]
            target_pdbid, target_chainid = domainid[1:5], domainid[5].upper()
            # scop chain id ga lowercase na reiha 8 rei sonzaisuruga koreraha toriaezu akirameru
        elif os.path.basename(target_path).startswith("PDP"):
            basename = os.path.splitext(os.path.basename(target_path))[0]
            assert basename[:4] == "PDP_", f"{path}\n{target_path}"
            domainid = basename[4:]
            if domainid.startswith("d"):
                if domainid in ["d2bez_1"]:
                    print(f"ERROR: cannot handle interchain scopdomain {domainid}")
                    return
                assert len(domainid) == 7, f"{path}\n{target_path}"
                # /work01/casp12/database/forte/pdp_ali.aseq/PDP_d2joza1.fas
                target_pdbid, target_chainid = domainid[1:5], domainid[5].upper()
                # scop chain id ga lowercase na reiha 8 rei sonzaisuruga koreraha toriaezu akirameru
            else:
                assert len(domainid) == 6, f"{path}\n{target_path}"
                # /work01/casp12/database/forte/pdp_ali.aseq/PDP_2VX8Aa.fas
                target_pdbid, target_chainid = domainid[0:4], domainid[4]
        else:
            assert False
        if target_chainid == "_":
            pdb_path = f"/home/tsukasa/work2/pdb/{target_pdbid[1:3].lower()}/pdb{target_pdbid.lower()}.ent"  \
                       if os.path.isfile(f"/home/tsukasa/work2/pdb/{target_pdbid[1:3].lower()}/pdb{target_pdbid.lower()}.ent") \
                       else f"/home/tsukasa/work2/pdb_obsolete/{target_pdbid[1:3].lower()}/pdb{target_pdbid.lower()}.ent" 
            target_chainid_l = [r.annotations["chain"] for r in Bio.SeqIO.parse(pdb_path, "pdb-seqres") if target_seq.replace("-", "") in r.seq]
            if len(target_chainid_l)!=1:
                target_chainid = next(Bio.SeqIO.parse(pdb_path, "pdb-seqres")).annotations["chain"]
                #rough guess
            else:
                target_chainid = target_chainid_l[0]
        zscore, rank = head1000dict[f"{query_id}-{meth}"][target_path]

        #check
        if not os.path.isfile(f"./pir/{query_id}-{meth}-{target_pdbid}_{target_chainid}.pir"):
            print(f"ERROR, rank, {rank} :  ./pir/{query_id}-{meth}-{target_pdbid}_{target_chainid}.pir not exists")
        else:
            if len(glob.glob(f"./model/{query_id}-{meth}-{target_pdbid}_{target_chainid}.*.pdb"))==0:
                print(f"ERROR, rank, {rank} : ./model/{query_id}-{meth}-{target_pdbid}_{target_chainid}.*.pdb")


if __name__ == '__main__':
    paths = sorted(glob.glob("./alignment/*.ali"))

    head1000paths = sorted(glob.glob("./head1000/*.head1000"))
    head1000dict = {os.path.basename(p).split(".")[0]: {line.split()[1]: (line.split()[2], i+1) for i, line in enumerate(open(p).read().splitlines()[1:]) } for p in head1000paths}
                     # query-meth                         sequence_path    ( z-score ,  rank )

    list(map(f,paths))
    # with multiprocessing.Pool(processes=18) as pool:
    #    pool.map(f,paths, chunksize=1)
