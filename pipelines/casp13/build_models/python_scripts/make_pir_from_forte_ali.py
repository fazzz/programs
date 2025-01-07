#forte ali to pir
import re
import os
import Bio
import Bio.PDB
import Bio.pairwise2
import Bio.SeqIO
import multiprocessing

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



def write_pir(pir_path, query_id, target_pdbid, target_chainid, query_seq, target_seq, zscore, rank, multimer=1):
    if f"{target_pdbid.lower()}_{target_chainid}" in ["1bv6_A"]:
        print(f"ERROR: do not pass FastMMCIFParser {target_pdbid}_{target_chainid}")
        return
    if f"{target_pdbid.lower()}_{target_chainid}" in ["1ffk_B", "1phs_A", "3jcr_A"]:
        print(f"ERROR: Polypeptide.PPBuilder do not parse {target_pdbid}_{target_chainid}")
        return   
    parser = Bio.PDB.FastMMCIFParser(QUIET=True)
    try:
        structure = parser.get_structure("X", f"/home/tsukasa/work2/mmCIF/{target_pdbid[1:3].lower()}/{target_pdbid.lower()}.cif")
    except FileNotFoundError:
        try:
            structure = parser.get_structure("X", f"/home/tsukasa/work2/mmCIF_obsolete/{target_pdbid[1:3].lower()}/{target_pdbid.lower()}.cif")
        except FileNotFoundError:
            print(f"ERROR FileNotFoundError: /home/tsukasa/work2/mmCIF/{target_pdbid[1:3].lower()}/{target_pdbid.lower()}.cif")
            print(f"ERROR FileNotFoundError: /home/tsukasa/work2/mmCIF_obsolete/{target_pdbid[1:3].lower()}/{target_pdbid.lower()}.cif")
            return
    cappb = Bio.PDB.Polypeptide.CaPPBuilder()
    tmpl_seq = Bio.Seq.Seq("")
    for pp in cappb.build_peptides(structure[0][target_chainid], aa_only=False):
        # print(pp.get_sequence())
        tmpl_seq += pp.get_sequence()
    if len(tmpl_seq) == 0 and len(structure[0][target_chainid].child_list) != 0:
        print(f"ERROR: Polypeptide.PPBuilder do not parse {target_pdbid}_{target_chainid}")
        return          
    tmpl_target_alignment = Bio.pairwise2.align.localms(tmpl_seq, target_seq, 20, -10, -4, -1, one_alignment_only=True)[0] # even if one_alignment_only=True, it returns list of one alignment. so need [0]
    assert len(tmpl_target_alignment[0]) == len(tmpl_target_alignment[1])
    aligned_tmpl_seq=tmpl_target_alignment[0]
    aligned_target_seq=tmpl_target_alignment[1]
    for i, (aligned_tmpl_r, aligned_target_r) in enumerate(zip(aligned_tmpl_seq,aligned_target_seq)):
        if aligned_tmpl_r=="-" or aligned_target_r=="-":
            continue
        else:
            # index on alignment, aligned_tmpl_reidue, index on tmpl residue seq, aligned_target_reidue, index on target residue seq, 
            first_match = i, aligned_tmpl_r, i-list(aligned_tmpl_seq[:i+1]).count("-"), aligned_target_r, i-list(aligned_target_seq[:i+1]).count("-")
            break
    for i, (aligned_tmpl_r, aligned_target_r) in reversed(list(enumerate(zip(aligned_tmpl_seq,aligned_target_seq)))):
        if aligned_tmpl_r=="-" or aligned_target_r=="-":
            continue
        else:
            last_match = i, aligned_tmpl_r, i-list(aligned_tmpl_seq[:i+1]).count("-"), aligned_target_r, i-list(aligned_target_seq[:i+1]).count("-")
            break
    if len(aligned_target_seq) > len(target_seq):
        positions_where_addtional_gap_inserted_by_alignment = []
        index_on_target_seq = 0
        for i, (aligned_tmpl_r, aligned_target_r) in enumerate(zip(aligned_tmpl_seq,aligned_target_seq)):
            if len(target_seq) <= index_on_target_seq:
                # end flanking gap tekina monowo soutei, koreizyou index_on_target_seq wo huyasenai joutai
                positions_where_addtional_gap_inserted_by_alignment.append(i)
                continue
            if aligned_target_r == target_seq[index_on_target_seq]:
                # start flanking gap tekina monowo soutei
                index_on_target_seq += 1
                continue
            else:
                positions_where_addtional_gap_inserted_by_alignment.append(i)
        # print("".join([str(i%10) for i in range(len(aligned_tmpl_seq))]))
        # print(aligned_tmpl_seq)
        # print(aligned_target_seq)
        # print(positions_where_addtional_gap_inserted_by_alignment)
        l_aligned_tmpl_seq = list(aligned_tmpl_seq)
        for i in  reversed(positions_where_addtional_gap_inserted_by_alignment):
            l_aligned_tmpl_seq.pop(i)
        aligned_tmpl_seq = "".join(l_aligned_tmpl_seq)
        # print(target_seq)
        # print(aligned_tmpl_seq)
    assert len(aligned_tmpl_seq)==len(target_seq)
    # first_id = structure[0][target_chainid].child_list[first_match[2]].id[1]
    first_id = "."
    # last_id = structure[0][target_chainid].child_list[last_match[2]].id[1]
    last_id = "."
    nl = '\n'
    pir = f"""
>P1;{query_id}
sequence:{query_id}: : : : : : :0.00:0.00
{(query_seq+"/{nl}")*multimer if multimer>1 else query_seq}
*
>P1;{target_pdbid}_{target_chainid}
structure:{target_pdbid.lower()}:{first_id}:{target_chainid}:{last_id}:{target_chainid}: : :-1.00:-1.00
{(aligned_tmpl_seq+"/{nl}")*multimer if multimer>1 else aligned_tmpl_seq}
*
C;   query:{query_seq}
C;  target:{target_seq}
C;template:{aligned_tmpl_seq}

C; z-score:{zscore}
C;    rank:{rank}
"""
    # print(pir)
    with open(pir_path, "w") as f:
        f.write(pir)

def write_pir_pir_check_totally_dependes_on_modeller(pir_path, query_id, target_pdbid, target_chainid, query_seq, target_seq, multimer=1):
    nl = '\n'
    pir = textwrap.dedent(f"""
>P1;{query_id}
sequence:{query_id}: : : : : : : :
{(query_seq+"/{nl}")*multimer if multimer>1 else query_seq}
*
>P1;{target_pdbid}_{target_chainid}
structure:{target_pdbid}_{target_chainid}: : : : : : : :
{(target_seq+"/{nl}")*multimer if multimer>1 else target_seq}
*
""")
    print(pir)
    # with open(pir_path, "w") as f:
    #     f.write(pir)

def f(path):
    alis = parse_forte_ali(open(path, "r").read(), path)
    print(path, len(alis))
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
        if os.path.isfile(f"./pir/{query_id}-{meth}-{target_pdbid}_{target_chainid}.pir"):
            print(f"ERROR : ./pir/{query_id}-{meth}-{target_pdbid}_{target_chainid}.pir already exists")

        zscore, rank = head1000dict[f"{query_id}-{meth}"][target_path]
        write_pir(f"./pir/{query_id}-{meth}-{target_pdbid}_{target_chainid}.pir", f"{query_id}-{meth}-{target_pdbid}_{target_chainid}", target_pdbid, target_chainid, query_seq, target_seq, zscore, rank)



if __name__ == '__main__':
    import glob
    paths = sorted(glob.glob("./alignment/*.ali"))
    # paths = reversed(paths)

    head1000paths = sorted(glob.glob("./head1000/*.head1000"))
    head1000dict = {os.path.basename(p).split(".")[0]: {line.split()[1]: (line.split()[2], i+1) for i, line in enumerate(open(p).read().splitlines()[1:]) } for p in head1000paths}
                     # query-meth                         sequence_path    ( z-score ,  rank )

    with multiprocessing.Pool(processes=18) as pool:
       pool.map(f,paths, chunksize=1)


    # for path in paths:
    #     print(path)
    #     alis = parse_forte_ali(open(path, "r").read(), path)
    #     query_id, meth = os.path.basename(path).split("_")[0], "_".join(os.path.basename(path).split("_")[1:])
    #     for ali in alis:
    #         query_seq = ali[0]
    #         target_seq = ali[1]
    #         target_len = ali[2]
    #         target_path = ali[3]
    #         score = ali[4]
    #         if os.path.basename(target_path).startswith(("1","2","3","4","5","6","7","8","9")): # pdbid
    #             target_id = "_".join(os.path.splitext(os.path.basename(target_path))[0].split("_")[:2])
    #         elif os.path.basename(target_path).startswith("d") or os.path.basename(target_path).startswith("g"):
    #             target_id = os.path.splitext(os.path.basename(target_path))[0]
    #             print(f"cannot handle now f{target_path}")
    #             continue
    #         elif os.path.basename(target_path).startswith("PDP"):
    #             print(f"cannot handle now f{target_path}")
    #             continue
    #         else:
    #             assert False
    #         write_pir(f"./pir/{meth}.{query_id}.{target_id}.pir", query_id, target_id, query_seq, target_seq)
