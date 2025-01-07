import os
import glob
import parse
import Bio
import Bio.pairwise2




pir_parse = """
>P1;{query_id}
sequence:{query_id}: : : : : : :0.00:0.00
{query_seq}
*
>P1;{target_pdbid}_{target_chainid}
structure:{target_pdbid_lower}:{first_id}:{target_chainid}:{last_id}:{target_chainid}: : :-1.00:-1.00
{aligned_tmpl_seq}
*
C;   query:{query_seq}
C;  target:{target_seq}
C;template:{aligned_tmpl_seq}

C; z-score:{zscore}
C;    rank:{rank}
"""
#pirs = ["pir/T0984_N39C48-f_hhpx_sspx-5e1j_A.pir","pir/T0984_N39C48-f_db_sspx3-6c96_A.pir"]
#pirs = sorted(glob.glob("pir/*1zat*.pir")+glob.glob("pir/*1ZAT*.pir")+glob.glob("pir/*5uwv*.pir")+glob.glob("pir/*5UWV*.pir"))
pirs = ["pir/T1023s3-f_db_hhp-5k0y_S.pir"]
pirs = sorted(pirs+glob.glob("pir/*2pmd*.pir")+glob.glob("pir/*2d74*.pir")+glob.glob("pir/*4m53*.pir")+glob.glob("pir/*1kjz*.pir")+glob.glob("pir/*3j81*.pir"))
#pirs = sorted(pirs+glob.glob("template/*.pir"))

#pirs = [ "pir/"+p for p in [
#"T1012-f_hhpx_sspx3-4fd5_A.pir","T1012-f_ssp_ssp-4fd4_A.pir","T1012-f_ssp_ssp3-4fd4_A.pir","T1012-f_sspx_sspx3-4fd4_A.pir","T1012-f_hhpx_px-4fd4_A.pir","T1012-f_ssp_pb-4fd4_A.pir","T1012-f_ssp_px-4fd4_A.pir","T1012-f_pdppx_ssp3-4M85_A.pir","T1012-f_ssp_ssp3-4fd5_A.pir","T1012-f_pdppx_pb-4M85_A.pir","T1012-f_sspx_pb-4fd4_A.pir","T1012-f_pdpp_hhx-4M85_A.pir","T1012-f_sspx_ssp-4fd4_A.pir","T1012-f_ssp_sspx3-4fd4_A.pir","T1012-f_hhpx_db-4fd4_A.pir","T1012-f_hhpx_sspx3-5k9n_A.pir","T1012-s_sspx_sspx3-2qec_A.pir","T1012-f_pdppx_sspx3-4M85_A.pir","T1012-s_ssp_sspx-2qec_A.pir","T1012-f_pdpp_ssp3-5K9N_A.pir","T1012-f_ssp_sspx-4fd7_A.pir","T1012-f_ssp_sspx3-4fd7_A.pir","T1012-s_ssp_px-4fd4_A.pir","T1012-s_dbx_sspx3-2qec_A.pir","T1012-f_ssp_ssp-4fd5_A.pir","T1012-f_sspx_sspx3-2qec_A.pir","T1012-s_ssp_sspx3-4fd4_A.pir","T1012-f_ssp_ssp3-4fd7_A.pir","T1012-s_ssp_ssp3-2qec_A.pir","T1012-f_pdpp_px-4M85_A.pir","T1012-s_pdpp_hhx-3ZJ0_A.pir","T1012-f_pdppx_ssp3-5K9N_A.pir","T1012-s_ssp_sspx3-2qec_A.pir","T1012-f_pdpp_sspx3-4M85_A.pir","T1012-f_pdpp_pb-4M85_A.pir","T1012-f_sspx_sspx3-4fd7_A.pir","T1012-f_hhpx_db-4fd5_A.pir","T1012-f_pdpp_ssp-4FD5_A.pir","T1012-f_hhpx_pb-4bmh_A.pir","T1012-s_sspx_ssp-2qec_A.pir","T1012-f_ssp_hhpx-4fd4_A.pir","T1012-f_hhpx_dbx-5k9n_A.pir","T1012-s_pdpp_ssp3-5K9N_A.pir","T1012-f_dbx_ssp-2qec_A.pir","T1012-f_sspx_ssp3-4fd4_A.pir","T1012-f_db_ssp-1u6m_A.pir","T1012-f_dbx_ssp3-2qec_A.pir","T1012-f_pdpp_hhpx-4M85_A.pir","T1012-f_ssp_ssp-2qec_A.pir","T1012-f_sspx_sspx-2qec_A.pir","T1012-f_hhpx_hhp-4m85_A.pir","T1012-f_hhpx_sspx-2qec_A.pir","T1012-f_sspx_ssp-2qec_A.pir","T1012-s_dbx_ssp3-2qec_A.pir","T1012-f_hhpx_sspx3-4fd4_A.pir","T1012-f_sspx_px-4fd4_A.pir","T1012-f_sspx_ssp3-4fd7_A.pir","T1012-f_hhpx_ssp-4m85_A.pir","T1012-s_pdppx_px-4M85_A.pir","T1012-s_sspx_sspx-2qec_A.pir","T1012-f_sspx_ssp3-4bmh_A.pir","T1012-f_ssp_sspx-4fd4_A.pir","T1012-s_pdpp_ssp-5K9N_A.pir","T1012-f_pdppx_px-4M85_A.pir","T1012-s_hhpx_ssp3-2qec_A.pir","T1012-s_ssp_px-2qec_A.pir","T1012-s_pdpp_pb-4M85_A.pir","T1012-f_hhpx_hhpx-4m85_A.pir","T1012-f_hhpx_hhx-4m85_A.pir","T1012-f_dbx_sspx3-2qec_A.pir","T1012-f_dbx_sspx-2qec_A.pir","T1012-f_hhpx_hhpx-4fd4_A.pir","T1012-s_pdppx_ssp3-5K9N_A.pir","T1012-s_sspx_pb-4fd4_A.pir","T1012-s_ssp_sspx3-4fd7_A.pir","T1012-f_hhpx_hhp-4fd4_A.pir","T1012-f_ssp_hhx-4fd4_A.pir","T1012-s_hhpx_sspx-2qec_A.pir","T1012-f_hhpx_hhp-3qb8_A.pir","T1012-s_ssp_hhp-2qec_A.pir","T1012-s_pdppx_pb-4M85_A.pir","T1012-f_hhpx_ssp3-4fd7_A.pir","T1012-f_hhpx_ssp3-2qec_A.pir","T1012-s_sspx_ssp3-2qec_A.pir","T1012-f_pdppx_hhx-1u6m_A.pir","T1012-s_sspx_ssp3-4fd4_A.pir","T1012-f_hhpx_px-4fd7_A.pir","T1012-f_hhpx_hhpx-3zj0_A.pir","T1012-f_db_sspx3-2qec_A.pir","T1012-s_ssp_ssp-2qec_A.pir","T1012-f_ssp_pb-2qec_A.pir","T1012-f_ssp_sspx3-2qec_A.pir","T1012-s_ssp_hhpx-2qec_A.pir","T1012-f_hhpx_hhx-4fd4_A.pir","T1012-f_sspx_px-4fd7_A.pir","T1012-f_hhpx_ssp3-5k9n_A.pir","T1012-f_ssp_ssp3-2qec_A.pir","T1012-s_ssp_hhx-2qec_A.pir","T1012-s_ssp_pb-4fd4_A.pir","T1012-s_ssp_pb-2qec_A.pir","T1012-f_sspx_ssp3-4fd5_A.pir","T1012-f_ssp_pb-4fd7_A.pir","T1012-f_pdppx_ssp3-3QB8_A.pir","T1012-f_hhpx_pb-4m85_A.pir","T1012-s_ssp_hhx-4fd5_A.pir","T1012-f_hhpx_px-3zj0_A.pir","T1012-f_sspx_ssp3-2qec_A.pir","T1012-f_db_sspx-2qec_A.pir","T1012-s_ssp_sspx3-4fd5_A.pir","T1012-s_hhpx_ssp3-4fd4_A.pir","T1012-s_sspx_hhx-2qec_A.pir","T1012-s_sspx_hhp-2qec_A.pir","T1012-s_pdppx_sspx3-4FD4_A.pir","T1012-f_pdppx_sspx-5K9N_A.pir","T1012-s_hhpx_ssp-2qec_A.pir","T1012-s_hhpx_ssp3-5k9n_A.pir","T1012-f_sspx_hhx-4fd4_A.pir","T1012-f_pdpp_ssp3-4FD5_A.pir","T1012-f_hhpx_db-1u6m_A.pir","T1012-s_ssp_ssp3-4fd7_A.pir","T1012-f_sspx_sspx-5k9n_A.pir","T1012-s_hhpx_ssp3-4m85_A.pir","T1012-f_ssp_sspx-2qec_A.pir","T1012-s_ssp_dbx-2qec_A.pir","T1012-s_ssp_ssp-4fd5_A.pir","T1012-f_sspx_px-2qec_A.pir","T1012-s_ssp_ssp3-4fd5_A.pir","T1012-f_hhpx_ssp3-3qb8_A.pir","T1012-f_ssp_px-4fd7_A.pir","T1012-f_sspx_ssp-4fd5_A.pir","T1012-f_hhpx_sspx3-4m85_A.pir","T1012-f_sspx_dbx-4fd5_A.pir","T1012-s_ssp_ssp3-4fd4_A.pir","T1012-s_pdpp_px-4M85_A.pir","T1012-f_dbx_pb-2qec_A.pir","T1012-f_ssp_hhp-4fd5_A.pir","T1012-f_sspx_hhp-3qb8_A.pir","T1012-s_ssp_hhpx-4fd5_A.pir","T1012-f_pdpp_ssp-5K9N_A.pir","T1012-s_sspx_ssp3-4fd7_A.pir","T1012-f_hhpx_dbx-4fd4_A.pir","T1012-f_sspx_hhx-2qec_A.pir","T1012-f_pdppx_ssp-4M85_A.pir","T1012-s_hhpx_sspx3-2qec_A.pir","T1012-f_hhpx_ssp3-4fd4_A.pir","T1012-s_hhpx_hhx-2qec_A.pir","T1012-f_hhpx_sspx3-3te4_A.pir","T1012-f_sspx_hhpx-2qec_A.pir","T1012-f_sspx_sspx-4fd4_A.pir","T1012-f_hhpx_ssp-2qec_A.pir","T1012-f_ssp_hhpx-3te4_A.pir","T1012-f_ssp_dbx-4fd4_A.pir","T1012-s_pdppx_hhpx-4FD5_A.pir","T1012-f_ssp_hhp-2qec_A.pir","T1012-s_dbx_px-2qec_A.pir","T1012-f_sspx_pb-4fd7_A.pir","T1012-f_dbx_px-1u6m_A.pir","T1012-f_hhpx_pb-2qec_A.pir","T1012-f_sspx_pb-2qec_A.pir","T1012-f_hhpx_hhpx-2qec_A.pir","T1012-f_sspx_ssp-3qb8_A.pir","T1012-f_hhpx_dbx-4m85_A.pir","T1012-f_pdppx_pb-5K9N_A.pir","T1012-s_pdpp_hhx-4M85_A.pir","T1012-f_ssp_hhp-3te4_A.pir","T1012-f_dbx_hhpx-1u6m_A.pir","T1012-f_pdppx_ssp-5K9N_A.pir","T1012-s_hhpx_hhx-3qb8_A.pir","T1012-f_sspx_px-4fd5_A.pir","T1012-f_hhpx_px-2qec_A.pir","T1012-s_sspx_ssp3-4fd5_A.pir","T1012-f_hhpx_hhpx-3qb8_A.pir","T1012-f_pdppx_sspx3-4FD4_A.pir","T1012-s_sspx_hhpx-2qec_A.pir","T1012-f_pdppx_sspx3-5K9N_A.pir","T1012-f_pdppx_hhpx-4M85_A.pir","T1012-f_pdppx_hhpx-4m85_A.pir","T1012-f_hhpx_pb-4fd5_A.pir","T1012-s_sspx_pb-2qec_A.pir","T1012-f_hhpx_dbx-3qb8_A.pir","T1012-f_pdpp_ssp-4FD4_A.pir","T1012-s_hhpx_px-2qec_A.pir","T1012-f_ssp_sspx-4fd5_A.pir","T1012-f_sspx_dbx-2qec_A.pir","T1012-f_pdppx_hhp-4FD5_A.pir","T1012-s_pdpp_ssp3-4FD7_A.pir","T1012-f_sspx_ssp-4fd7_A.pir","T1012-f_sspx_ssp-4bmh_A.pir","T1012-f_pdpp_db-4FD4_A.pir","T1012-f_pdpp_ssp3-4M85_A.pir","T1012-f_ssp_hhp-4fd4_A.pir","T1012-s_hhpx_hhpx-2qec_A.pir","T1012-f_sspx_px-3te4_A.pir","T1012-s_ssp_hhx-4fd4_A.pir",]
#]


for fp1 in pirs:
    for fp2 in pirs:
        if fp1 == fp2 : continue
        if fp1.split("-")[-1] == fp2.split("-")[-1]: continue
        if fp1.split("-")[-1] != "5k0y_S.pir": continue
        print(fp1, fp2)
        if os.path.basename(fp1).split("-")[-1] == os.path.basename(fp2).split("-")[-1]: continue
        # fp1 = "pir/T0949_N23-f_hhpx_hhpx-1cc3_A.pir"
        # fp2 = "pir/T0949_N23-s_hhpx_hhx-1ezl_A.pir"

        r1 = parse.parse(pir_parse,open(fp1).read())
        r2 = parse.parse(pir_parse,open(fp2).read())

        aligned_qs1,aligned_qs2,_,_,_ = Bio.pairwise2.align.globalms(r1["query_seq"], r2["query_seq"], 10, -1, -.5, -.1, one_alignment_only=True)[0]


        if len(aligned_qs1) > len(r1["query_seq"]):
            positions_where_addtional_gap_inserted_by_alignment = []
            index_on_qs = 0
            for i, aligned_qs1_r in enumerate(aligned_qs1):
                if len(r1["query_seq"]) <= index_on_qs:
                    # end flanking gap tekina monowo soutei, koreizyou index_on_qs wo huyasenai joutai
                    positions_where_addtional_gap_inserted_by_alignment.append(i)
                    continue
                if aligned_qs1_r == r1["query_seq"][index_on_qs]:
                    # start flanking gap tekina monowo soutei
                    index_on_qs += 1
                    continue
                else:
                    positions_where_addtional_gap_inserted_by_alignment.append(i)
            l_aligned_tmpl_seq = list(r1["aligned_tmpl_seq"])
            for i in  positions_where_addtional_gap_inserted_by_alignment:
                l_aligned_tmpl_seq.insert(i, "-")
            aligned_tmpl_seq1 = "".join(l_aligned_tmpl_seq)
        else:
            aligned_tmpl_seq1 = r1["aligned_tmpl_seq"]


        if len(aligned_qs2) > len(r2["query_seq"]):
            positions_where_addtional_gap_inserted_by_alignment = []
            index_on_qs = 0
            for i, aligned_qs2_r in enumerate(aligned_qs2):
                if len(r2["query_seq"]) <= index_on_qs:
                    # end flanking gap tekina monowo soutei, koreizyou index_on_qs wo huyasenai joutai
                    positions_where_addtional_gap_inserted_by_alignment.append(i)
                    continue
                if aligned_qs2_r == r2["query_seq"][index_on_qs]:
                    # start flanking gap tekina monowo soutei
                    index_on_qs += 1
                    continue
                else:
                    positions_where_addtional_gap_inserted_by_alignment.append(i)
            l_aligned_tmpl_seq = list(r2["aligned_tmpl_seq"])
            for i in  positions_where_addtional_gap_inserted_by_alignment:
                l_aligned_tmpl_seq.insert(i, "-")
            aligned_tmpl_seq2 = "".join(l_aligned_tmpl_seq)
        else:
            aligned_tmpl_seq2  = r2["aligned_tmpl_seq"]

        aligned = aligned_qs1 if aligned_qs1.count("-") < aligned_qs2.count("-") else aligned_qs2
        print(fp1,fp2, aligned)

        query_id_1 = os.path.splitext(os.path.basename(fp1))[0]
        query_id_2 = os.path.splitext(os.path.basename(fp2))[0]

        query_id = "-".join([query_id_1.split("-")[0], *query_id_1.split("-")[1:],*query_id_2.split("-")[1:]])

        pir = f"""
>P1;{query_id}
sequence:{query_id}: : : : : : :0.00:0.00
{aligned}
*
>P1;{r1["target_pdbid"]}_{r1["target_chainid"]}
structure:{r1["target_pdbid"].lower()}:.:.:.:.: : :-1.00:-1.00
{aligned_tmpl_seq1}
*
>P1;{r2["target_pdbid"]}_{r2["target_chainid"]}
structure:{r2["target_pdbid"].lower()}:.:.:.:.: : :-1.00:-1.00
{aligned_tmpl_seq2}
*
C;    query:{aligned}
C;template1:{aligned_tmpl_seq1}
C;template2:{aligned_tmpl_seq2}

C; z-score:{r1["zscore"]}
C;    rank:{r1["rank"]}
C; z-score:{r2["zscore"]}
C;    rank:{r2["rank"]}
"""

        pir_path="pir2/"+query_id+".pir"
        with open(pir_path, "w") as f:
            f.write(pir)
