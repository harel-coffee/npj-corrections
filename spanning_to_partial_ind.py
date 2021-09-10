import os
from scipy.stats import kruskal
import pysam

#MIN_AVERAGE_PHRED = 20  # (up to 37)
EDGE_NUC_NUM = 6  # recommended 3-10
MIDDLE_NUC_NUM = 6  # recommended 3-10

DIFF_NUC_ALLOWED_EDGE = 1  # always less then EDGE_NUC_NUM
DIFF_NUC_ALLOWED_MIDDLE = 1  # always less then EDGE_NUC_NUM


samfile = pysam.AlignmentFile("/groups/nshomron/tomr/projects/cffdna/runs/fam03/S03.recal.sorted.hg38.bam", "rb")
path_mom = "/groups/nshomron/tomr/projects/cffdna/runs/fam03/M03.recal.sorted.hg38.bam"
path_plasma = "/groups/nshomron/tomr/projects/cffdna/runs/fam03/S03.recal.sorted.hg38.bam"


# path_mom_simulation = "/groups/nshomron/hadasvol/projects/simulation/art/parents41/art_sim/FM41sim.srt.bam"
# path_plasma_simulation = "/groups/nshomron/hadasvol/projects/simulation/art/S41/art_sim/Ssim.chrX.srt.bam"


def create_all_possibilities():
    edge1 = "C"
    edge2 = "G"
    edge3 = "G"
    for i in range(1, EDGE_NUC_NUM):
        if (i % 3 == 1):
            edge1 += "G"
            edge2 += "G"
            edge3 += "C"
        elif (i % 3 == 2):
            edge1 += "G"
            edge2 += "C"
            edge3 += "G"
        else:
            edge1 += "C"
            edge2 += "G"
            edge3 += "G"
    return edge1, edge2, edge3


def diff_nuc(check, edge1, edge2, edge3):
    a = sum(check[i] != edge1[i] for i in range(len(check)))
    b = sum(check[i] != edge2[i] for i in range(len(check)))
    c = sum(check[i] != edge3[i] for i in range(len(check)))
    # print(f"min diff of {min(a,b,c)}")
    return min(a, b, c)


def count_cgg_left(readseq, edge1, edge2, edge3):
    cgg_loc = 0
    while (diff_nuc(readseq[cgg_loc:cgg_loc + MIDDLE_NUC_NUM], edge1, edge2, edge3) <= DIFF_NUC_ALLOWED_MIDDLE):
        cgg_loc += 1
    return cgg_loc + MIDDLE_NUC_NUM - DIFF_NUC_ALLOWED_MIDDLE



def count_cgg_right(readseq, edge1, edge2, edge3):
    cgg_loc = 0
    while (diff_nuc(readseq[-1 * MIDDLE_NUC_NUM - cgg_loc: -1 * cgg_loc], edge1, edge2,
                    edge3) <= DIFF_NUC_ALLOWED_MIDDLE):
        cgg_loc += 1
    return cgg_loc + MIDDLE_NUC_NUM - DIFF_NUC_ALLOWED_MIDDLE





def informative_vs_spanning(samfile, HEALTHY_MAX_CGG_NUC_REP, MIN_MAPPING_QUAL, MIN_LENGTH, MIN_AVERAGE_PHRED, verbose):
    bad_phred_edge = 0
    tot_cnt = 0
    used_cnt = 0
    clean_cnt = 0
    clean_weight = 0
    spanning_cnt = 0
    spanning_weight = 0
    informative_right_partial_cnt = 0
    informative_left_partial_cnt = 0
    noninformative_left_partial_cnt = 0
    noninformative_right_partial_cnt = 0
    partial_weight = 0
    # used for self check
    short_cnt = 0
    bad_map_cnt = 0
    edge1, edge2, edge3 = create_all_possibilities()
    for read in samfile.fetch("chrX", 147912050, 147912110):
        tot_cnt += 1
        # only use long reads (small ones can't be spanning so they dont add data)
        if (len(read.seq) >= MIN_LENGTH and read.mapping_quality >= MIN_MAPPING_QUAL and sum( 
                read.query_qualities[:EDGE_NUC_NUM]) / EDGE_NUC_NUM >= MIN_AVERAGE_PHRED and sum( 
            read.query_qualities[-1 * EDGE_NUC_NUM:]) / EDGE_NUC_NUM >= MIN_AVERAGE_PHRED):
            used_cnt += 1
            left_is_cgg = False
            right_is_cgg = False
            # check if left part is clean
            if diff_nuc(read.query_sequence[:EDGE_NUC_NUM], edge1, edge2, edge3) <= DIFF_NUC_ALLOWED_EDGE:
                left_is_cgg = True
            # check if right part is clean
            if diff_nuc(read.query_sequence[-1 * EDGE_NUC_NUM:], edge1, edge2, edge3) <= DIFF_NUC_ALLOWED_EDGE:
                right_is_cgg = True
            if (left_is_cgg and right_is_cgg):
                clean_cnt += 1
                clean_weight += read.mapping_quality * (sum(read.query_qualities[:EDGE_NUC_NUM]) / EDGE_NUC_NUM + sum(
                    read.query_qualities[-1 * EDGE_NUC_NUM:]) / EDGE_NUC_NUM)
                if (verbose):
                    print("clean: " + str(read.query_sequence))
            elif not (left_is_cgg or right_is_cgg):
                spanning_cnt += 1
                spanning_weight += read.mapping_quality * (
                        sum(read.query_qualities[:EDGE_NUC_NUM]) / EDGE_NUC_NUM + sum(
                    read.query_qualities[-1 * EDGE_NUC_NUM:]) / EDGE_NUC_NUM)
                if (verbose):
                    print("spanning: " + str(read.query_sequence))
            # right partial
            elif (right_is_cgg):
                cgg_count = count_cgg_right(read.query_sequence, edge1, edge2, edge3)
                if (cgg_count > HEALTHY_MAX_CGG_NUC_REP):
                    informative_right_partial_cnt += 1
                else:
                    noninformative_right_partial_cnt += 1
                if (verbose):
                    print(f"right partial ({cgg_count} nucleotides): " + str(read.query_sequence))
            elif (left_is_cgg):
                cgg_count = count_cgg_left(read.query_sequence, edge1, edge2, edge3)
                if (cgg_count > HEALTHY_MAX_CGG_NUC_REP):
                    informative_left_partial_cnt += 1
                else:
                    noninformative_left_partial_cnt += 1
                if (verbose):
                    print(f"left partial ({cgg_count} nucleotides): " + str(read.query_sequence))
        else:
            if (not len(read.seq) >= MIN_LENGTH):
                short_cnt += 1
            if (not read.mapping_quality >= MIN_MAPPING_QUAL):
                bad_map_cnt += 1
            if (not (sum(read.query_qualities[:EDGE_NUC_NUM]) / EDGE_NUC_NUM >= MIN_AVERAGE_PHRED and sum(
                    read.query_qualities[-1 * EDGE_NUC_NUM:]) / EDGE_NUC_NUM >= MIN_AVERAGE_PHRED)):
                bad_phred_edge += 1
    if (verbose):
        print(f"threw out: \n{short_cnt} for being too short \n{bad_map_cnt} for bad mapping\n{bad_phred_edge} for bad "
              f"phred on edge")
    if (verbose):
        print(
            f"found {clean_cnt} clean, {noninformative_right_partial_cnt} right noninformative, {noninformative_left_partial_cnt} left noninformative.\n\n")
    if (verbose):
        print(
        f"found {spanning_cnt} spanning, {informative_right_partial_cnt} right informative, {informative_left_partial_cnt} left informative.\n\n")
    return (used_cnt, spanning_cnt, informative_right_partial_cnt, informative_left_partial_cnt)



fam_41_path = "/groups/nshomron/hadasvol/projects/simulation/art/fragileXsim/fam41_model/plasma"
hiseq_path = "/groups/nshomron/hadasvol/projects/simulation/art/fragileXsim/HiSeqX/plasma"



def run_in_path(directory, fam, min_map_qual, min_length, min_average_phred, sick_len, verbose):
    healthy_ans = []
    sick_ans = []
    used_files = 0
    if (verbose):
        print(f"spanning, rightinformative, leftinformative, firstrep, secondrep, firstperc, secondperc, healthy")
    for filename in os.listdir(directory):
        if (filename.endswith(".bam")):
            path = os.path.join(directory, filename)
            try:
                sam_file = pysam.AlignmentFile(path, "rb")
                if (fam):
                    first_rep = int(filename[filename.find("STR")+3:filename.find("-")])
                    first_perc = int(filename[filename[:].find("-")+1:filename[7:].find(".") +7])
                    second_rep = int(filename[filename[6:].find("STR")+9:filename[12:].find("-")+12])
                else:
                    first_rep = int(filename[filename.find("STR")+3:filename.find("-")])
                    first_perc = int(filename[filename[:].find("-")+1:filename[9:].find(".") +9])
                    second_rep = int(filename[filename[10:].find("STR")+13:filename[16:].find("-")+16])
                if ((first_rep < 41 or second_rep < 41) and (first_rep > 40 or second_rep > 40) and first_perc != 50 and sick_len in (0,first_rep, second_rep)):
                    used_files += 1
                    used_cnt, sp, ri, le = informative_vs_spanning(sam_file, min(first_rep, second_rep), min_map_qual, min_length, min_average_phred, False)
                    healthy = True
                    if ((first_rep < second_rep and first_perc < 50) or (first_rep > second_rep and first_perc > 50)):
                        healthy = False
                        if ((ri + le) > 0):
                            healthy_ans.append(sp / ((ri + le)/2))
                        else:
                            healthy_ans.append(0)
                    else:
                        if ((ri + le) > 0):
                            sick_ans.append(sp / ((ri + le)/2))
                        else:
                            sick_ans.append(0)
                    if (verbose):
                        print(f"{sp},{ri},{le},{first_rep},{second_rep}, {first_perc}, {100-first_perc}, {healthy}")
            except:
                continue
        else:
            continue
    return used_files, healthy_ans, sick_ans



def try_kruskal(h, s):
    try:
       stat, p = kruskal(h,s)
    except:
        stat, p = 0, 0
    return stat, p
    

def run_with_kruskel_both_paths(map_qual, min_len, sick_len):
    if(map_qual):
        min_map_qual = 0
        print("dataset,min_map_qual,kruskal,p,used_cnt")
        while (min_map_qual < 70):
            used_cnt, h, s = run_in_path(fam_41_path, True, min_map_qual, 0, 20, 0, False)
            stat, p = try_kruskal(h,s)
            print("family," + str(min_map_qual) + "," + str(stat) + "," + str(p) + "," + str(used_cnt))
            min_map_qual += 10    
        min_map_qual = 0
        while (min_map_qual < 70):  
            used_cnt, h, s = run_in_path(hiseq_path, False, min_map_qual, 0, 20, 0, False)
            stat, p = try_kruskal(h,s)
            print("hiseq," + str(min_map_qual) + "," + str(stat) + "," + str(p) + "," + str(used_cnt))
            min_map_qual += 10
        print("-----------------------------------------------")
    if(min_len):
        min_len = 30
        print("dataset,min_len,kruskal,p,used_cnt")
        while (min_len < 160):
            used_cnt, h, s = run_in_path(fam_41_path, True, 0, min_len, 20, 0, False)
            stat, p = try_kruskal(h,s)
            print("family," + str(min_len) + "," + str(stat) + "," + str(p) + "," + str(used_cnt))
            min_len += 10    
        min_len = 30
        while (min_len < 160):
            used_cnt, h, s = run_in_path(hiseq_path, False, 0, min_len, 20, 0, False)
            stat, p = try_kruskal(h,s)
            print("hiseq," + str(min_len) + "," + str(stat) + "," + str(p) + "," + str(used_cnt))
            min_len += 10
        print("-----------------------------------------------")
    if(sick_len):
        print("dataset,sick_len,kruskal,p,used_cnt")
        sick_len = 70
        used_cnt, h, s = run_in_path(fam_41_path, False, 0, 0, 20, sick_len, False)
        stat, p = try_kruskal(h,s)
        print("family," + str(sick_len) + "," + str(stat) + "," + str(p) + "," + str(used_cnt))
        sick_len = 120
        used_cnt, h, s = run_in_path(fam_41_path, False, 0, 0, 20, sick_len, False)
        stat, p = try_kruskal(h,s)
        print("family," + str(sick_len) + "," + str(stat) + "," + str(p) + "," + str(used_cnt))
        sick_len = 200
        used_cnt, h, s = run_in_path(fam_41_path, False, 0, 0, 20, sick_len, False)
        stat, p = try_kruskal(h,s)
        print("family," + str(sick_len) + "," + str(stat) + "," + str(p) + "," + str(used_cnt))
        
        sick_len = 70
        used_cnt, h, s = run_in_path(hiseq_path, False, 0, 0, 20, sick_len, False)
        stat, p = try_kruskal(h,s)
        print("hiseq," + str(sick_len) + "," + str(stat) + "," + str(p) + "," + str(used_cnt))
        sick_len = 120
        used_cnt, h, s = run_in_path(hiseq_path, False, 0, 0, 20, sick_len, False)
        stat, p = try_kruskal(h,s)
        print("hiseq," + str(sick_len) + "," + str(stat) + "," + str(p) + "," + str(used_cnt))
        sick_len = 200
        used_cnt, h, s = run_in_path(hiseq_path, False, 0, min_len, 20, sick_len, False)
        stat, p = try_kruskal(h,s)
        print("hiseq," + str(sick_len) + "," + str(stat) + "," + str(p) + "," + str(used_cnt))


path_sim_1 = "/groups/nshomron/hadasvol/projects/simulation/art/fragileXsim/fam41_model/plasma/S41.STR200-40.STR20-60.bam"
#
sam_sim = pysam.AlignmentFile(path_sim_1, "rb")
sam_plasma = pysam.AlignmentFile(path_plasma, "rb")
# print_res(path_plasma_simulation, path_mom_simulation)
# result =create_all_possibilities()
# print(result)