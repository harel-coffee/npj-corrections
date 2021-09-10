module load python/miniconda3-4.5.12-pytorch
/share/apps/python/miniconda3-4.5.12-pytorch/bin/python




MIN_LENGTH = 140
MIN_AVERAGE_PHRED = 20 #(up to 37)
MIN_MAPPING_QUAL = 20 #(up to 60, pysam recommends 50) 
EDGE_NUC_NUM = 6 #recommended 3-10
DIFF_NUC_ALLOWED = 1 #always less then EDGE_NUC_NUM



import pysam
samfile = pysam.AlignmentFile("/groups/nshomron/tomr/projects/cffdna/runs/fam03/S03.recal.sorted.hg38.bam", "rb")
path_mom = "/groups/nshomron/tomr/projects/cffdna/runs/fam03/M03.recal.sorted.hg38.bam"
path_plasma = "/groups/nshomron/tomr/projects/cffdna/runs/fam03/S03.recal.sorted.hg38.bam"

path_mom_simulation = "/groups/nshomron/hadasvol/projects/simulation/art/parents41/art_sim/FM41sim.srt.bam"
path_plasma_simulation = "/groups/nshomron/hadasvol/projects/simulation/art/S41/art_sim/Ssim.chrX.srt.bam"



def create_all_possibilities():
    edge1 = "C"
    edge2 = "G"
    edge3 = "G"
    for i in range(1,EDGE_NUC_NUM):
        if (i%3 == 1):
            edge1 += "G"
            edge2 += "G"
            edge3 += "C"
        elif (i%3 == 2):
            edge1 += "G"
            edge2 += "C"
            edge3 += "G"
        else:
            edge1 += "C"
            edge2 += "G"
            edge3 += "G"
    return edge1, edge2, edge3

  
def diff_nuc(check,edge1, edge2, edge3):
    a =  sum( check[i] != edge1[i] for i in range(len(check)) )
    b =  sum( check[i] != edge2[i] for i in range(len(check)) )
    c =  sum( check[i] != edge3[i] for i in range(len(check)) )
    #print(f"min diff of {min(a,b,c)}")
    return min(a,b,c)



def clean_to_span(samfile): 
    tot_cnt = 0
    used_cnt = 0
    clean_cnt = 0
    clean_weight = 0
    spanning_cnt = 0
    spanning_weight = 0
    partial_cnt = 0
    partial_weight = 0
    
    #used for self check
    short_cnt = 0
    bad_map_cnt = 0
    bad_phred_edge = 0
    
    edge1, edge2, edge3 = create_all_possibilities()
    
    for read in samfile.fetch("chrX",147912050, 147912110):
        tot_cnt += 1
        #only use long reads (small ones can't be spanning so they dont add data)
        if (len(read.seq) >= MIN_LENGTH and read.mapping_quality >= MIN_MAPPING_QUAL and sum(read.query_qualities[:EDGE_NUC_NUM]) / EDGE_NUC_NUM >= MIN_AVERAGE_PHRED and sum(read.query_qualities[-1 * EDGE_NUC_NUM:]) / EDGE_NUC_NUM >= MIN_AVERAGE_PHRED):
            used_cnt +=1
            clean = 0
            #check if left part is clean
            if (diff_nuc(read.query_sequence[:EDGE_NUC_NUM], edge1, edge2, edge3) <= DIFF_NUC_ALLOWED):
                clean+=1
            #check if right part is clean
            if (diff_nuc(read.query_sequence[-1 * EDGE_NUC_NUM:], edge1, edge2, edge3) <= DIFF_NUC_ALLOWED):
                clean+=1
            if (clean == 2):
                clean_cnt += 1
                clean_weight +=  read.mapping_quality * (sum(read.query_qualities[:EDGE_NUC_NUM]) / EDGE_NUC_NUM + sum(read.query_qualities[-1 * EDGE_NUC_NUM:]) / EDGE_NUC_NUM)
                print("clean: " + str(read.query_sequence))
            elif (clean == 0): 
                spanning_cnt += 1
                spanning_weight += read.mapping_quality * (sum(read.query_qualities[:EDGE_NUC_NUM]) / EDGE_NUC_NUM + sum(read.query_qualities[-1 * EDGE_NUC_NUM:]) / EDGE_NUC_NUM)
                print("spanning: " + str(read.query_sequence))
            else:
                partial_cnt += 1
                partial_weight +=  read.mapping_quality * (sum(read.query_qualities[:EDGE_NUC_NUM])  / EDGE_NUC_NUM  + sum(read.query_qualities[-1 * EDGE_NUC_NUM:]) / EDGE_NUC_NUM)
                print("partial: " + str(read.query_sequence))
        else:
            if (not len(read.seq) >= MIN_LENGTH ):
                short_cnt +=1
            if (not read.mapping_quality >= MIN_MAPPING_QUAL):
                bad_map_cnt +=1
            if (not (sum(read.query_qualities[:EDGE_NUC_NUM]) / EDGE_NUC_NUM >= MIN_AVERAGE_PHRED and sum(read.query_qualities[-1 * EDGE_NUC_NUM:]) / EDGE_NUC_NUM >= MIN_AVERAGE_PHRED)):
                bad_phred_edge += 1
    print(f"threw out: \n{short_cnt} for being too short \n{bad_map_cnt} for bad mapping\n{bad_phred_edge} for bad phred on edge")
    print(f"found {clean_cnt} clean, {spanning_cnt} spanning and {partial_cnt} partial.")
    if (partial_cnt != 0):
        return (spanning_cnt / partial_cnt,spanning_weight / partial_weight, used_cnt, tot_cnt)
    return (1,clean_cnt + spanning_cnt, used_cnt, tot_cnt)

    


def print_res(path_plasma, path_mom):
    sam_mom = pysam.AlignmentFile(path_mom, "rb")
    sam_plasma = pysam.AlignmentFile(path_plasma, "rb")
    (a,b,c,d) = clean_to_span(sam_mom)
    (e,f,g,h) = clean_to_span(sam_plasma)
    print(f"mom's clean to spanning ratio is {b} ({a} when weighted) out of {c} good samples out of {d} samples).\n plasma's clean to spanning ratio is {f} ({e} when weighted) out of {g} good samples out of {h} samples).")
    print(f"Parameter is {(a/e)} .The bigger it is (when bigger than 1), the bigger the chance of healthy child )")
    return a/e


print_res(path_plasma, path_mom)
print_res(path_plasma_simulation, path_mom_simulation)