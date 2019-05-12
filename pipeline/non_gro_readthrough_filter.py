import sys,getopt,time,re
def filterRT_Fus(fr,list_filter):
    dict_total={}
    list_filter_fus=[]
    for i in  range(len(fr)):
        par="_".join(fr[i][:3])
        dict_total[par]=fr[i]
    for j in range(len(fr)):
        if abs(float(fr[j][2])-float(fr[j][1]))<=adj_level:
            list_filter_fus.append(fr[j])
    split_gene=neigborGene()
    list_filter_t1=filter_neigborGene_adjacent(list_filter_fus,split_gene)
    list_filter_final=list_filter_t1
    for y in range(len(list_filter_final)):
        par2="_".join(map(str,list_filter_final[y][:3]))
        del dict_total[par2]
    for k in dict_total.keys():
        s1=dict_total[k]
        list_filter.append(s1)
    return list_filter

def filter_neigborGene_adjacent(list_filter_fus,split_gene):
    readth_list=[]
    for x1 in range(len(list_filter_fus)):
        pp=list_filter_fus[x1][0]
        single_list=split_gene[pp]
        TSE_l=min(float(list_filter_fus[x1][1]),float(list_filter_fus[x1][2]))
        TSE_r=max(float(list_filter_fus[x1][1]),float(list_filter_fus[x1][2]))
        strand_l=list_filter_fus[x1][5]
        strand_r=list_filter_fus[x1][6]
        left_gene=False
        right_gene=False
        between_gene=[]
        for x2  in range(len(single_list)):
            gene_l=min(float(single_list[x2][3]),float(single_list[x2][4]))
            gene_r=max(float(single_list[x2][3]),float(single_list[x2][4]))
            gene_strand=single_list[x2][6]
            if strand_l==strand_r:  #when the strand of two sites are same, find the interval gene in the same strand
                if strand_l==gene_strand:
                    if min(TSE_r,gene_r)-max(TSE_l,gene_l)>0:
                        between_gene.append(single_list[x2])
                    if TSE_l >= gene_l and TSE_l <= gene_r :
                        left_gene = True
                    if TSE_r >= gene_l and TSE_r <= gene_r:
                        right_gene = True
            else:
                if min(TSE_r, gene_r) - max(TSE_l, gene_l) > 0:
                    between_gene.append(single_list[x2])
                if TSE_l >= gene_l and TSE_l <= gene_r:
                    left_gene = True
                if TSE_r >= gene_l and TSE_r <= gene_r:
                    right_gene = True
        if len(between_gene)<=1: # the gene number before uniq
            readth_list.append(list_filter_fus[x1] + ['none'])
        else:
            between_gene_uniq=[between_gene[0]]
            gene_side_o=float(between_gene[0][4])
            list_uniqGene=[between_gene[0][7]]
            for r1 in  range(1,len(between_gene)):
                gene_side=float(between_gene[r1][4])
                if between_gene[r1][7] not in list_uniqGene:
                    list_uniqGene.append(between_gene[r1][7])
                    if gene_side > gene_side_o:
                        between_gene_uniq.append(between_gene[r1])
                        gene_side_o=float(between_gene[r1][4])
            n_b=len(between_gene_uniq)
            if  n_b < 1 : ##the gene nubere after uniq
                readth_list.append(list_filter_fus[x1] + ['none'])
            else:
                if left_gene and right_gene:
                    if n_b < 4:
                        readth_list.append(list_filter_fus[x1]+['gene-gene']+list_uniqGene)
                elif left_gene or right_gene:
                    if n_b < 3:
                        readth_list.append(list_filter_fus[x1]+['gene-none']+list_uniqGene)
                else:
                    if n_b <2:
                        readth_list.append(list_filter_fus[x1]+['none-none']+list_uniqGene)
    return readth_list



def neigborGene():
    split_gene={}
    fr_gene = [i.strip().split("\t") for i in open(genefile).readlines()]
    for s in range(len(fr_gene)):
        p1=fr_gene[s][0]
        if split_gene.has_key(p1):
            split_gene[p1].append(fr_gene[s])
        else:
            split_gene[p1]=[fr_gene[s]]
    return split_gene

def filterReadthrough():
    fr = [i.strip().split("\t") for i in open(targetfile).readlines()]
    list_filter = []
    for i in range(1,len(fr)):
        xr = fr[i]
        if xr[5] == xr[6] == "+":
            if float(xr[1]) < float(xr[2]):
                if abs(float(xr[2]) - float(xr[1])) > 1000000:
                    list_filter.append(xr)
            else:
                list_filter = filterRT_Fus([xr], list_filter)
        elif xr[5] == xr[6] == "-":
            if float(xr[1]) > float(xr[2]):
                if abs(float(xr[1]) - float(xr[2])) > 1000000:
                    list_filter.append(xr)
            else:
                list_filter = filterRT_Fus([xr], list_filter)
        else:
            list_filter = filterRT_Fus([xr], list_filter)
    list_filter2 = []
    list_filter2.append("\t".join(["Chromosome", "5_JunctionSites", "3_JunctionSites", "SplicingType", "Surpported JunctionsReads", "5_JunctionSites_strand",
         "3_JunctionSites_strand", "5_Junction_coordinates", "3_Junction_coordinates", "5_parent_gene", "3_parent_gene", "SV potential"]))
    for k in list_filter:
        s1 = "\t".join(k)
        list_filter2.append(s1)
    if len(list_filter2) > 0:
        lw1 = open(prefix + "_nonGRO_RTfiltered.txt", 'w+')
        lw1.writelines("\n".join(list_filter2) + "\n")
        lw1.close()


def main(argv):
    global targetfile
    global genefile
    genefile=''
    targetfile=''
    try:
        opts, args = getopt.getopt(argv, "hi:j:", ["ifile=","jfile="])
    except getopt.GetoptError:
        print 'python non_gro_readthrough_filter.py -i <typeII erro filtered file>  -j  <genefile> '
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print 'python non_gro_readthrough_filter.py -i <typeII erro filtered file>  -j  <genefile> '
            sys.exit()
        elif opt in ("-i", "--ifile"):
            targetfile = arg
        elif opt in ("-j", "-jfile"):
            genefile = arg
    t1=time.time()
    global adj_level
    adj_level=200000
    global prefix
    prefix = targetfile
    filterReadthrough()
    t2=time.time()
    print "The time needed in non_gro_readthrough_readthrough_filter for fusions is:", t2 - t1
if __name__ == "__main__":
    main(sys.argv[1:])
