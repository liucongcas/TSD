import sys,getopt,time

def Annotation_sites():
    fr_filtered = [r1.strip().split("\t") for r1 in open(final_filtered).readlines()]
    fr_ins = [r1.strip().split("\t") for r1 in open(genefile).readlines()]
    final_list=[];final_list2=[]
    for xr in fr_filtered:
        if xr[3]=="+":
            final_list.append([xr[0], xr[1], xr[2],xr[3],xr[4],"+", "+","minus","plus"])
        elif xr[3]=="-":
            final_list.append([xr[0],xr[2],xr[1],xr[3],xr[4],"-","-","plus","minus"])
        else:
            final_list.append(get_strand(xr))
    dict_gene=split_file(fr_ins)
    for yr in final_list:
        chr_list=dict_gene[yr[0]]
        yr_l=float(yr[1])
        yr_r=float(yr[2])
        yr_l_strand=yr[5]
        yr_r_strand=yr[6]
        fac_l=False
        fac_r=False
        left_gene="none"
        right_gene="none"
        for gr in chr_list:
            gene_s=float(gr[3])
            gene_e=float(gr[4])
            gene_strand=gr[6]
            gene_name = gr[7]
            if gene_strand == yr_l_strand and yr_l >=gene_s and yr_l <=gene_e:
                left_gene=gene_name
                fac_l=True
            if gene_strand == yr_r_strand and yr_r >=gene_s and yr_r <=gene_e:
                right_gene=gene_name
                fac_r=True
            if fac_l and fac_r:
                break
        final_list2.append(yr+[left_gene,right_gene])
    final_list_in=Annotation_sites_SV(final_list2)
    lw=open(final_filtered+"_annotated",'w+')
    tt_list=[]
    tt_list.append("\t".join(["Chromosome","5_JunctionSites","3_JunctionSites","SplicingType","Surpported JunctionsReads","5_JunctionSites_strand",
                    "3_JunctionSites_strand","5_Junction_coordinates","3_Junction_coordinates","5_parent_gene","3_parent_gene","SV potential"]))
    for xx in final_list_in:
        tt_list.append("\t".join(map(str,xx)))
    lw.writelines("\n".join(tt_list)+"\n")
    lw.close()
def split_file(fr_ins):
    dict_t={}
    for ss in fr_ins:
        par=ss[0]
        if par in dict_t:
            dict_t[par].append(ss)
        else:
            dict_t[par]=[ss]
    return dict_t
def get_strand(xr):
    pp=xr[6][:2]
    tem_list=[]
    if xr[3]=="ff":
        if pp=="GT" or pp=="GC":
            tem_list=[xr[0],xr[1],xr[2],xr[3],xr[4],"+","+","minus","plus"]
        elif pp=="CT":
            tem_list=[xr[0],xr[2],xr[1],xr[3],xr[4],"-", "-","plus","minus"]
    elif xr[3]=="fr":
        if pp=="GT" or pp=="GC":
            tem_list=[xr[0],xr[1],xr[2],xr[3],xr[4],"+","-","minus","minus"]
        elif pp=="CT":
            tem_list=[xr[0],xr[2],xr[1],xr[3],xr[4], "+", "-","minus","minus"]
    elif xr[3]=="rf":
        if pp=="GT" or pp=="GC":
            tem_list=[xr[0],xr[1],xr[2],xr[3],xr[4],"-","+","plus","plus"]
        elif pp=="CT":
            tem_list=[xr[0],xr[2],xr[1],xr[3],xr[4], "-", "+","plus","plus"]
    elif xr[3]=="rr":
        if pp=="GT" or pp=="GC":
            tem_list=[xr[0],xr[1],xr[2],xr[3],xr[4],"-","-","plus","minus"]
        elif pp=="CT":
            tem_list=[xr[0],xr[2],xr[1],xr[3] ,xr[4],"+", "+","minus","plus"]
    return tem_list


def Annotation_sites_SV(final_list2):
    final_list_in=[];dict_list={};list_re=[]
    dict_list2={}
    for xr in final_list2:
        par=xr[0]+"_"+xr[2]
        dict_list[par]=xr
        dict_list2[par]=[xr[0],float(xr[2]) - binz, float(xr[2]) + binz]
    for xr in dict_list:
        binsize=dict_list2[xr][1:]
        chr_x=dict_list2[xr][0]
        for yr in final_list2:
            if chr_x == yr[0] and float(yr[1]) >=binsize[0] and float(yr[1]) <=binsize[1]:
                if binsize[0] >=float(dict_list[xr][1]) and binsize[1] <=float(yr[2]):
                    pass
                elif binsize[1] <=float(dict_list[xr][1]) and binsize[0] >=float(yr[2]):
                    pass
                else:
                    pp1 = yr + ["Potentially derived from SVs"]
                    pp2 = dict_list[xr] + ["Potentially derived from SVs"]
                    if yr not in list_re:
                        final_list_in.append(pp1)
                        list_re.append(yr)
                    if dict_list[xr] not in list_re:
                        final_list_in.append(pp2)
                        list_re.append(dict_list[xr])
    final_list_in, list_re=annotation_sites_SV_gene(final_list2,list_re,final_list_in)
    for yr3 in final_list2:
        if yr3 not in list_re:
            final_list_in.append(yr3 + ["No potential signals"])
    return  final_list_in
def annotation_sites_SV_gene(final_list2,list_re,final_list_in):
    dict_ingene = {}
    dict_g=get_nearby_gene()
    for xr2 in final_list2:
        if xr2[10] !="none" and xr2[9]!="none":
            par = xr2[0] + "_" + xr2[10]
            dict_ingene[par] = xr2
    for yr2 in final_list2:
        gene_m=yr2[9]
        gene_r=yr2[10]
        if gene_m !="none" and gene_r!="none":
            par2 = yr2[0] + "_" + yr2[9]
            if par2 in dict_ingene:
                gene_1=dict_ingene[par2][9]
                fac=True
                start_1=float(dict_g[gene_1][3])
                end_1=float(dict_g[gene_1][4])
                start_m=float(dict_g[gene_m][3])
                end_m = float(dict_g[gene_m][4])
                start_r=float(dict_g[gene_r][3])
                end_r=float(dict_g[gene_r][4])
                if (end_1 < start_m and end_m < start_r) or (start_1 > end_m and start_m > end_r):
                    fac=False
                if fac:
                    pp1 = yr2 + ["Potentially derived from SVs"]
                    pp2 = dict_ingene[par2] + ["Potentially derived from SVs"]
                    if yr2 not in list_re:
                        final_list_in.append(pp1)
                        list_re.append(yr2)
                    if dict_ingene[par2] not in list_re:
                        final_list_in.append(pp2)
                        list_re.append(dict_ingene[par2])
    return final_list_in,list_re
def get_nearby_gene():
    fr_gene = [r1.strip().split("\t") for r1 in open(genefile).readlines()]
    dict_g={}
    len_g=len(fr_gene)
    for g in range(len_g):
        par_g=fr_gene[g][7]
        dict_g[par_g]=fr_gene[g]
    return dict_g

def main(argv):
    t1 = time.time()
    global genefile
    global final_filtered
    genefile = ''
    final_filtered=''
    try:
        opts, args = getopt.getopt(argv, "hi:a:o:", ["ifile=", "afile=","odir="])
    except getopt.GetoptError:
        print'python DNAlevel_NCL_filter.py  -i <final_filtered> -a <genefile>'
        sys.exit(3)
    for opt, arg in opts:
        if opt == '-h':
            print'python DNAlevel_NCL_filter.py -i <final_filtered> -a <genefile>'
            sys.exit()
        elif opt in ("-i", "--ifile"):
            final_filtered = arg
        elif opt in ("-a", "--afile"):
            genefile = arg
    global binz
    binz=5000
    Annotation_sites()
    t2 = time.time()
    print"The time needed in anntation for TSEs is:", t2 - t1

if __name__ == "__main__":
    main(sys.argv[1:])
