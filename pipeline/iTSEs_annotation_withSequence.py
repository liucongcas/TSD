import sys,getopt,time

def Annotation_sites():
    fr_filtered = [r1.strip().split("\t") for r1 in open(final_filtered).readlines()]
    fr_ins = [r1.strip().split("\t") for r1 in open(genefile).readlines()]
    final_list=[];final_list2=[]
    for xr in fr_filtered:
        if xr[3]=="+":
            final_list.append([xr[0], xr[1], xr[2],xr[3],xr[4],"+", "+","minus","plus",xr[5],xr[8]])
        elif xr[3]=="-":
            x1 = reverse(xr[8])
            x2 = reverse(xr[5])
            final_list.append([xr[0],xr[2],xr[1],xr[3],xr[4],"-","-","plus","minus",x1,x2])
        else:
            final_list.append(get_strand(xr))
    dict_gene=split_file(fr_ins)
    for yr in final_list:
        #print yr
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
    lw=open(final_filtered+"_annotated.txt",'w+')
    tt_list=[]
    tt_list.append("\t".join(["Chromosome","5_JunctionSites","3_JunctionSites","SplicingType","Surpported JunctionsReads","5_JunctionSites_strand",
                    "3_JunctionSites_strand","5_Junction_coordinates","3_Junction_coordinates","5_parent_gene","3_parent_gene"]))
    for xx in final_list2:
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
            tem_list=[xr[0],xr[1],xr[2],xr[3],xr[4],"+","+","minus","plus",xr[5],xr[8]]
        elif pp=="CT":
            x1=reverse(xr[8])
            x2=reverse(xr[5])
            tem_list=[xr[0],xr[2],xr[1],xr[3],xr[4],"-", "-","plus","minus",x1,x2]
    elif xr[3]=="fr":
        if pp=="GT" or pp=="GC":
            tem_list=[xr[0],xr[1],xr[2],xr[3],xr[4],"+","-","minus","minus",xr[5],xr[8]]
        elif pp=="CT":
            x1 = reverse(xr[8])
            x2 = reverse(xr[5])
            tem_list=[xr[0],xr[2],xr[1],xr[3],xr[4], "+", "-","minus","minus",x1,x2]
    elif xr[3]=="rf":
        if pp=="GT" or pp=="GC":
            tem_list=[xr[0],xr[1],xr[2],xr[3],xr[4],"-","+","plus","plus",xr[5],xr[8]]
        elif pp=="CT":
            x1 = reverse(xr[8])
            x2 = reverse(xr[5])
            tem_list=[xr[0],xr[2],xr[1],xr[3],xr[4], "-", "+","plus","plus",x1,x2]
    elif xr[3]=="rr":
        if pp=="GT" or pp=="GC":
            tem_list=[xr[0],xr[1],xr[2],xr[3],xr[4],"-","-","plus","minus",xr[5],xr[8]]
        elif pp=="CT":
            x1 = reverse(xr[8])
            x2 = reverse(xr[5])
            tem_list=[xr[0],xr[2],xr[1],xr[3] ,xr[4],"+", "+","minus","plus",x1,x2]
    return tem_list

def reverse(line):
    transline = line[::-1].replace('A', 't').replace('T', 'a').replace('G', 'c').replace('C', 'g').upper()
    return transline


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
    Annotation_sites()
    t2 = time.time()
    print"The time needed in anntation for TSEs is:", t2 - t1

if __name__ == "__main__":
    main(sys.argv[1:])