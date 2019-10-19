import sys,getopt,time

def TSnearSV():
    fr_filtered = [r1.strip().split("\t") for r1 in open(GRO_filteredfile).readlines()]
    fr_ins = [r1.strip().split("\t") for r1 in open(SVs_annotation_vcffile).readlines()]
    dict_ins = {};list_near2 = []
    fr_ins2=fr_ins
    for q2 in range(len(fr_ins2)):
        par = fr_ins[q2][0]
        if par in dict_ins:
            dict_ins[par].append(fr_ins[q2])
        else:
            dict_ins[par] = [fr_ins[q2]]
    for i in range(1,len(fr_filtered)):
        s_chr = fr_filtered[i][0]
        if dict_ins.has_key(s_chr):
            list_tmp1 = dict_ins[s_chr]
            s_l = float(fr_filtered[i][1])
            s_r = float(fr_filtered[i][2])
            for j in range(len(list_tmp1)):
                t_l = float(list_tmp1[j][1])
                t_r = float(list_tmp1[j][2])
                if abs(s_l - t_l) < SV_dis and abs(s_r - t_r) < SV_dis:
                        list_near2.append(fr_filtered[i])
                        break
    list_final2 = [i for i in fr_filtered if i not in list_near2]
    list_final3 = sorted(list_final2)
    tt_list=[]
    for iterm1 in list_final3:
        s1 = '\t'.join(iterm1)
        tt_list.append(s1)
    if len(tt_list)>0:
        lw1 = open(outputdir+"/final_iTSEs_removedSVs.txt", 'w+')
        lw1.writelines("\n".join(tt_list) + "\n")
        lw1.close()



def main(argv):
    t1 = time.time()
    global outputdir
    global SVs_annotation_vcffile
    global GRO_filteredfile
    global SV_dis
    SVs_annotation_vcffile = ''
    GRO_filteredfile=''
    SV_dis=''
    outputdir=''
    try:
        opts, args = getopt.getopt(argv, "hi:a:s:o:", ["ifile=", "afile=","sfile=","odir="])
    except getopt.GetoptError:
        print'python DNAlevel_NCL_filter.py  -i <GRO_filteredfile> -a <SVs_annotation_vcffile> -s <SV_dis> -o <outputdir>'
        sys.exit(4)
    for opt, arg in opts:
        if opt == '-h':
            print'python DNAlevel_NCL_filter.py -i <GRO_filteredfile> -a <SVs_annotation_vcffile> -s <SV_dis> -o <outputdir>'
            sys.exit()
        elif opt in ("-i", "--ifile"):
            GRO_filteredfile = arg
        elif opt in ("-a", "--afile"):
            SVs_annotation_vcffile = arg
        elif opt in ("-s", "--sfile"):
            SV_dis = arg
        elif opt in ("-o", "--odir"):
            outputdir = arg
    SV_dis=float(SV_dis)
    TSnearSV()
    t2 = time.time()
    print"The time needed in DNAlevel_NCL_filter for TSEs is:", t2 - t1

if __name__ == "__main__":
    main(sys.argv[1:])
