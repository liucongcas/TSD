import sys, getopt,time
import subprocess
def rm_inter_chrom(inputfile,filename):   #for fusion-files to delete the inter-chromosome fusions
    fr = [r.strip().split("\t") for r in open(inputfile).readlines()]
    list_filter=[]
    for pr in range(len(fr)):
        p1 = fr[pr][0].split("-")[0]
        p2 = fr[pr][0].split("-")[1]
        if p1 == p2:
            if float(fr[pr][4]) < min_sup and float(fr[pr][5]) == 0 and float(fr[pr][6]) == 0:
                pass
            else:
                tmp=[p1]+fr[pr][1:5]
                list_filter.append(tmp)
    tt_list=[]
    for iterm in list_filter:
        s = '\t'.join(iterm)
        tt_list.append(s)
    lw1 = open(filename, 'w+')
    lw1.writelines("\n".join(tt_list)  + "\n")
    lw1.close()
def prefilter(inputfile,genefile,gapfile,repfile1,repfile2,repswitch):
     # filter fusions located in inter-chromosomes
    inputfile2 = outputdir +"/"+inputfile_name+"_rminter.txt"
    rm_inter_chrom(inputfile, inputfile2)
    # filter out fusions smaller than 5k.
    subprocess.call("sed '/random/d' " + inputfile2+" | sed '/chrM/d' | awk '{if($3-$2>5000){print$0}}' > "+ inputfile2+"_del_random_filtered.txt",shell = True)
    # filter the fusions totally located in genes and get the finall list of fusions which have replicate support.
    fr_input = [r1.strip().split("\t") for r1 in open(inputfile2+"_del_random_filtered.txt").readlines()]
    fr_Gene = [r2.strip().split("\t") for r2 in open(genefile).readlines()]
    fr_gapfile=[r2.strip().split("\t") for r2 in open(gapfile).readlines()]
    n1 = len(fr_input);n2=len(fr_gapfile);filtered_list=[]
    input_list=[];n_gene = len(fr_Gene);finall_list = [];intraGene_list = []
    for pp1 in range(n1):  # get the fusions totally located in genes
        for pp2 in range(n_gene):
            if fr_input[pp1][0] == fr_Gene[pp2][0]:
                if float(fr_input[pp1][1]) >= float(fr_Gene[pp2][3]) and float(fr_input[pp1][2]) <= float(fr_Gene[pp2][4]):
                    intraGene_list.append(fr_input[pp1])
                    break
    for p in range(n1): #filter the fusions totally located in genes.
        if fr_input[p] not in intraGene_list and float(fr_input[p][1])>0:
            input_list.append(fr_input[p])
    n3=len(input_list)
    for q in range(n3):
        junc_l_s=float(input_list[q][1])-2000
        junc_l_e=float(input_list[q][1])+2000
        junc_r_s=float(input_list[q][2])-2000
        junc_r_e=float(input_list[q][2])+2000
        junc_chr=input_list[q][0]
        for q1 in range(n2):
            gap_chr=fr_gapfile[q1][0]
            gap_s=float(fr_gapfile[q1][1])
            gap_e=float(fr_gapfile[q1][2])
            if junc_chr==gap_chr:
                if (min(gap_e,junc_l_e)-max(gap_s,junc_l_s))>=0 or (min(gap_e,junc_r_e)-max(gap_s,junc_r_s))>=0:
                    filtered_list.append(input_list[q])
                    break
    for q2 in input_list:
        if q2 not in filtered_list:
            finall_list.append(q2)
    if repswitch =="Y":
        finall_list3 = get_rep_support(repfile1, repfile2, finall_list) #get the finall list of juctions which have replicate support
    else:
        finall_list3=finall_list
    write_name=outputdir+"/"+inputfile_name
    t2_list=[]
    for iterm in finall_list3:
        s = '\t'.join(iterm)
        t2_list.append(s)
    if len(t2_list)>0:
        lw = open(write_name + "_prefiltered.txt", 'w+')
        lw.writelines("\n".join(t2_list)+"\n")
        lw.close()
def get_rep_support(repfile1,repfile2,input_list): # get replicate support for the filtered juction
    repfile1_name = [i for i in repfile1.split("/") if i != ""][-1]
    repfile2_name = [i for i in repfile2.split("/") if i != ""][-1]
    repfile1_name2 = outputdir + "/" + repfile1_name + "_rminter.txt"
    repfile2_name2 = outputdir + "/" + repfile2_name + "_rminter.txt"
    rm_inter_chrom(repfile1,repfile1_name2)
    rm_inter_chrom(repfile2,repfile2_name2)
    rep1 = [t1.strip().split("\t") for t1 in open(repfile1_name2).readlines()]
    rep2 = [t2.strip().split("\t") for t2 in open(repfile2_name2).readlines()]
    m1 = len(rep1);m2 = len(rep2);m3 = len(input_list)
    rep1_dict = {};comm_dict = {};listfilter = []
    for i in range(m1):
        par1 = rep1[i][0]
        if rep1_dict.has_key(par1):
            rep1_dict[par1].append([float(rep1[i][1]),float(rep1[i][2])])
        else:
            rep1_dict[par1] =[[float(rep1[i][1]), float(rep1[i][2])]]
    for j in range(m2):
        t1_chr = rep2[j][0]
        t1_l = float(rep2[j][1])
        t1_r = float(rep2[j][2])
        if rep1_dict.has_key(t1_chr):
            list_tmp=rep1_dict[t1_chr]
            for j2 in list_tmp:
                s1_l = j2[0]
                s1_r = j2[1]
                if abs(s1_l - t1_l) <= 6 and abs(s1_r - t1_r) <= 6:
                    if comm_dict.has_key(t1_chr):
                        comm_dict[t1_chr].append([t1_l,t1_r])
                        break
                    else:
                        comm_dict[t1_chr]=[[t1_l, t1_r]]
    for q in range(m3):
        t_chr = input_list[q][0]
        t_l = float(input_list[q][1])
        t_r = float(input_list[q][2])
        if comm_dict.has_key(t_chr):
            list_temp2=comm_dict[t_chr]
            for p in list_temp2:
                s_l =p[0]
                s_r =p[1]
                if abs(s_l-t_l)<=6 and abs(s_r-t_r)<=6:
                    listfilter.append(input_list[q])
                    break
    return listfilter
def main(argv):
    global outputdir
    global min_sup
    inputfile = ''
    genefile= ''
    gapfile= ''
    repfile1 = ''
    repfile2 = ''
    repswitch = ''
    min_sup=""
    outputdir = ''
    try:
        opts, args = getopt.getopt(argv,"hi:g:w:r:t:p:m:o:",["ifile=","gfile=","wfile","r1file=","r2file=","pfile=","mfile=","ofile="])
    except getopt.GetoptError:
        print 'pre_filter_fusions.py -i <inputfile> -g <genefile> -w <repswitch>  -r <repfile1> -t <repfile2> -p <gapfile> -m <min_sup> -o <outputdir>'
        sys.exit(8)
    for opt, arg in opts:
        if opt == '-h':
            print 'pre_filter_fusions.py -i <inputfile> -g <genefile> -w <repswitch> -r <repfile1> -t <repfile2> -p <gapfile> -m <min_sup> -o <outputdir>'
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-g", "--gfile"):
            genefile = arg
        elif opt in ("-w","--wfile"):
            repswitch=arg
        elif opt in ("-r", "--r1file"):
            repfile1 = arg
        elif opt in ("-t", "--r2file"):
            repfile2 = arg
        elif opt in ("-p", "--pfile"):
            gapfile = arg
        elif opt in ("-m","--mfile"):
            min_sup=arg
        elif opt in ("-o", "--odir"):
            outputdir = arg
    t1=time.time()
    min_sup=float(min_sup)
    global inputfile_name
    inputfile_name = [i for i in inputfile.split("/") if i != ""][-1]
    prefilter(inputfile, genefile, gapfile,repfile1, repfile2,repswitch)
    subprocess.call("rm "+outputdir +"/"+"*_rminter.txt*",shell=True)
    t2=time.time()
    print "The time needed in pre_filter for fusions is:",t2-t1

if __name__ == "__main__":
   main(sys.argv[1:])

