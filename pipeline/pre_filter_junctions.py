import sys, getopt,time
import subprocess

def  prefilter(inputfile2, genefile, gapfile,repfile1,repfile2,repswitch):
    # filter juctions smaller than 5k.
    subprocess.call("sed '/random/d' " + outputdir +"/"+inputfile2+" |  sed '/chrM/d' | awk '{if($3-$2>5000){print$0}}' > "+ outputdir +"/"+ inputfile2+"_del_random_filtered.txt",shell = True)
    fr_input = [r1.strip().split("\t") for r1 in open(outputdir +"/"+inputfile2+"_del_random_filtered.txt").readlines()]
    fr_Gene = [r2.strip().split("\t") for r2 in open(genefile).readlines()]
    fr_gapfile = [r2.strip().split("\t") for r2 in open(gapfile).readlines()]
    n1 = len(fr_input);n2 = len(fr_gapfile);filtered_list = [];input_list = []
    n_gene=len(fr_Gene);finall_list=[]
    intraGene_list=[]
    for pp1 in range(n1): # get the juctions totally located in genes
        for pp2 in range(n_gene):
            if fr_input[pp1][0]==fr_Gene[pp2][0]:
                if float(fr_input[pp1][1]) >= float(fr_Gene[pp2][3]) and float(fr_input[pp1][2])<=float(fr_Gene[pp2][4]):
                    intraGene_list.append(fr_input[pp1])
                    break
    for p in range(n1): #filter the juctions totally located in genes.
        if fr_input[p] not in intraGene_list and float(fr_input[p][1])>0:
            input_list.append(fr_input[p])
    n3 = len(input_list)
    for q in range(n3): # filter the junctions located near a gap
        junc_l_s = float(input_list[q][1]) - 2000
        junc_l_e = float(input_list[q][1]) + 2000
        junc_r_s = float(input_list[q][2]) - 2000
        junc_r_e = float(input_list[q][2]) + 2000
        junc_chr = input_list[q][0]
        for q1 in range(n2):
            gap_chr = fr_gapfile[q1][0]
            gap_s = float(fr_gapfile[q1][1])
            gap_e = float(fr_gapfile[q1][2])
            if junc_chr == gap_chr:
                if (min(gap_e, junc_l_e) - max(gap_s, junc_l_s)) >= 0 or (min(gap_e, junc_r_e) - max(gap_s, junc_r_s)) >= 0:
                    filtered_list.append(input_list[q])
                    break
    for q2 in input_list:
        if q2 not in filtered_list:
            finall_list.append(q2)
    if repswitch =="Y":
        finall_list2 = get_rep_support(repfile1, repfile2, finall_list) #get the finall list of juctions which have replicate support
    else:
        finall_list2=finall_list
    tt_list=[]
    for iterm in finall_list2:
        s = '\t'.join(iterm)
        tt_list.append(s)
    if len(tt_list) >0:
        lw = open(outputdir + "/" + inputfile_name + "_prefiltered.txt", 'w+')
        lw.writelines("\n".join(tt_list)+"\n")
        lw.close()
def get_rep_support(repfile1,repfile2,input_list): #get replicate support for the filtered juction
    repfile1_name = [i for i in repfile1.split("/") if i != ""][-1]
    repfile2_name = [i for i in repfile2.split("/") if i != ""][-1]
    subprocess.call(softwaredir+"/tools/bed_to_juncs < " + repfile1 + " > " + outputdir + "/" + repfile1_name + "_changed", shell=True)
    subprocess.call(softwaredir+"/tools/bed_to_juncs < " + repfile2 + " > " + outputdir + "/" + repfile2_name + "_changed", shell=True)
    rep1 = [t1.strip().split("\t") for t1 in open(outputdir + "/" + repfile1_name+"_changed").readlines()]
    rep2 = [t2.strip().split("\t") for t2 in open(outputdir + "/" + repfile2_name+"_changed").readlines()]
    m1 = len(rep1);m2 = len(rep2);m3 = len(input_list)
    rep1_dict = {};comm_dict = {};listfilter = []
    for i in range(m1):
        par1 = rep1[i][0]
        if rep1_dict.has_key(par1):
            a=rep1_dict[par1]
            a.append([float(rep1[i][1]),float(rep1[i][2])])
            rep1_dict[par1]=a
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
                if abs(s1_l - t1_l) <= 3 and abs(s1_r - t1_r) <= 3:
                    if comm_dict.has_key(t1_chr):
                        b=comm_dict[t1_chr]
                        b.append([t1_l,t1_r])
                        comm_dict[t1_chr]=b
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
    global inputfile
    global outputdir
    global softwaredir
    inputfile = ''
    genefile= ''
    repswitch=''
    repfile1=''
    repfile2=''
    gapfile = ''
    outputdir=''
    softwaredir = ""
    try:
        opts, args = getopt.getopt(argv,"hi:g:w:r:t:p:o:s:",["ifile=","gfile=","wfile","r1file=","r2file=","pfile=","odir=","sfile="])
    except getopt.GetoptError:
        print 'pre_filter_junctions.py -i <inputfile> -g <genefile> -w <repswitch> -r <repfile1> -t <repfile2> -p <gapfile> -o <outputdir> -s <softwaredir>'
        sys.exit(7)
    for opt, arg in opts:
        if opt == '-h':
            print 'pre_filter_junctions.py -i <inputfile> -g <genefile> -w <repswitch> -r <repfile1> -t <repfile2> -p <gapfile> -o <outputdir> -s <softwaredir> '
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
        elif opt in ("-o", "--odir"):
            outputdir = arg
        elif opt in ("-s", "--sfile"):
            softwaredir = arg
    t1=time.time()
    global inputfile_name
    inputfile_name=[i for i in inputfile.split("/") if i !=""][-1]
    subprocess.call(softwaredir+"/tools/bed_to_juncs < " + inputfile + " > " + outputdir +"/"+ inputfile_name + "_changed",shell = True)
    subprocess.call("sed '1d' " + inputfile + " | awk '{print $11"'"\t"'"$12"'"\t"'"$4}' > "+ outputdir +"/"+"t1 ",shell = True)
    subprocess.call("paste " + outputdir +"/"+ inputfile_name + "_changed "+ outputdir +"/"+ "t1 > " + outputdir +"/"+ inputfile_name + "_changed.txt",shell = True)
    subprocess.call("rm "+ outputdir +"/"+"t1",shell = True)
    prefilter(inputfile_name+"_changed.txt", genefile, gapfile, repfile1, repfile2,repswitch)
    subprocess.call("rm "+ outputdir +"/"+"*_changed", shell=True)
    subprocess.call("rm " + outputdir +"/"+ inputfile_name + "_changed.txt*", shell=True)
    t2 = time.time()
    print "The time needed in prefiltering for junctions is:", t2 - t1


if __name__ == "__main__":
   main(sys.argv[1:])

