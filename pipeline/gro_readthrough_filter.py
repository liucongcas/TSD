#!/user/bin/python

import subprocess,os,re
import time,sys,getopt
def readthroughfiltered(juncfile, step,window,readthroughfiledir,total_len):
    fr_junc = [r1.strip().split("\t") for r1 in open(juncfile).readlines()]
    nj=len(fr_junc);jun_dict={};filterd_list=[]
    for pj in range(1,nj):
        par=fr_junc[pj][0]
        if jun_dict.has_key(par):
            jun_dict[par].append(fr_junc[pj])
        else:
            jun_dict[par]=[fr_junc[pj]]
    chr_list=jun_dict.keys()
    readthrough_files=os.listdir(readthroughfiledir)
    for p1 in chr_list:
        p2=p1+"_"
        file_plus=[i for i in readthrough_files if re.search(p2,i) and re.search('plus',i)][0]
        file_minus=[i for i in readthrough_files if re.search(p2,i) and re.search('minus',i)][0]
        fr_plus=[r1.strip().split("\t") for r1 in open(readthroughfiledir+"/"+file_plus).readlines()]
        fr_minus=[r1.strip().split("\t") for r1 in open(readthroughfiledir+"/"+file_minus).readlines()]
        tmp_list = jun_dict[p1]
        for pr in tmp_list:
            junc_l=min(int(float(pr[1])),int(float(pr[2])))
            junc_r=max(int(float(pr[1])),int(float(pr[2])))
            read_l=int(junc_l/step)
            read_r=int(junc_r/step)
            windowCover_list_plus,stepCover_list_plus = get_cover_list(fr_plus,read_l,read_r,step,window)
            windowCover_list_minus,stepCover_list_minus = get_cover_list(fr_minus,read_l,read_r,step,window)
            zero_list_plus = get_window_zero_list(windowCover_list_plus,junc_l,junc_r,window,total_len)
            zero_list_minus =get_window_zero_list(windowCover_list_minus,junc_l,junc_r,window,total_len)
            pr_plus=len(zero_list_plus)
            pr_minus=len(zero_list_minus)
            if pr_plus != 0 and pr_minus != 0:
                filterd_list.append(pr)
            elif pr_plus != 0 and pr_minus == 0:
                zero_list_minus = get_step_zero_list(stepCover_list_minus, step,total_len)
                if  len(zero_list_minus) != 0:
                    filterd_list.append(pr)
            elif pr_plus == 0 and pr_minus != 0:
                zero_list_plus = get_step_zero_list(stepCover_list_plus, step,total_len)
                if len(zero_list_plus) != 0:
                    filterd_list.append(pr )
            elif pr_plus != 0 and pr_minus != 0:
                zero_list_plus=get_step_zero_list(stepCover_list_plus,step,total_len)
                zero_list_minus=get_step_zero_list(stepCover_list_minus,step,total_len)
                if len(zero_list_plus) != 0 and len(zero_list_minus) != 0:
                    filterd_list.append(pr)
    tt_list=[]
    tt_list.append("\t".join(["Chromosome", "5_JunctionSites", "3_JunctionSites", "SplicingType", "Surpported JunctionsReads","5_JunctionSites_strand",
                              "3_JunctionSites_strand", "5_Junction_coordinates", "3_Junction_coordinates", "5_parent_gene", "3_parent_gene","SV potential"]))
    for it in filterd_list:
        ss = '\t'.join(map(str, it))
        tt_list.append(ss)
    if len(tt_list)>0:
        lf = open(juncfile + "_GRO_RTfiltered.txt", 'w+')
        lf.writelines("\n".join(tt_list) + "\n")
        lf.close()


def get_step_zero_list(stepCover_list,step,total_len): ##find zero transcription fragement in step size
    n_t = len(stepCover_list);nonStep_list=[];zero_list=[]
    for qt in range(n_t):
        if stepCover_list[qt] != 0:
            nonStep_list.append(qt)
    if nonStep_list[0] != 0:
        nonStep_list.insert(0, -1)
    if nonStep_list[-1]<=n_t-1:
        nonStep_list.append(n_t-1)
    zero_list_tmp = [nonStep_list[i + 1] - nonStep_list[i] - 1 for i in range(len(nonStep_list) - 1)]
    for q4 in zero_list_tmp:
        if q4 >= total_len/float(step):
            zero_list.append(q4 * step)
    return zero_list

def get_cover_list(fr_file,read_l,read_r,step,window): ##get the read coverage of everry window from jun_l to jun_r
    windowCover_list = [];q_step=window/step;stepCover_list=[]
    for q in range(read_l,read_r+1,q_step):
        if q <read_r+1-q_step:
            binC = 0
            for q in range(q,q+q_step):
                binC+=float(fr_file[q][2])
                stepCover_list.append(float(fr_file[q][2]))
            windowCover_list.append(binC)
    return windowCover_list,stepCover_list

def get_window_zero_list(windowCover_list,junc_l,junc_r,window,total_len):
    zero_list=[];n_b=len(windowCover_list)
    if sum(windowCover_list)==0:  ##if all of the windows is zero
        zero_list.append(junc_r-junc_l)
    else:
        non_list=[]
        index_x=False
        for q2 in range(n_b):
            if windowCover_list[q2] != 0:
                non_list.append(q2)
                if windowCover_list[q2]<1:
                    index_x=True
        if len(non_list)<n_b-1: ##filtered all the bin has reads
            if non_list[0]!=0:
                non_list.insert(0,-1)
            if  non_list[-1]<n_b-1:
                non_list.append(n_b-1)
            nn=len(non_list)
            if index_x:
                zero_list_tmp1=[]
                zero_list_tmp2=[]
                for q3 in range(nn-1):
                    inx=non_list[q3]
                    if windowCover_list[inx] <= 1:
                        ex_l=non_list[q3-1]
                        ex_r=non_list[q3+1]
                        if inx-ex_l-1>=2 and ex_r-inx-1>=2:  ##4kbin is zero
                            zero_list_tmp1.append([ex_l,ex_r])
                zero_list_tmp = [[non_list[i],non_list[i + 1]] for i in range(nn - 1)]
                for pp1 in zero_list_tmp:
                    if  pp1[1]-pp1[0]-1 >= total_len/float(window):
                        zero_list_tmp2.append(pp1)
                zero_list_tmp_total=sorted(zero_list_tmp1+zero_list_tmp2, key=lambda gene: gene[1])
                zero_list_tmp_total=mergeList(zero_list_tmp_total)
                if len(zero_list_tmp_total)!=0:
                    zero_list=[(i[1]-i[0]-1)*window for i in zero_list_tmp_total]
            else:
                zero_list_tmp_ad= [non_list[i + 1] - non_list[i] - 1 for i in range(nn - 1)]
                for q4 in zero_list_tmp_ad:
                    if q4>=total_len/float(window):
                        zero_list.append(q4*window)
    return  zero_list

def mergeList(list_for_merge):
    merge_x = 0
    merge_y = 0
    merge_list = []
    for merge_i in range(len(list_for_merge)):
        p_merge = min(list_for_merge[merge_i][1], merge_y) - max(list_for_merge[merge_i][0], merge_x)
        if merge_i < len(list_for_merge) - 1:
            if p_merge >= 0:
                merge_x = min(list_for_merge[merge_i][0], merge_x)
                merge_y = max(list_for_merge[merge_i][1], merge_y)
            else:
                merge_list.append([merge_x, merge_y])
                merge_x = list_for_merge[merge_i][0]
                merge_y = list_for_merge[merge_i][1]
        else:
            if p_merge >= 0:
                merge_x = min(list_for_merge[merge_i][0], merge_x)
                merge_y = max(list_for_merge[merge_i][1], merge_y)
                merge_list.append([merge_x, merge_y])
            else:
                merge_list.append([merge_x, merge_y])
                merge_list.append(list_for_merge[merge_i])
    if [0,0] in merge_list:
        merge_list.remove([0,0])
    return merge_list


def main(argv):
    t1 = time.time()
    juncfile = ''
    readthroughfiledir = ''
    total_len=''
    try:
        opts, args = getopt.getopt(argv, "hi:d:t:",["ifile=","ddir=","tfile="])
    except getopt.GetoptError:
        print 'python gro_readthrough_readthrouth_filter -i <typeII erro filtered file> -d <readthroughfiledir> -t <total_len>'
        sys.exit(3)
    for opt, arg in opts:
        if opt == '-h':
            print 'python gro_readthrough_readthrouth_filter -i <typeII erro filtered file> -d <readthroughfiledir> -t <total_len>'
            sys.exit()
        elif opt in ("-i", "--ifile"):
            juncfile = arg
        elif opt in ("-d", "--ddir"):
            readthroughfiledir = arg
        elif opt in ("-t", "--tfile"):
            total_len = arg
    step=200
    window=2000
    total_len=int(total_len)
    readthroughfiltered(juncfile, step,window,readthroughfiledir,total_len)
    t2=time.time()
    print "The time needed in gro_readthrough_readthrough_filter for fusions is:",t2-t1
if __name__ == "__main__":
   main(sys.argv[1:])









