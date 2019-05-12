
#!/usr/bin/python2.7
# -*- coding: UTF-8 -*-

import sys, getopt,re

def red_depth(wigfile,step):
    fr_wig = [r1.strip().split("\t") for r1 in open(wigfile).readlines()]
    nw=len(fr_wig);wg_dict={}
    for pw in range(nw):
        par=fr_wig[pw][0]
        if wg_dict.has_key(par):
            wg_dict[par].append(fr_wig[pw])
        else:
            wg_dict[par]=[fr_wig[pw]]
    chr_list=wg_dict.keys()
    for pw1 in chr_list:
        tmp_list=wg_dict[pw1]
        n_r=len(tmp_list)
        len_s=int(float(tmp_list[n_r-1][2])/step)+1
        sim_list =[]
        index_r=0
        for index_s in range(len_s):
            if index_r<n_r:
                s1=step*index_s
                s2=step*(index_s+1)
                r1=float(tmp_list[index_r][1])
                r2=float(tmp_list[index_r][2])
                r_v=float(tmp_list[index_r][3])
                s_v=0
                pp = min(r2, s2) - max(r1, s1)
                if r1 == 0 and s2 <= r2:
                    if pp > 0 and s2 <r2:
                        r_v_t = pp * r_v
                        s_v = (s_v + r_v_t)/step
                        sim_list.append([str(s1), str(s2), str(s_v)])
                        continue
                    elif pp > 0 and s2 ==r2:
                        r_v_t = pp * r_v
                        s_v = (s_v + r_v_t) / step
                        sim_list.append([str(s1), str(s2), str(s_v)])
                        index_r = index_r + 1
                        continue
                else:
                    if pp>0:
                        if r2 <= s2:
                            s_v_a,index_r=calSimValue(tmp_list,index_r,s1,s2,r1,r2,r_v,n_r)
                            s_v=s_v_a/step
                            sim_list.append([str(s1), str(s2), str(s_v)])
                            continue
                        elif s2>r1 and s2<r2:
                            s_v= pp * r_v/step
                            sim_list.append([str(s1), str(s2), str(s_v)])
                            continue
        filename="la"
        if re.search('plus', wigfile,re.IGNORECASE):
            filename = pw1 + '_plus_strand_GROseq_' + str(step) + 'bp.txt'
        elif re.search('minus', wigfile,re.IGNORECASE):
            filename = pw1 + '_minus_strand_GROseq_' + str(step) + 'bp.txt'
        lf = open(filename, 'w+')
        tt_list=[]
        for it in sim_list:
            ss = '\t'.join(it)
            tt_list.append(ss)
        lf.writelines("\n".join(tt_list) + "\n")
        lf.close()


def calSimValue(tmp_list,index_r,s1,s2,r1,r2,r_v,n_r):
    s_v=0
    while r2 <= s2:
        pp = min(r2, s2) - max(r1, s1)
        s_v_t = pp * r_v
        s_v = s_v + s_v_t
        if index_r<n_r-1:
            index_r = index_r + 1
            r1 = float(tmp_list[index_r][1])
            r2 = float(tmp_list[index_r][2])
            r_v = float(tmp_list[index_r][3])
        else:
            index_r = index_r + 1
            break
    if s2>r1 and s2<r2:
        pp = min(r2, s2) - max(r1, s1)
        s_v_t = pp * r_v
        s_v = s_v + s_v_t
    return s_v,index_r

def main(argv):
    GRO_bedgraphfile = ''
    step_par= ''
    try:
        opts, args = getopt.getopt(argv, "hi:b:",["ifile=","bfile="])
    except getopt.GetoptError:
        print 'prepare_GROseq_data -i GRO_bedgraphfile -b step_par'
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print 'prepare_GROseq_data -i GRO_bedgraphfile -b step_par'
            sys.exit()
        elif opt in ("-i", "--ifile"):
            GRO_bedgraphfile = arg
        elif opt in ("-b", "--bfile"):
            step_par = arg
    if re.search('plus', GRO_bedgraphfile,re.IGNORECASE):
        print "This is plus strand"
        step=int(step_par)
        red_depth(GRO_bedgraphfile,step)
    if re.search('minus', GRO_bedgraphfile,re.IGNORECASE):
        print "This is minus strand"
        step = int(step_par)
        red_depth(GRO_bedgraphfile, step)
if __name__ == "__main__":
   main(sys.argv[1:])




