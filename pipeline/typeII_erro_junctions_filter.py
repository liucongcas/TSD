
import time,sys,getopt,os

def bl8_jun_filter(bl8_format,junfile,genefile,criterion_len,prefix,gap_j):
    fr_jun_bl8 = [r1.strip().split("\t") for r1 in open(bl8_format).readlines()]
    fr_jun_type1 = [r2.strip().split("\t") for r2 in open(junfile).readlines()]
    fr_genebody=[r3.strip().split("\t") for r3 in open(genefile).readlines()]
    nf=len(fr_jun_bl8);type2_erro_list=[];jun_dict={}
    ns= len(fr_jun_type1);seg_dict={}
    ne=len(fr_genebody);gene_plus_dict={};gene_minus_dict={}
    for b_j in range(nf):
        par_j = fr_jun_bl8[b_j][0].split("*")[0]
        if jun_dict.has_key(par_j):
            jun_dict[par_j].append(fr_jun_bl8[b_j])
        else:
            jun_dict[par_j]=[fr_jun_bl8[b_j]]
    for b_e in range(ne):
        par_e = [fr_genebody[b_e][0], float(fr_genebody[b_e][3]), float(fr_genebody[b_e][4]), fr_genebody[b_e][7],fr_genebody[b_e][8]]
        par_chr=fr_genebody[b_e][0]
        if fr_genebody[b_e][6]=="+":
            if gene_plus_dict.has_key(par_chr):
                gene_plus_dict[par_chr].append(par_e)
            else:
                gene_plus_dict[par_chr]=[par_e]
        elif fr_genebody[b_e][6]=="-":
            if gene_minus_dict.has_key(par_chr):
                gene_minus_dict[par_chr].append(par_e)
            else:
                gene_minus_dict[par_chr]=[par_e]
    for b_s in range(ns):
        par_s=fr_jun_type1[b_s][0]+"_"+fr_jun_type1[b_s][1]+"_"+fr_jun_type1[b_s][2]
        seg_dict[par_s]=[fr_jun_type1[b_s][5],fr_jun_type1[b_s][3]]
    for vx in jun_dict.keys():
        jun_list=jun_dict[vx]
        jun_l_l=float(jun_dict[vx][0][0].split("*")[1].split("~")[0])
        jun_l_r=float(jun_dict[vx][0][0].split("*")[1].split("~")[1])
        jun_l_mid = (jun_l_l+jun_l_r) / 2
        jun_r_l=float(jun_dict[vx][0][0].split("*")[1].split("~")[2])
        jun_r_r=float(jun_dict[vx][0][0].split("*")[1].split("~")[3][:-2])
        jun_r_mid = (jun_r_l +jun_r_r) / 2
        chr_w=seg_dict[vx][1]
        nchr = jun_dict[vx][0][0].split("_")[0]
        leftfactor = False
        rightfactor = False
        leftfactor, rightfactor, fus_list_els = find_primary_map(jun_list,jun_l_mid, jun_r_mid,nchr,gap_j,leftfactor,rightfactor)
        if leftfactor and rightfactor:
            for pj in jun_list:
                old_len=len(type2_erro_list)
                type2_erro_list=filter_one_jun(pj,nchr,jun_l_mid, jun_r_mid, jun_l_l, jun_l_r, jun_r_l, jun_r_r, type2_erro_list, chr_w, gene_plus_dict, gene_minus_dict, vx,criterion_len,gap_j)
                if len(type2_erro_list)> old_len:
                    break
        else:
            type2_erro_list.append([vx] + jun_list[1]+["couldnot mapping back"])
    lf = open(prefix+"_wrong_bl8",'w+')
    t1_list=[]
    for it in type2_erro_list:
        ss = '\t'.join(it)
        t1_list.append(ss)
    lf.writelines("\n".join(t1_list) + "\n")
    lf.close()
    return type2_erro_list
def find_primary_map(jun_list,jun_l_mid, jun_r_mid,nchr,gap_j,leftfactor,rightfactor):
    jun_list_pre=[];gap_j=gap_j+20
    for pj_1 in jun_list:
        q_l = float(pj_1[6])
        q_r = float(pj_1[7])
        len_seq = float(pj_1[3])
        map_identity = float(pj_1[2])
        t_l_1 = min(float(pj_1[8]), float(pj_1[9]))
        t_r_1 = max(float(pj_1[8]), float(pj_1[9]))
        mid_t_1=(t_l_1+t_r_1)/2
        if map_identity > 90 and abs(len_seq -100)<=30 and pj_1[1]==nchr:
            if abs(mid_t_1-jun_l_mid)<=gap_j or abs(mid_t_1-jun_r_mid)<=gap_j:
                if q_r >= 100 - gap_j and q_r <= 100 + gap_j:
                    jun_list_pre.append(pj_1)
                    leftfactor=True
                elif q_l >= 100-gap_j and q_l <= 100+gap_j:
                    jun_list_pre.append(pj_1)
                    rightfactor = True
    jun_list_els= [kk for kk in jun_list if kk not in jun_list_pre]
    return leftfactor,rightfactor,jun_list_els

def filter_one_jun(pj,nchr,jun_l_mid,jun_r_mid,jun_l_l,jun_l_r,jun_r_l,jun_r_r,type2_erro_list,chr_w,gene_plus_dict,gene_minus_dict,vx,criterion_len,gap_j):
    q_l = float(pj[6]);q_r = float(pj[7]);t_l = min(float(pj[8]),float(pj[9]));t_r = max(float(pj[8]),float(pj[9]))
    mid_t = (t_l + t_r) / 2;len_seq = float(pj[3]);map_identity = float(pj[2])
    chr_p = pj[1]
    if abs(mid_t - jun_l_mid) > 50 and abs(mid_t - jun_r_mid) > 50 and map_identity > 85:
        if (q_r >= 100-gap_j and q_r <= 100+gap_j) or (q_l >= 100-gap_j and q_l <= 100+gap_j):
            if map_identity ==100 and len_seq > criterion_len:
                type2_erro_list.append([vx] + pj+['SameIdentityfalse_1'])
            elif map_identity >= 95 and len_seq >= 100-gap_j-20 and len_seq<=100+gap_j+20:
                type2_erro_list.append([vx] + pj + ['SameIdentityfalse_2'])
            elif len_seq >20 and gene_plus_dict.has_key(chr_p) and gene_minus_dict.has_key(chr_p):
                type2_erro_list = filter_same_gene(vx, pj,nchr,t_l, t_r, gene_plus_dict, gene_minus_dict,type2_erro_list,chr_w,jun_l_l,jun_l_r,jun_r_l, jun_r_r)
        elif q_l <= 60 and q_r >= 140 and len_seq > criterion_len and map_identity > 95:
            type2_erro_list.append([vx] + pj+["Co_linear explanation"])
    return type2_erro_list

def filter_same_gene(vx, pj,nchr,t_l, t_r,gene_plus_dict, gene_minus_dict, type2_erro_list,chr_w,jun_l_l,jun_l_r,jun_r_l, jun_r_r):
    gene_jun_l=gene_jun_r = 0;gene_q = 1;chr_t=pj[1];gene_list=[]
    if chr_w == "+":
        gene_list = gene_plus_dict[chr_t]
    elif chr_w == "-":
        gene_list = gene_minus_dict[chr_t]
    for p_e1 in gene_list:
        if nchr == p_e1[0] and (min(jun_l_r,p_e1[2])-max(jun_l_l,p_e1[1]))>20:
            gene_jun_l = p_e1[3:]
            break
        else:
            gene_jun_l = ['none1','none1']
    for p_e2 in gene_list:
        if nchr == p_e2[0] and (min(jun_r_r,p_e2[2])-max(jun_r_l,p_e2[1]))>20:
            gene_jun_r = p_e2[3:]
            break
        else:
            gene_jun_r = ['none2', 'none2']
    for p_e3 in gene_list:
        if chr_t == p_e3[0] and (min(t_r,p_e3[2])-max(t_l,p_e3[1]))>20:
            gene_q = p_e3[3:]
            break
        else:
            gene_q=['none3', 'none3']
    if gene_q[0]==gene_jun_l[0] or gene_q[0]==gene_jun_r[0] or gene_jun_l[0]==gene_jun_r[0]:
        type2_erro_list.append([vx] + pj+['SameGeneFalse'])
    else:
        if gene_q[1] != 'none':
            if gene_q[1] == gene_jun_l[1] or gene_q[1] == gene_jun_r[1] or gene_jun_l[1]==gene_jun_r[1]:
                if gene_jun_r[1] != "none":
                    type2_erro_list.append([vx] + pj + ['SameGeneFalse'])
    return type2_erro_list

def co_jun_filter(bl8file,typeI_junfile,prefix):
    fr_bl8 = [r1.strip().split("\t") for r1 in open(bl8file).readlines()]
    fr_jun = [r1.strip().split("\t") for r1 in open(typeI_junfile).readlines()]
    n_jun=len(fr_jun);all_jun={}
    list_jun=[];filtered_jun=[]
    for ap in range(n_jun):
        pp=fr_jun[ap][0]+'_'+fr_jun[ap][1]+"_"+fr_jun[ap][2]
        all_jun[pp]=1
    for bp in fr_bl8:
        pp1=bp[0]
        if all_jun.has_key(pp1):
            all_jun[pp1]=all_jun[pp1]+1
    for (k,v) in all_jun.items():
        if v==1:
            list_jun.append(k)
    for ap2 in range(n_jun):
        pp3=fr_jun[ap2][0]+'_'+fr_jun[ap2][1]+"_"+fr_jun[ap2][2]
        if pp3 in list_jun:
            filtered_jun.append(fr_jun[ap2])
    t2_list=[]
    for it2 in filtered_jun:
        ss2 = '\t'.join(it2[:5]+it2[-4:])
        t2_list.append(ss2)
    if len(t2_list)>0:
        lf2 = open(prefix + "_typeIIfiltered.txt", 'w+')
        lf2.writelines("\n".join(t2_list) + "\n")
        lf2.close()
    return all_jun,list_jun

def main(argv):
    t1 = time.time()
    typeI_junfile =''
    bl8_format=''
    genefile=''
    try:
        opts, args = getopt.getopt(argv, "hi:b:e:", ["ifile=","bfile=","efile="])
    except getopt.GetoptError:
        print 'python typeII_erro_junctions_filter.py -i <typeI erro filtered junctions file> -b <bl8_format file> -e genefile '
        sys.exit(3)
    for opt, arg in opts:
        if opt == '-h':
            print 'python typeII_erro_junctions_filter.py -i <typeI erro filtered junctions file> -b <bl8_format file> -e genefile'
            sys.exit()
        elif opt in ("-i", "--ifile"):
            typeI_junfile = arg
        elif opt in ("-b", "--bfile"):
            bl8_format = arg
        elif opt in ("-e", "--efile"):
            genefile = arg
    prefix=typeI_junfile[:-4]
    gap_j=4
    bl8_jun_filter(bl8_format, typeI_junfile, genefile, 80,prefix,gap_j)
    co_jun_filter(prefix+"_wrong_bl8",typeI_junfile,prefix)
    os.system("rm " + prefix + "*exon*")
    os.system("rm " + prefix + "*bl8*")
    t2=time.time()
    print "The time needed in typeII error filter for junctions is:",t2-t1

if __name__ == "__main__":
   main(sys.argv[1:])


