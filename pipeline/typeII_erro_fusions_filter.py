import time,sys,getopt,os

def bl8_fus_filter(bl8_fusfile,typeI_fusfile,genefile,criterion_len,prefix,gap_j):
    fr_fus_bl8 = [r1.strip().split("\t") for r1 in open(bl8_fusfile).readlines()]
    fr_fus_type1 = [r2.strip().split("\t") for r2 in open(typeI_fusfile).readlines()]
    fr_gene=[r3.strip().split("\t") for r3 in open(genefile).readlines()]
    nf=len(fr_fus_bl8);type2_erro_list=[];fus_dict={}
    ns= len(fr_fus_type1);seg_dict={}
    ne=len(fr_gene);gene_plus_dict={};gene_minus_dict={}
    for b_j in range(nf):
        par_j = fr_fus_bl8[b_j][0].split("*")[0]
        if fus_dict.has_key(par_j):
            fus_dict[par_j].append(fr_fus_bl8[b_j])
        else:
            fus_dict[par_j]=[fr_fus_bl8[b_j]]
    for b_e in range(ne):
        par_e = [fr_gene[b_e][0], float(fr_gene[b_e][3]), float(fr_gene[b_e][4]), fr_gene[b_e][7],fr_gene[b_e][8]]
        par_chr=fr_gene[b_e][0]
        if fr_gene[b_e][6]=="+":
            if gene_plus_dict.has_key(par_chr):
                gene_plus_dict[par_chr].append(par_e)
            else:
                gene_plus_dict[par_chr]=[par_e]
        elif fr_gene[b_e][6]=="-":
            if gene_minus_dict.has_key(par_chr):
                gene_minus_dict[par_chr].append(par_e)
            else:
                gene_minus_dict[par_chr]=[par_e]
    for b_s in range(ns):
        par_s=fr_fus_type1[b_s][0]+"_"+fr_fus_type1[b_s][1]+"_"+fr_fus_type1[b_s][2]
        seg_dict[par_s]=[fr_fus_type1[b_s][8],fr_fus_type1[b_s][9]]
    for vx in fus_dict.keys():
        fus_list=fus_dict[vx]
        fus_l_l=float(fus_dict[vx][0][0].split("*")[1].split("~")[0])
        fus_l_r=float(fus_dict[vx][0][0].split("*")[1].split("~")[1])
        fus_l_mid = (fus_l_l+fus_l_r) / 2
        fus_r_l=float(fus_dict[vx][0][0].split("*")[1].split("~")[2])
        fus_r_r=float(fus_dict[vx][0][0].split("*")[1].split("~")[3][:-3])
        fus_r_mid = (fus_r_l +fus_r_r) / 2
        nchr = fus_dict[vx][0][0].split("_")[0]
        leftfactor=False
        rightfactor=False
        leftfactor, rightfactor, fus_list_els=find_primary_map(fus_list,fus_l_mid, fus_r_mid,nchr,gap_j,leftfactor,rightfactor)
        if leftfactor and rightfactor:
            for pj in fus_list_els:
                old_len=len(type2_erro_list)
                type2_erro_list=filter_one_fus(pj,nchr,fus_l_mid, fus_r_mid, fus_l_l, fus_l_r, fus_r_l, fus_r_r, type2_erro_list,gene_plus_dict, gene_minus_dict, vx,criterion_len,gap_j)
                if len(type2_erro_list)> old_len:
                    break
        else:
            type2_erro_list.append([vx] + fus_list[0]+["couldnot mapping back"])
    lf = open(prefix+"_wrong_bl8",'w+')
    t1_list=[]
    for it in type2_erro_list:
        ss = '\t'.join(it)
        t1_list.append(ss)
    lf.writelines("\n".join(t1_list) + "\n")
    lf.close()
    return type2_erro_list

def find_primary_map(fus_list,fus_l_mid, fus_r_mid,nchr,gap_j,leftfactor,rightfactor):
    fus_list_pre=[];gap_j=gap_j+20
    for pj_1 in fus_list:
        q_l = float(pj_1[6])
        q_r = float(pj_1[7])
        len_seq = float(pj_1[3])
        map_identity = float(pj_1[2])
        t_l_1 = min(float(pj_1[8]), float(pj_1[9]))
        t_r_1 = max(float(pj_1[8]), float(pj_1[9]))
        mid_t_1=(t_l_1+t_r_1)/2
        if map_identity > 90 and abs(len_seq -100)<=30 and pj_1[1]==nchr:
            if abs(mid_t_1-fus_l_mid)<=gap_j or abs(mid_t_1-fus_r_mid)<=gap_j:
                if q_r >= 100 - gap_j and q_r <= 100 + gap_j:
                    fus_list_pre.append(pj_1)
                    leftfactor=True
                elif q_l >= 100-gap_j and q_l <= 100+gap_j:
                    fus_list_pre.append(pj_1)
                    rightfactor = True
    fus_list_els= [kk for kk in fus_list if kk not in fus_list_pre]
    return leftfactor,rightfactor,fus_list_els


def filter_one_fus(pj,nchr,fus_l_mid,fus_r_mid,fus_l_l,fus_l_r,fus_r_l,fus_r_r,type2_erro_list,gene_plus_dict,gene_minus_dict,vx,criterion_len,gap_j):
    q_l = float(pj[6]);q_r = float(pj[7]);t_l = min(float(pj[8]),float(pj[9]));t_r = max(float(pj[8]),float(pj[9]))
    chr_p = pj[1]
    mid_t = (t_l + t_r) / 2;len_seq = float(pj[3]);map_identity = float(pj[2])
    if abs(mid_t - fus_l_mid) > 50 and abs(mid_t - fus_r_mid) > 50 and map_identity > 85:
        if (q_r >= 100-gap_j and q_r <= 100+gap_j) or (q_l >= 100-gap_j and q_l <= 100+gap_j):
            if map_identity ==100 and len_seq > criterion_len:
                type2_erro_list.append([vx] + pj+['SameIdentityfalse_1'])
            elif map_identity >= min_score_filter and len_seq >= 100-gap_j-18 and len_seq<=100+gap_j+18:
                type2_erro_list.append([vx] + pj + ['SameIdentityfalse_2'])
            elif len_seq >20 and gene_plus_dict.has_key(chr_p) and gene_minus_dict.has_key(chr_p):
                type2_erro_list = filter_same_gene(vx, pj,nchr,t_l, t_r, gene_plus_dict, gene_minus_dict,type2_erro_list,fus_l_l,fus_l_r,fus_r_l, fus_r_r)
        elif q_l <= 60 and q_r >= 140 and len_seq > criterion_len and map_identity > 95:
            type2_erro_list.append([vx] + pj+["Co_linear explanation"])
    return type2_erro_list

def filter_same_gene(vx, pj,nchr,t_l, t_r,gene_plus_dict, gene_minus_dict, type2_erro_list,fus_l_l,fus_l_r,fus_r_l, fus_r_r):
    gene_fus_l=gene_fus_r = 0;gene_q = 1;chr_t=pj[1]
    gene_list = gene_plus_dict[chr_t]+ gene_minus_dict[chr_t]
    for p_e1 in gene_list:
        if nchr == p_e1[0] and (min(fus_l_r,p_e1[2])-max(fus_l_l,p_e1[1]))>20:
            gene_fus_l = p_e1[3:]
            break
        else:
            gene_fus_l = ['none1','none1']
    for p_e2 in gene_list:
        if nchr == p_e2[0] and (min(fus_r_r,p_e2[2])-max(fus_r_l,p_e2[1]))>20:
            gene_fus_r = p_e2[3:]
            break
        else:
            gene_fus_r = ['none2', 'none2']
    for p_e3 in gene_list:
        if chr_t == p_e3[0] and (min(t_r,p_e3[2])-max(t_l,p_e3[1]))>20:
            gene_q = p_e3[3:]
            break
        else:
            gene_q=['none3', 'none3']
    if gene_q[0]==gene_fus_l[0] or gene_q[0]==gene_fus_r[0] or gene_fus_r[0]==gene_fus_l[0]:
        type2_erro_list.append([vx] + pj+['SameGeneFalse'])
    else:
        if gene_q[1]!='none':
            if gene_q[1]==gene_fus_l[1] or gene_q[1]==gene_fus_r[1] or gene_fus_r[1]==gene_fus_l[1] :
                if gene_fus_r[1]!="none":
                    type2_erro_list.append([vx] + pj + ['SameGeneFalse'])
    return type2_erro_list

def co_fus_filter(bl8file,typeI_fusfile,prefix):
    fr_bl8 = [r1.strip().split("\t") for r1 in open(bl8file).readlines()]
    fr_fus = [r1.strip().split("\t") for r1 in open(typeI_fusfile).readlines()]
    n_fus=len(fr_fus);all_fus={}
    list_fus=[];filtered_fus=[]
    for ap in range(n_fus):
        pp=fr_fus[ap][0]+'_'+fr_fus[ap][1]+"_"+fr_fus[ap][2]
        all_fus[pp]=1
    for bp in fr_bl8:
        pp1=bp[0]
        if all_fus.has_key(pp1):
            all_fus[pp1]=all_fus[pp1]+1
    for (k,v) in all_fus.items():
        if v==1:
            list_fus.append(k)
    for ap2 in range(n_fus):
        pp3=fr_fus[ap2][0]+'_'+fr_fus[ap2][1]+"_"+fr_fus[ap2][2]
        if pp3 in list_fus:
            filtered_fus.append(fr_fus[ap2])
    t2_list = []
    for it2 in filtered_fus:
        ss2 = '\t'.join(it2[:5]+it2[-4:])
        t2_list.append(ss2)
    if len(t2_list)>0:
        lf2 = open(prefix + "_typeIIfiltered.txt", 'w+')
        lf2.writelines("\n".join(t2_list) + "\n")
        lf2.close()
    return all_fus,list_fus

def main(argv):
    global min_score_filter
    t1 = time.time()
    typeI_fusfile =''
    bl8_format=''
    genefile=''
    min_score_filter=''
    try:
        opts, args = getopt.getopt(argv, "hi:b:e:m:", ["ifile=","bfile=","efile=","mfile="])
    except getopt.GetoptError:
        print'python typeII_erro_fusions_filter -i <typeI erro filtered fusions file> -b <bl8_format file> -e <genefile> -m <min_score_filter> '
        sys.exit(4)
    for opt, arg in opts:
        if opt == '-h':
            print'python typeII_erro_fusions_filter -i <typeI erro filtered fusions file> -b <bl8_format file> -e <genefile> -m <min_score_filter> '
            sys.exit()
        elif opt in ("-i", "--ifile"):
            typeI_fusfile = arg
        elif opt in ("-b", "--bfile"):
            bl8_format = arg
        elif opt in ("-e", "--efile"):
            genefile = arg
        elif opt in ("-m", "--mfile"):
            min_score_filter = arg
    min_score_filter=float(min_score_filter)
    prefix=typeI_fusfile[:-4]
    gap_j=4
    bl8_fus_filter(bl8_format,typeI_fusfile, genefile,80,prefix,gap_j)
    co_fus_filter(prefix+"_wrong_bl8",typeI_fusfile,prefix)
    os.system("rm " + prefix + "*exon*")
    os.system("rm " + prefix + "*bl8*")
    t2=time.time()
    print"The time needed in typeII error filter for fusions is:",t2-t1

if __name__ == "__main__":
   main(sys.argv[1:])


