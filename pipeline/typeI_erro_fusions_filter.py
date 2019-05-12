
import sys, getopt,time,Bio
from Bio import pairwise2
#get the total massage of fusion.out after pre-filter
def fusionseq(raw_fusion,fusionfilteredfile):
    fr_raw = [r1.strip().split("\t") for r1 in open(raw_fusion).readlines()]
    fr_filtered = [r1.strip().split("\t") for r1 in open(fusionfilteredfile).readlines()]
    mm1=len(fr_raw);mm2=len(fr_filtered);filter_dict={};fusion_list=[]
    for i in range(mm2):
        par1=fr_filtered[i][0]+"_"+fr_filtered[i][1]+"_"+fr_filtered[i][2]
        filter_dict[par1]=1
    for j in range(mm1):
        p1= fr_raw[j][0].split("-")[0]
        p2= fr_raw[j][0].split("-")[1]
        if p1==p2:
            par2=p1+"_"+fr_raw[j][1]+"_"+fr_raw[j][2]
            if filter_dict.has_key(par2):
                pp=[p1]+fr_raw[j][1:]
                fusion_list.append(pp)
    return fusion_list
# get the trans-splicing pattern fusion and filtering the TypeI erro
def fusion_filter(fusion_list,file_prefix):
    nf=len(fusion_list);fusion_pattern = ["GT-AG","GC-AG","CT-AC","CT-GC"]
    typeIErr_filter_fusion=[]
    for i in range(nf):
        length=15
        shift_l,shift_r=find_real_site(fusion_list[i])
        total_seq=str.upper(fusion_list[i][14].split(" ")[0]+fusion_list[i][14].split(" ")[1]+fusion_list[i][16].split(" ")[0]+fusion_list[i][16].split(" ")[1])
        le=total_seq[50+shift_l-length:50+shift_l]
        li=total_seq[50+shift_l:50+shift_l+length]
        ri = total_seq[150+shift_r-length:150+shift_r]
        re = total_seq[150+shift_r:150+shift_r+length]
        pattern=li[:2]+"-"+ri[-2:]
        if pattern in fusion_pattern:
            map_score1 = pairwise2.align.globalxx(le, ri, score_only=True)
            map_score2 = pairwise2.align.globalxx(li, re, score_only=True)
            if map_score1 <= 12 and map_score2 <= 12:
                typeIErr_filter_fusion.append(fusion_list[i]+[le,li,ri,re])
    typeIErr_filter_fusion2 = mer(typeIErr_filter_fusion)
    t1_list=[]
    for it in typeIErr_filter_fusion2:
        ss = '\t'.join(it)
        t1_list.append(ss)
    if len(t1_list) >0:
        lw = open(file_prefix + "_typeIfiltered.txt", 'w+')
        lw.writelines("\n".join(t1_list) + "\n")
        lw.close()
    return typeIErr_filter_fusion
def mer(finall_list): #filtered replicate fusions
    dict_select={};select_list=[]
    for i in range(len(finall_list)):
        t_chr= finall_list[i][0]
        t_l=float(finall_list[i][1])
        t_r=float(finall_list[i][2])
        rep=False
        for j in dict_select.keys():
            s_chr=j.split("_")[0]
            if s_chr==t_chr:
                s_l =float(j.split("_")[1])
                s_r =float(j.split("_")[2])
                if abs(s_l-t_l)<=6 and abs(s_r-t_r)<=6: #select rep fusions with more repport reads
                    if float(dict_select[j][4]) >=float(finall_list[i][4]):
                        rep=True
                        break
                    else:
                        del dict_select[j]
                        par_t=t_chr+"_"+finall_list[i][1]+"_"+finall_list[i][2]
                        dict_select[par_t]=finall_list[i]
                        rep = True
                        break
        if not rep:
            par=t_chr+"_"+finall_list[i][1]+"_"+finall_list[i][2]
            dict_select[par]=finall_list[i]
    for k in  dict_select.keys():
        select_list.append(dict_select[k])
    return select_list
##find the real splicing site
def find_real_site(rp):
    splic_l=list(rp[14].split(' ')[0][-6:]+rp[14].split(' ')[1][:10])
    splic_r=list(rp[16].split(' ')[0][-6:]+rp[16].split(' ')[1][:10])
    standard_l=6;standard_r=5;shift_l=0;shift_r=0
    fac=False
    for i_l in range(2,15):
        patern_l = str.upper(splic_l[i_l] + splic_l[i_l + 1])
        if patern_l == 'GT'or patern_l=="GC":
            index_l=i_l
            shift_l = index_l - standard_l
            patern_r=str.upper(splic_r[standard_r + shift_l-1]+splic_r[standard_r+shift_l])
            if patern_r == 'AG':
                index_r = standard_r+shift_l
                if shift_l>0:
                    shift_seq_l="".join(splic_l[standard_l:index_l])
                    shift_seq_r="".join(splic_r[standard_r+1:index_r+1])
                    if shift_seq_l ==shift_seq_r:
                        fac = True
                        break
                elif shift_l<0:
                    shift_seq_l = "".join(splic_l[index_l:standard_l])
                    shift_seq_r = "".join(splic_r[index_r + 1:standard_r + 1])
                    if shift_seq_l == shift_seq_r:
                        fac = True
                        break
                elif shift_l==0:
                    break
        elif patern_l == 'CT':
            index_l=i_l
            shift_l = index_l - standard_l
            patern_r=str.upper(splic_r[standard_r + shift_l-1]+splic_r[standard_r+shift_l])
            if patern_r == 'AC'or patern_r=="GC":
                index_r = standard_r+shift_l
                if shift_l>0:
                    shift_seq_l="".join(splic_l[standard_l:index_l])
                    shift_seq_r="".join(splic_r[standard_r+1:index_r+1])
                    if shift_seq_l ==shift_seq_r:
                        fac = True
                        break
                elif shift_l<0:
                    shift_seq_l = "".join(splic_l[index_l:standard_l])
                    shift_seq_r = "".join(splic_r[index_r + 1:standard_r + 1])
                    if shift_seq_l == shift_seq_r:
                        fac = True
                        break
                elif shift_l==0:
                    break
    if fac:
        shift_l=shift_l
        shift_r=shift_l
    else:
        shift_l=0
        fac2 = False
        fusion_l =  str.upper(rp[14].split(" ")[1][:4])
        fusion_r =  str.upper(rp[16].split(" ")[0][-4:])
        for q in range(3):
            p_l = fusion_l[q:q + 2]
            if p_l == 'GT' or p_l == 'GC':
                for p in range(3):
                    p_r = fusion_r[p:p + 2]
                    if p_r == 'AG':
                        shift_l=q
                        shift_r=p-2
                        fac2 = True
                        break
            elif p_l == 'CT':
                for p in range(3):
                    p_r = fusion_r[p:p + 2]
                    if p_r == 'AC' or p_r == "GC":
                        shift_l = q
                        shift_r=p-2
                        fac2=True
                        break
            if fac2:
                break
    return shift_l,shift_r
def main(argv):
    t1 = time.time()
    raw_fusion = ''
    fusionfilteredfile = ''
    try:
        opts, args = getopt.getopt(argv, "hr:i:", ["rfile=", "ifile="])
    except getopt.GetoptError:
        print'python typeI_erro_fusions_filter.py -r <raw_fusion> -i <prefiltered fusions file>  '
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print'python typeI_erro_fusions_filter.py -r <raw_fusion> -i <prefiltered fusions file>  '
            sys.exit()
        elif opt in ("-r", "--rfile"):
            raw_fusion = arg
        elif opt in ("-i", "--ifile"):
            fusionfilteredfile = arg
    file_prefix=fusionfilteredfile[:-4]
    fusion_list=fusionseq(raw_fusion, fusionfilteredfile)
    fusion_filter(fusion_list,file_prefix)
    t2 = time.time()
    print"The time needed in typeI error filter for fusions is:", t2 - t1

if __name__ == "__main__":
    main(sys.argv[1:])










