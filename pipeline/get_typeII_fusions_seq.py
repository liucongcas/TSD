import sys, getopt
import subprocess
def fusion_adding(fusionTypeIfile,file_prefix):
    fr_fusion = [r1.strip().split("\t") for r1 in open(fusionTypeIfile).readlines()]
    nf=len(fr_fusion);fus_leftexon=[];fus_rightexon=[]
    for pr in range(nf):
        l_length=100
        r_length=100
        f_l = int(float(fr_fusion[pr][1]))
        f_r = int(float(fr_fusion[pr][2]))
        shift=find_real_site(fr_fusion[pr])
        if fr_fusion[pr][3] == 'ff':
            f_l=f_l+shift
            f_r=f_r+shift
            l_e_s=f_l-l_length+1
            l_e_e=f_l+1
            r_e_s=f_r
            r_e_e=f_r+r_length
            fus_leftexon.append([fr_fusion[pr][0]] + [l_e_s, l_e_e, fr_fusion[pr][3], fr_fusion[pr][8], '+'])
            fus_rightexon.append([fr_fusion[pr][0]] + [r_e_s, r_e_e, fr_fusion[pr][3], fr_fusion[pr][9], '+'])
        if fr_fusion[pr][3] == 'fr':
            f_l = f_l + shift
            f_r = f_r - shift
            l_e_s = f_l - l_length + 1
            l_e_e = f_l + 1
            r_e_s = f_r - r_length + 1
            r_e_e = f_r + 1
            fus_leftexon.append([fr_fusion[pr][0]] + [l_e_s, l_e_e, fr_fusion[pr][3], fr_fusion[pr][8], '+'])
            fus_rightexon.append([fr_fusion[pr][0]] + [r_e_s, r_e_e, fr_fusion[pr][3], fr_fusion[pr][9],'-'])
        if fr_fusion[pr][3] == 'rf':
            f_l = f_l - shift
            f_r = f_r + shift
            l_e_s = f_l
            l_e_e = f_l + l_length
            r_e_s = f_r
            r_e_e = f_r + r_length
            fus_leftexon.append([fr_fusion[pr][0]] + [l_e_s, l_e_e, fr_fusion[pr][3], fr_fusion[pr][8],'-'])
            fus_rightexon.append([fr_fusion[pr][0]] + [r_e_s, r_e_e, fr_fusion[pr][3], fr_fusion[pr][9],'+'])
        if fr_fusion[pr][3] == 'rr':
            f_l = f_l - shift
            f_r = f_r - shift
            l_e_s = f_l
            l_e_e = f_l + l_length
            r_e_s = f_r - r_length + 1
            r_e_e = f_r + 1
            fus_leftexon.append([fr_fusion[pr][0]] + [l_e_s, l_e_e, fr_fusion[pr][3], fr_fusion[pr][8],'-'])
            fus_rightexon.append([fr_fusion[pr][0]] + [r_e_s, r_e_e, fr_fusion[pr][3], fr_fusion[pr][9],'-'])
    lf1 = open(file_prefix+"_leftexon.bed", 'w+')
    t1_list=[]
    for lt1 in fus_leftexon:
        s1 = '\t'.join(map(str,lt1))
        t1_list.append(s1)
    lf1.writelines("\n".join(t1_list) + "\n")
    lf1.close()
    lf4 = open(file_prefix+"_rightexon.bed", 'w+')
    t2_list = []
    for lt4 in fus_rightexon:
        s4 = '\t'.join(map(str,lt4))
        t2_list.append(s4)
    lf4.writelines("\n".join(t2_list) + "\n")
    lf4.close()
    return fr_fusion, nf
def find_real_site(rp):
    splic_l=list(rp[14].split(' ')[0][-6:]+rp[14].split(' ')[1][:10])
    splic_r=list(rp[16].split(' ')[0][-6:]+rp[16].split(' ')[1][:10])
    standard_l=6;standard_r=5;shift_l=0;fac=False
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
            else:
                shift_l=0
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
    else:
        shift_l=0
    return shift_l

def add_fus_sequence(fr_fusion,nf,fastafile,file_prefix):
    subprocess.call("bedtools getfasta -fi " +fastafile+ " -bed "+file_prefix+ "_leftexon.bed -s -fo "+file_prefix+"_leftexon.fa",shell = True)
    subprocess.call("bedtools getfasta -fi " +fastafile+ " -bed "+file_prefix+"_rightexon.bed -s -fo "+file_prefix+"_rightexon.fa",shell = True)
    fr_le = [r1.strip().split("\t") for r1 in open(file_prefix+"_leftexon.fa").readlines()]
    fr_re = [r1.strip().split("\t") for r1 in open(file_prefix+"_rightexon.fa").readlines()]
    fus_le_b=[r1.strip().split("\t") for r1 in open(file_prefix+"_leftexon.bed").readlines()]
    fus_re_b = [r1.strip().split("\t") for r1 in open(file_prefix+"_rightexon.bed").readlines()]
    fus_segment=[]
    for qf in range(nf):
        ff1=[">"+fr_fusion[qf][0]+"_"+fr_fusion[qf][1]+"_"+fr_fusion[qf][2] + "*" + fus_le_b[qf][1]+"~"+fus_le_b[qf][2]+"~"+ fus_re_b[qf][1]+"~"+fus_re_b[qf][2]+"_" + fr_fusion[qf][3]]
        ff2=['.'.join(fr_le[qf * 2 + 1])+'.'.join(fr_re[qf * 2 + 1])]
        fus_segment.append(ff1)
        fus_segment.append(ff2)
    tt_list=[]
    for it in fus_segment:
        ss='\t'.join(it)
        tt_list.append(ss)
    lf = open(file_prefix + "_exonseq.fa", 'w+')
    lf.writelines("\n".join(tt_list) + "\n")
    lf.close()

def main(argv):
    global softwaredir
    fusionTypeIfile = ''
    fastafile = ''
    fasta2bitfile=''
    softwaredir=""
    try:
        opts, args = getopt.getopt(argv, "hi:f:b:s:", ["ifile=","ffile=","bfile=","sfile="])
    except getopt.GetoptError:
        print 'python get_typeII_fusions_seq.py -i <typeI erro filtered fusions file>  -f <fastafile> -b <fasta2bitfile> -s <softwaredir>'
        sys.exit(4)
    for opt, arg in opts:
        if opt == '-h':
            print 'python get_typeII_fusions_seq.py -i <typeI erro filtered fusions file> -f <fastafile> -b <fasta2bitfile> -s <softwaredir> '
            sys.exit()
        elif opt in ("-i", "--ifile"):
            fusionTypeIfile = arg
        elif opt in ("-f", "--ffile"):
            fastafile = arg
        elif opt in ("-b", "--bfile"):
            fasta2bitfile = arg
        elif opt in ("-s", "--sfile"):
            softwaredir = arg
    file_prefix = fusionTypeIfile[:-4]
    fr_fusion, nf=fusion_adding(fusionTypeIfile,file_prefix)
    add_fus_sequence(fr_fusion, nf,fastafile,file_prefix)
    subprocess.call(softwaredir+"/tools/blat -stepSize=5 -repMatch=2253 -minScore=20 -minIdentity=80 -out=blast8 "+ fasta2bitfile + " "+file_prefix+"_exonseq.fa "+file_prefix+"_falsemapping_bl8.out",shell = True)
if __name__ == "__main__":
   main(sys.argv[1:])










