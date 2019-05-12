import sys, getopt
import subprocess

def junc_filter(juncfile,prefix):
    fr_junc = [r1.strip().split("\t") for r1 in open(juncfile).readlines()]
    n=len(fr_junc);junc_leftexon=[];junc_rightexon=[]
    for i in range(n):
        v1=100
        v2=100
        if fr_junc[i][3] == '+':
            ple1=str(int(float(fr_junc[i][1]) +1 -v1))
            ple2=str(int(float(fr_junc[i][1]) +1))
            pre1=str(int(float(fr_junc[i][2])))
            pre2=str(int(float(fr_junc[i][2])+v2))
            junc_leftexon.append([fr_junc[i][0]] + [ple1, ple2, fr_junc[i][7], fr_junc[i][4], fr_junc[i][3]])
            junc_rightexon.append([fr_junc[i][0]] + [pre1, pre2, fr_junc[i][7], fr_junc[i][4], fr_junc[i][3]])
        elif fr_junc[i][3]=='-':
            ple1 = str(int(float(fr_junc[i][1]) + 1 - v2))
            ple2 = str(int(float(fr_junc[i][1]) + 1))
            pre1 = str(int(float(fr_junc[i][2])))
            pre2 = str(int(float(fr_junc[i][2]) + v1))
            junc_leftexon.append([fr_junc[i][0]] + [ple1,ple2,fr_junc[i][7],fr_junc[i][4],fr_junc[i][3]])
            junc_rightexon.append([fr_junc[i][0]] + [pre1,pre2,fr_junc[i][7],fr_junc[i][4],fr_junc[i][3]])
    lw1 = open(prefix+"_leftexon.bed", 'w+')
    t1_list=[]
    for iterm1 in junc_leftexon:
        s1 = '\t'.join(iterm1)
        t1_list.append(s1)
    lw1.writelines("\n".join(t1_list) + "\n")
    lw1.close()
    lw2 = open(prefix+"_rightexon.bed", 'w+')
    t2_list = []
    for iterm2 in junc_rightexon:
        s2 = '\t'.join(iterm2)
        t2_list.append(s2)
    lw2.writelines("\n".join(t2_list) + "\n")
    lw2.close()
    return fr_junc, n
def add_junc_sequence(fr_junc,n,juncfile,fastafile):
    subprocess.call("bedtools getfasta -fi "+fastafile+ " -bed "+ juncfile + "_leftexon.bed -s -fo "+juncfile+"_leftexon.fa",shell = True)
    subprocess.call("bedtools getfasta -fi "+fastafile+ " -bed  "+ juncfile + "_rightexon.bed -s -fo "+juncfile+ "_rightexon.fa",shell = True)
    fr_le = [r1.strip().split("\t") for r1 in open(juncfile+"_leftexon.fa").readlines()]
    fr_re = [r1.strip().split("\t") for r1 in open(juncfile+ "_rightexon.fa").readlines()]
    jun_l_b = [r1.strip().split("\t") for r1 in open(juncfile+"_leftexon.bed").readlines()]
    jun_r_b = [r1.strip().split("\t") for r1 in open(juncfile+"_rightexon.bed").readlines()]
    seqlist=[];junc_segment=[]
    for p in range(n):
        if fr_junc[p][3]=='+':
            par1=fr_le[p * 2 + 1]+fr_re[p * 2 + 1]
            seqtmp = fr_junc[p][:8] + par1
            seqlist.append(seqtmp)
            jf1 = [">" + fr_junc[p][0] + "_" + fr_junc[p][1] + "_" + fr_junc[p][2]+"*"+ jun_l_b[p][1]+"~"+ jun_l_b[p][2]+ "~"+jun_r_b[p][1]+"~"+ jun_r_b[p][2]+"_"+fr_junc[p][3]]
            jf2 = ['.'.join(fr_le[p * 2 + 1]) + '.'.join(fr_re[p * 2 + 1])]
            junc_segment.append(jf1)
            junc_segment.append(jf2)
        elif fr_junc[p][3] == '-':
            par2 =  fr_re[p * 2 + 1]+fr_le[p * 2 + 1]
            seqtmp = fr_junc[p][:8] + par2
            seqlist.append(seqtmp)
            jf1 = [">" + fr_junc[p][0] + "_" + fr_junc[p][1] + "_" + fr_junc[p][2] + "*"+ jun_r_b[p][1]+"~"+ jun_r_b[p][2]+ "~"+jun_l_b[p][1]+"~"+ jun_l_b[p][2]+"_"+fr_junc[p][3]]
            jf2 = ['.'.join(fr_re[p * 2 + 1]) + '.'.join(fr_le[p * 2 + 1])]
            junc_segment.append(jf1)
            junc_segment.append(jf2)
    lj2 = open(juncfile+"_exonseq.fa", 'w+')
    tt_list=[]
    for itt in junc_segment:
        sj2 = '\t'.join(itt)
        tt_list.append(sj2)
    lj2.writelines("\n".join(tt_list) + "\n")
    lj2.close()

def main(argv):
    global softwaredir
    candidatefile = ''
    fastafile = ''
    fasta2bitfile=''
    softwaredir = ""

    try:
        opts, args = getopt.getopt(argv, "hi:f:b:s:", ["ifile=","ffile=","bfile=","sfile="])
    except getopt.GetoptError:
        print 'python get_typeII_junctions_seq.py -i <typeI erro filtered junctions file>  -f <fastafile> -b <fasta2bitfile> -s <softwaredir>  '
        sys.exit(4)
    for opt, arg in opts:
        if opt == '-h':
            print 'python get_typeII_junctions_seq.py -i <typeI erro filtered junctions file> -f <fastafile> -b <fasta2bitfile> -s <softwaredir>  '
            sys.exit()
        elif opt in ("-i", "--ifile"):
            candidatefile = arg
        elif opt in ("-f", "--ffile"):
            fastafile = arg
        elif opt in ("-b", "--bfile"):
            fasta2bitfile = arg
        elif opt in ("-s", "--sfile"):
            softwaredir = arg

    prefix=candidatefile[:-4]
    fr_junc, n=junc_filter(candidatefile,prefix)
    add_junc_sequence(fr_junc, n, prefix, fastafile)
    subprocess.call(softwaredir+"/tools/blat -stepSize=5 -repMatch=2253 -minScore=20 -minIdentity=80 -out=blast8 " + fasta2bitfile + " " + prefix + "_exonseq.fa " + prefix + "_falsemapping_bl8.out",shell = True)



if __name__ == "__main__":
   main(sys.argv[1:])






