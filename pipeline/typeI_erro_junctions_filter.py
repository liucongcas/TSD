
import sys, getopt,os,time
from Bio import pairwise2

def addseq(seq_list,fr_junc,m):
    seqlist=[]
    for p in range(m):
        seqtmp = fr_junc[p] + seq_list[p*2+1]
        seqlist.append(seqtmp)
    return seqlist

def add_total_seq(juncfile,fastafile):
    os.system("awk '{print $1"'"\t"'"$2-14"'"\t"'"$2+1"'"\t"'"$3"'"\t"'"$5"'"\t"'"$4}' " + juncfile + " > " + juncfile + "_leftexon.bed")
    os.system("awk '{print $1"'"\t"'"$2+1"'"\t"'"$2+16"'"\t"'"$3"'"\t"'"$5"'"\t"'"$4}' " + juncfile + " > " + juncfile + "_leftintron.bed")
    os.system("awk '{print $1"'"\t"'"$3"'"\t"'"$3+15"'"\t"'"$2"'"\t"'"$5"'"\t"'"$4}' " + juncfile + " > " + juncfile + "_rightexon.bed")
    os.system("awk '{print $1"'"\t"'"$3-15"'"\t"'"$3"'"\t"'"$2"'"\t"'"$5"'"\t"'"$4}' " + juncfile + " > " + juncfile + "_rightintron.bed")
    os.system("bedtools getfasta -fi " + fastafile + " -bed " + juncfile + "_leftexon.bed -s -fo " + juncfile + "_leftexon.fa")
    os.system("bedtools getfasta -fi " + fastafile + " -bed " + juncfile + "_leftintron.bed -s -fo " + juncfile + "_leftintron.fa")
    os.system( "bedtools getfasta -fi " + fastafile + " -bed " + juncfile + "_rightexon.bed -s -fo " + juncfile + "_rightexon.fa")
    os.system( "bedtools getfasta -fi " + fastafile + " -bed " + juncfile + "_rightintron.bed -s -fo " + juncfile + "_rightintron.fa")
    fr_junc = [r.strip().split("\t") for r in open(juncfile).readlines()]
    fr_le = [r1.strip().split("\t") for r1 in open(juncfile + "_leftexon.fa").readlines()]
    fr_li = [r2.strip().split("\t") for r2 in open(juncfile + "_leftintron.fa").readlines()]
    fr_re = [r1.strip().split("\t") for r1 in open(juncfile + "_rightexon.fa").readlines()]
    fr_ri = [r1.strip().split("\t") for r1 in open(juncfile + "_rightintron.fa").readlines()]
    m = len(fr_junc)
    add_leseq = addseq(fr_le, fr_junc, m)
    add_liseq = addseq(fr_li, add_leseq, m)
    add_riseq = addseq(fr_ri, add_liseq, m)
    add_reseq = addseq(fr_re, add_riseq, m)
    return add_reseq
def junc_filter(fr_junc,prefix):
    n=len(fr_junc);filter_junc=[];typeIErr_filter_junc=[]
    junc_pattern1 = ["GT-AG", "GC-AG", "AT-AC"]
    junc_pattern2 = ["AG-GT", "AG-GC", "AC-AT"]
    for i in range(n):
        if fr_junc[i][3]=='+':
            par1 = str.upper(fr_junc[i][9][0:2])+"-"+str.upper(fr_junc[i][10][13:15])
            if par1 in junc_pattern1:
                filter_junc.append(fr_junc[i])
        elif fr_junc[i][3] == '-':
            par2 = str.upper(fr_junc[i][9][13:15]) + "-" + str.upper(fr_junc[i][10][0:2])
            if par2 in junc_pattern2:
                filter_junc.append(fr_junc[i])
    m = len(filter_junc)
    for j in range(m):
        map_score1 = pairwise2.align.globalxx(filter_junc[j][8], filter_junc[j][10], score_only=True)
        map_score2 = pairwise2.align.globalxx(filter_junc[j][9], filter_junc[j][11], score_only=True)
        if map_score1 <=12 and map_score2<=12:
            typeIErr_filter_junc.append(filter_junc[j])
    t1_list=[]
    for iterm2 in typeIErr_filter_junc:
        s2 = '\t'.join(iterm2)
        t1_list.append(s2)
    if len(t1_list)>0:
        lw2 = open(prefix + "_typeIfiltered.txt", 'w+')
        lw2.writelines("\n".join(t1_list) + "\n")
        lw2.close()

def main(argv):
    t1 = time.time()
    juncfile = ''
    fastafile = ''
    global outputdir
    outputdir = ''
    try:
        opts, args = getopt.getopt(argv, "hi:f:o:", ["ifile=","ffile=","odir="])
    except getopt.GetoptError:
        print 'python typeI_erro_junctions_filter.py  -i <prefiltered junctions file> -f <fastafile> -o <outputdir>'
        sys.exit(3)
    for opt, arg in opts:
        if opt == '-h':
            print 'python typeI_erro_junctions_filter.py -i <prefiltered junctions file>  -f <fastafile> -o <outputdir>'
            sys.exit()
        elif opt in ("-i", "--ifile"):
            juncfile = arg
        elif opt in ("-f", "--ffile"):
            fastafile = arg
        elif opt in ("-o", "--odir"):
            outputdir = arg
    prefix=juncfile[:-4]
    add_reseq=add_total_seq(juncfile, fastafile)
    junc_filter(add_reseq,prefix)
    os.system("rm " + juncfile + "_*exon.fa")
    os.system("rm " + juncfile + "_*intron.fa")
    os.system("rm " + juncfile + "_*exon.bed")
    os.system("rm " + juncfile + "_*intron.bed")
    t2=time.time()
    print "The time needed in typeI error filter for junctions is:",t2-t1

if __name__ == "__main__":
   main(sys.argv[1:])






