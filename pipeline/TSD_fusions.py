import re,os,getopt,sys
import subprocess

def filterjunc(configurefile,outputdir):
    fr_con = [r1.strip().split("\t") for r1 in open(configurefile).readlines()]
    genefile = ''
    gapfile = ''
    fastafile = ''
    fasta2bitfile = ''
    fusions_inputfile = ''
    fusions_repswitch = ''
    fusions_repfile1 = ''
    fusions_repfile2 = ''
    min_sup=''
    min_score_filter=''
    Softwaredir = ''
    for pr in fr_con:
        if re.search("=", pr[0]):
            vx=pr[0].split("=")[0].split("#")[0].strip(" ")
            vy=pr[0].split("=")[1].split("#")[0].strip(" ")
            if re.search("genefile",vx):
                genefile=vy
            elif re.search("gapfile",vx):
                gapfile=vy
            elif re.search("fastafile",vx):
                fastafile=vy
            elif re.search("fasta2bitfile",vx):
                fasta2bitfile=vy
            elif re.search("fusions_inputfile",vx):
                fusions_inputfile=vy
            elif re.search("fusions_repswitch",vx):
                fusions_repswitch=vy
            elif re.search("fusions_repfile1",vx):
                fusions_repfile1=vy
            elif re.search("fusions_repfile2",vx):
                fusions_repfile2=vy
            elif re.search("min_score_filter",vx):
                min_score_filter=vy
            elif re.search("min_sup",vx):
                min_sup=vy
            elif re.search("Softwaredir",vx):
                Softwaredir=vy
    print "******Filtering with TSD-fusions******"
    subprocess.call("python "+Softwaredir+"/pipeline/pre_filter_fusions.py -i " + fusions_inputfile + " -g  " + genefile +\
                    " -w " + fusions_repswitch + " -r " + fusions_repfile1 + " -t " + fusions_repfile2 + " -p " + gapfile+ " -m "+min_sup+" -o "+outputdir,shell = True)
    fusions_inputfile2 = outputdir + "/" + [x for x in fusions_inputfile.split("/") if x != ""][-1]
    if os.path.isfile(fusions_inputfile2+ "_prefiltered.txt"):
        subprocess.call("python "+Softwaredir+"/pipeline/typeI_erro_fusions_filter.py -r " + fusions_inputfile \
                        + " -i " + fusions_inputfile2 + "_prefiltered.txt",shell = True)
        if os.path.isfile(fusions_inputfile2 + "_prefiltered_typeIfiltered.txt"):
            subprocess.call("python " + Softwaredir + "/pipeline/get_typeII_fusions_seq.py -i " + fusions_inputfile2 + \
                            "_prefiltered_typeIfiltered.txt -f " + fastafile + " -b " + fasta2bitfile + " -s " + Softwaredir, shell=True)
            subprocess.call("python " + Softwaredir + "/pipeline/typeII_erro_fusions_filter.py -i " + fusions_inputfile2 \
                            + "_prefiltered_typeIfiltered.txt -b " + fusions_inputfile2 + \
                "_prefiltered_typeIfiltered_falsemapping_bl8.out -e " + genefile + " -m " + min_score_filter,shell=True)
        else:
            print "There remains no candidate after TypeI erro filter step"
    else:
        print "There remains no candidate after Pre-filer step"

def main(argv):
    configurefile=""
    outputdir=""
    try:
        opts, args = getopt.getopt(argv, "hi:o:", ["ifile=","odir="])
    except getopt.GetoptError:
        print 'python TSD_fusions.py -i <configurefile> -o <outputdir> '
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print 'python TSD_fusions.py -i <configurefile>  -o <outputdir>'
            sys.exit()
        elif opt in ("-i", "--ifile"):
            configurefile = arg
        elif opt in ("-o", "--odir"):
            outputdir = arg
    filterjunc(configurefile,outputdir)

if __name__ == "__main__":
    main(sys.argv[1:])