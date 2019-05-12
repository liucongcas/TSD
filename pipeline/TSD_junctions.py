import re,getopt,sys,os
import subprocess


def filterjunc(configurefile,outputdir):
    fr_con = [r1.strip().split("\t") for r1 in open(configurefile).readlines()]
    genefile = ''
    gapfile = ''
    fastafile = ''
    fasta2bitfile = ''
    junctions_inputfile = ''
    junctions_repswitch = ''
    junctions_repfile1 = ''
    junctions_repfile2 = ''
    Softwaredir=''
    for pr in fr_con:
        if re.search("=",pr[0]):
            vx = pr[0].split("=")[0].split("#")[0].strip(" ")
            vy = pr[0].split("=")[1].split("#")[0].strip(" ")
            if re.search("genefile",vx):
                genefile=vy
            elif re.search("gapfile",vx):
                gapfile=vy
            elif re.search("fastafile",vx):
                fastafile=vy
            elif re.search("fasta2bitfile",vx):
                fasta2bitfile=vy
            elif re.search("junctions_inputfile",vx):
                junctions_inputfile=vy
            elif re.search("junctions_repswitch",vx):
                junctions_repswitch=vy
            elif re.search("junctions_repfile1",vx):
                junctions_repfile1=vy
            elif re.search("junctions_repfile2",vx):
                junctions_repfile2=vy
            elif re.search("Softwaredir",vx):
                Softwaredir=vy

    print "******Filtering with TSD-junctions******"
    subprocess.call("python "+Softwaredir+"/pipeline/pre_filter_junctions.py -i " + junctions_inputfile + " -g  " + genefile +\
                    " -w "+ junctions_repswitch  + " -r "+ junctions_repfile1+" -t "+junctions_repfile2 +\
                    " -p " + gapfile+ " -o "+outputdir+" -s " + Softwaredir,shell = True)
    junctions_inputfile2=outputdir+"/"+[x for x in junctions_inputfile.split("/") if x !=""][-1]
    if os.path.isfile(junctions_inputfile2+ "_prefiltered.txt"):
        subprocess.call("python "+Softwaredir+"/pipeline/typeI_erro_junctions_filter.py -i " + \
                        junctions_inputfile2+"_prefiltered.txt  -f "+ fastafile,shell = True)
        if os.path.isfile(junctions_inputfile2 + "_prefiltered_typeIfiltered.txt"):
            subprocess.call("python "+Softwaredir+"/pipeline/get_typeII_junctions_seq.py -i " + junctions_inputfile2 \
                            + "_prefiltered_typeIfiltered.txt -f" + fastafile + " -b " + fasta2bitfile +" -s " + Softwaredir,shell = True)
            subprocess.call("python "+Softwaredir+"/pipeline/typeII_erro_junctions_filter.py -i " + junctions_inputfile2 \
                            + "_prefiltered_typeIfiltered.txt -b " + junctions_inputfile2 + \
            "_prefiltered_typeIfiltered_falsemapping_bl8.out -e " + genefile,shell = True)
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
        print 'python TSD_junctions.py -i <configurefile> -o <outputdir> '
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print 'python TSD_junctions.py -i <configurefile>  -o <outputdir>'
            sys.exit()
        elif opt in ("-i", "--ifile"):
            configurefile = arg
        elif opt in ("-o", "--odir"):
            outputdir = arg
    filterjunc(configurefile,outputdir)

if __name__ == "__main__":
    main(sys.argv[1:])
