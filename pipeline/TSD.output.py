import os,getopt,sys
import subprocess,re
def combin_results():
    subprocess.call("python ./pipeline/TSD_junctions.py -i "+configurefile+" -o "+ outputdir,shell = True)
    subprocess.call("python ./pipeline/TSD_fusions.py -i "+configurefile+" -o "+ outputdir,shell = True)
    filelist=os.listdir(outputdir)
    inputfile=filelist[0].split("_")[0]
    merge_typeIIFiltered(filelist, inputfile)
    merge_read_throughFiltered(filelist, inputfile)
    choosingfile_removedSV = [file for file in filelist if re.search("removedSVs", file)]
    if len(choosingfile_removedSV)!=0:
        merge_SVfiltered(inputfile, choosingfile_removedSV)
        subprocess.call("rm " + outputdir + "/*_RTfiltered_removedSVs.txt", shell=True)
    #subprocess.call("rm " + outputdir + "/*_RTfiltered.txt", shell=True)
    #subprocess.call("rm " + outputdir + "/*_prefiltered*", shell=True)
def merge_typeIIFiltered(filelist,inputfile):
    fr_total_typeII = []
    choosingfile_typeIIFiltered = [file for file in filelist if re.search("typeIIfiltered.txt", file)]
    for file in choosingfile_typeIIFiltered:
        if file.split("_")[-1] =="typeIIfiltered.txt":
            fr_total_typeII += [r1.strip().split("\t") for r1 in open(outputdir + "/" + file).readlines()]
    lw = open(outputdir + "/" + inputfile + "_tmp", 'w+')
    for iterm in fr_total_typeII:
        s = '\t'.join(iterm)
        lw.writelines(s + "\n")
    lw.close()
    subprocess.call("sort -k1.4,1 -k2n,2 -k3n,3 " + outputdir + "/" + inputfile + "_tmp | cut -f 1-5 > " + outputdir + "/total_typeIIfiltered.txt",shell=True)
    subprocess.call("rm " + outputdir + "/" + inputfile + "_tmp", shell=True)
def merge_read_throughFiltered(filelist,inputfile):
    fr_total_read_through = []
    choosingfile_read_throughFiltered = [file for file in filelist if re.search("RTfiltered.txt", file)]
    for file in choosingfile_read_throughFiltered:
        fr_total_read_through += [r1.strip().split("\t") for r1 in open(outputdir + "/" + file).readlines()]
    lw = open(outputdir + "/" + inputfile + "_tmp", 'w+')
    for iterm in fr_total_read_through:
        s = '\t'.join(iterm)
        lw.writelines(s + "\n")
    lw.close()
    subprocess.call("sort -k1.4,1 -k2n,2 -k3n,3 " + outputdir + "/" + inputfile + "_tmp | cut -f 1-5 | sed '/^\s*$/d' > " + outputdir + "/final_iTSEs.txt",shell=True)
    subprocess.call("rm " + outputdir + "/" + inputfile + "_tmp", shell=True)
def merge_SVfiltered(inputfile,choosingfile_removedSV):
    fr_total_typeII = []
    for file in choosingfile_removedSV:
        fr_total_typeII += [r1.strip().split("\t") for r1 in open(outputdir + "/" + file).readlines()]
    lw = open(outputdir + "/" + inputfile + "_tmp", 'w+')
    for iterm in fr_total_typeII:
        s = '\t'.join(iterm)
        lw.writelines(s + "\n")
    lw.close()
    subprocess.call("sort -k1.4,1 -k2n,2 -k3n,3 " + outputdir + "/" + inputfile + "_tmp | cut -f 1-5 | sed '/^\s*$/d' > " + outputdir + "/final_iTSEs_removedSVs.txt",shell=True)
    subprocess.call("rm " + outputdir + "/" + inputfile + "_tmp", shell=True)


def main(argv):
    global configurefile
    global outputdir
    configurefile=""
    outputdir=""
    try:
        opts, args = getopt.getopt(argv, "hi:o:", ["ifile=","odir="])
    except getopt.GetoptError:
        print 'python TSD.py -i <configurefile> -o <outputdir> '
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print 'python TSD.py -i <configurefile>  -o <outputdir>'
            sys.exit()
        elif opt in ("-i", "--ifile"):
            configurefile = arg
        elif opt in ("-o", "--odir"):
            outputdir = arg
    if not os.path.exists(outputdir):
        os.mkdir(outputdir)
    combin_results()

if __name__ == "__main__":
    main(sys.argv[1:])
