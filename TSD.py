import os,getopt,sys
import subprocess,re
def combin_results():
    fr_con = [r1.strip().split("\t") for r1 in open(configurefile).readlines()]
    global Softwaredir
    global genefile
    global junctions_inputfile
    global fusions_inputfile
    global SV_dis
    global blank_len
    global delSV_switch
    global grofiledir
    global gro_switch
    junctions_inputfile=''
    fusions_inputfile=''
    gro_switch = ''
    grofiledir = ''
    delSV_switch = ''
    SVVCF = ''
    blank_len=''
    Softwaredir = ''
    SV_dis=''
    for pr in fr_con:
        if re.search("=", pr[0]):
            vx = pr[0].split("=")[0].split("#")[0].strip(" ")
            vy = pr[0].split("=")[1].split("#")[0].strip(" ")
            if re.search("genefile", vx):
                genefile = vy
            elif re.search("junctions_inputfile",vx):
                junctions_inputfile=vy
            elif re.search("fusions_inputfile",vx):
                fusions_inputfile=vy
            elif re.search("gro_switch", vx):
                gro_switch = vy
            elif re.search("grofiledir", vx):
                grofiledir = vy
            elif re.search("delSV_switch", vx):
                delSV_switch = vy
            elif re.search("SVVCF", vx):
                SVVCF = vy
            elif re.search("blank_len", vx):
                blank_len = vy
            elif re.search("Softwaredir", vx):
                Softwaredir = vy
            elif re.search("SV_dis", vx):
                SV_dis = vy
    subprocess.call("python "+Softwaredir+"/pipeline/TSD_junctions.py -i "+configurefile+" -o "+ outputdir,shell = True)
    subprocess.call("python "+Softwaredir+"/pipeline/TSD_fusions.py -i "+configurefile+" -o "+ outputdir,shell = True)
    filelist=os.listdir(outputdir)
    inputfile="iTSEs_candidates"
    choosingfile_typeIIFiltered = [file for file in filelist if re.search("typeIIfiltered.txt", file)]
    if len(choosingfile_typeIIFiltered)>0:
        merge_typeIIFiltered(filelist, inputfile)
        annotation(outputdir + "/" + inputfile)
        file_For_Readthrough=outputdir + "/" +inputfile + "_annotated"
        if gro_switch == "Y":
            subprocess.call("python " + Softwaredir + "/pipeline/gro_readthrough_filter.py -i " + file_For_Readthrough \
                            +" -d " + grofiledir + " -t "+blank_len, shell=True)
            file_For_SV = file_For_Readthrough + "_GRO_RTfiltered.txt"
            if delSV_switch == "Y":
                if os.path.isfile(file_For_SV):
                    subprocess.call("python " + Softwaredir + "/pipeline/DNAlevel_NCL_filter.py -i " + file_For_SV\
                                    + " -a " + SVVCF + " -s " + SV_dis+ " -o " + outputdir, shell=True)
                else:
                    print "There remains no candidate after read-through filter step"
            else:
                subprocess.call("mv " + file_For_SV + "final_iTSEs.txt", shell=True)
        else:
            subprocess.call("python " + Softwaredir + "/pipeline/non_gro_readthrough_filter.py -i " + file_For_Readthrough + " -j " + genefile, shell=True)
            file_For_SV = file_For_Readthrough + "_nonGRO_RTfiltered.txt"
            if delSV_switch == "Y":
                if os.path.isfile(file_For_SV):
                    subprocess.call("python " + Softwaredir + "/pipeline/DNAlevel_NCL_filter.py -i " + file_For_SV\
                                    + " -a " + SVVCF + " -s " + SV_dis + " -o " + outputdir, shell=True)
                else:
                    print "There remains no candidate after read-through filter step"
            else:
                subprocess.call("mv " + file_For_SV + "final_iTSEs.txt", shell=True)
    else:
        print "There remains no candidate after TypeII erro filter step"

    if os.path.isfile(outputdir + "/final_iTSEs"):
        subprocess.call("rm " + outputdir + "/*_RTfiltered.txt", shell=True)
        subprocess.call("rm " + outputdir + "/*_prefiltered*", shell=True)

def merge_typeIIFiltered(filelist,inputfile):
    fr_typeIIfiltered = []
    choosingfile_typeIIFiltered = [file for file in filelist if re.search("typeIIfiltered.txt", file)]
    for file in choosingfile_typeIIFiltered:
        fr_typeIIfiltered += [r1.strip().split("\t") for r1 in open(outputdir + "/" + file).readlines()]
    lw = open(outputdir + "/" + inputfile, 'w+')
    t2_list = []
    for iterm in fr_typeIIfiltered:
        s = '\t'.join(iterm)
        t2_list.append(s)
    lw.writelines("\n".join(t2_list) + "\n")
    lw.close()
def annotation(annfile):
    subprocess.call("python " + Softwaredir + "/pipeline/iTSEs_annotation.py -i "+annfile+" -a "+ genefile,shell=True)


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
