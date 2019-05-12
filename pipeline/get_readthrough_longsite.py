import sys,getopt,time,os,commands
def get_time_mem():
    fr_r = [r1.strip().split("\t") for r1 in open(readthroughfile).readlines()]
    fr_g = [r1.strip().split("\t") for r1 in open(genefile).readlines()]
    list_total=[];dict_g={};list_la=[]
    for x in fr_g:
        par=x[7]
        dict_g[par]=[x[0],x[3],x[4],x[6]]
    for xr in fr_r:
        x1=xr[0].split("-")[0]
        x2=xr[0].split("-")[-1]
        if x1 in dict_g and x2 in dict_g:
            list_total.append([xr[0],dict_g[x1][0],dict_g[x1][1],dict_g[x1][3],dict_g[x2][0],dict_g[x2][2],dict_g[x2][3]])
        else:
            list_la.append(xr)
    dict_gene = splitfile(fr_g)
    list_total2=[]
    for yr in list_total:
        chr_r=yr[1]
        g1 = yr[0].split("-")[0]
        g2 = yr[0].split("-")[-1]
        s_r=float(yr[2])
        e_r=float(yr[5])
        tmp=[]
        list_gene=dict_gene[chr_r]
        for y2 in list_gene:
            s_g=float(y2[3])
            e_g=float(y2[4])
            if min(e_r,e_g) > max(s_r,s_g) and y2[7] !=g1 and y2[7]!=g2:
                tmp.append(y2)
        if tmp !=[]:
            print tmp
            q=''
            for p in tmp:
                q=q+"_".join(p)
            list_total2.append(yr+[q])
        else:
            list_total2.append(yr + ["nointermedianGene"])
    list_total3=[]
    for j in list_total2:
            list_total3.append("\t".join(map(str,j)))
    lw=open(readthroughfile+"_readthroughsites","w+")
    lw.writelines("\n".join(list_total3)+"\n")
    list_la2 = []
    for j in list_la:
        print j
        list_la2.append(j[0])
    lw = open(readthroughfile + "_notdefined", "w+")
    lw.writelines("\n".join(list_la2) + "\n")

def splitfile(fr_g):
    dict_s={}
    for s in fr_g:
        par=s[0]
        if par in dict_s:
            dict_s[par].append(s)
        else:
            dict_s[par]=[s]
    return dict_s

def main(argv):
    t1 = time.time()
    global readthroughfile
    global genefile
    genefile = ''
    readthroughfile=''
    try:
        opts, args = getopt.getopt(argv, "ha:o:", ["afile=","odir="])
    except getopt.GetoptError:
        print'python DNAlevel_NCL_filter.py -a <genefile> -o <readthroughfile>'
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print'python DNAlevel_NCL_filter.py -a <genefile> -o <readthroughfile>'
            sys.exit()
        elif opt in ("-a", "--afile"):
            genefile = arg
        elif opt in ("-o", "--odir"):
            readthroughfile = arg
    get_time_mem()
    t2 = time.time()
    print"The time needed in DNAlevel_NCL_filter for TSEs is:", t2 - t1

if __name__ == "__main__":
    main(sys.argv[1:])