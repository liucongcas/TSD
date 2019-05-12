ARGS=`getopt -o hvw:i:x:y:o:s:m:l:r:p: --long help,version,softwaredir:index:read1:read2:outputdir:min_sup:max_mismatch:library_type:RSL: -- "$@"`
if [ $? != 0 ] ; then echo "Terminating..." >&2 ; exit 1 ; fi
eval set -- "$ARGS"
while true;do
    case "$1" in
        -w|--softwaredir)
	softwaredir=$2
        shift 2    ;;
	-i|--index)
        index=$2
        shift 2    ;;
	-x|--read1)
	read1=$2
	shift 2	;;
	-y|--read2)
        read2=$2
        shift 2 ;;
        -o|--outputdir)
	outputdir=$2
        shift 2    ;;
	-s|--min_sup)
	min_sup=$2
        shift 2    ;;
	-m|--max_mismatch)
	max_mismatch=$2
        shift 2    ;;
	-l|--library_type)
	library_type=$2
        shift 2    ;;
	-r|--RSL)
	RSL=$2
        shift 2    ;;
	-p|--threads)
	threads=$2
	shift 2    ;;
        -v|--version)
            echo "TSD-V1.0"
            ;;
        -h|--help)
		echo "sh TSD.tophat.sh \n -w <softwaredir> ;required \n -i <index> ;required \n -x <read1> ;required \n -y <read2> ;required \n -o <outputdir>; default:./tophat-results \n -s <min_sup> ;default:4 \n -m <max_mismatch> ; For RNA-seq data at high read-depth (over 20-fold), 0 is suggested. Default:1 \n -l <library_type> ;default:fr-unstranded  \n -r <RSL> ; the span length of junction read over the junction site (RSL). For RNA-seq data with read-length (less than 75bp), 15 is suggested. Default:25\n -p <threads> ; default:4"
	    exit 1 ;
            shift
            ;;
        --)
            shift
            break
            ;;
        *) 
            echo "Unknown attribute:{$1}"
            exit 1
            ;;
    esac
done


RSL=${RSL:-25}
min_sup=${min_sup:-4}
max_mismatch=${max_mismatch:-1}
outputdir=${outputdir:-./tophat-results}
library_type=${library_type:-fr-unstranded}
threads=${threads:-4}


if [ ! $softwaredir ]; then
    software=./tools/tophat
else
    software=${softwaredir}/tools/tophat
fi

echo "$software --bowtie1 -p $threads --b2-very-sensitive --mate-std-dev 50 --fusion-search --fusion-min-dist 500000 --fusion-anchor-length ${RSL} --fusion-read-mismatches $max_mismatch --library-type $library_type -o $outputdir $index  $read1 $read2"
$software --bowtie1 -p $threads --b2-very-sensitive --mate-std-dev 50 --fusion-search --fusion-min-dist 500000 --fusion-anchor-length ${RSL} --fusion-read-mismatches $max_mismatch --library-type $library_type -o $outputdir $index  $read1 $read2

