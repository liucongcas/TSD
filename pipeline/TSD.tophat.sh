ARGS=`getopt -o hvw:i:x:y:o:s:m:l:r:p: --long help,version,softwaredir:reference:read1:read2:outputdir:min_sup:min_mismatch:library_type:RSL: -- "$@"`
if [ $? != 0 ] ; then echo "Terminating..." >&2 ; exit 1 ; fi
eval set -- "$ARGS"
while true;do
    case "$1" in
        -w|--softwaredir)
	softwaredir=$2
        shift 2    ;;
	-i|--reference)
        reference=$2
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
	-m|--min_mismatch)
	min_mismatch=$2
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
            echo "sh TSD.tophat.sh -w <softwaredir> ;required \n -i <reference> ;required \n -x <read1> ;required \n -y <read2> ;required \n -o <outputdir>; default:./tophat-results \n -s <min_sup> ;default:4 \n -m <min_mismatch> ;default:1 \n -l <library_type> ;default:fr-unstranded  \n -r <RSL> ;default:25\n -p <threads> ; default:4"
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
min_mismatch=${min_mismatch:-1}
outputdir=${outputdir:-./tophat-results}
library_type=${library_type:-fr-unstranded}
threads=${threads:-4}


if [ ! $softwaredir ]; then
    software=./tools/tophat
else
    software=${softwaredir}/tools/tophat
fi

$software --bowtie1 -p $threads --b2-very-sensitive --mate-std-dev 50 --fusion-search --fusion-min-dist 500000 --fusion-anchor-length ${RSL} --fusion-read-mismatches $min_mismatch --library-type $library_type -o $outputdir $reference  $read1 $read2

