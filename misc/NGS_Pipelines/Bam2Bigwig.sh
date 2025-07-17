#!/bin/bash

#/ Convert BAM files to Bigwigs

set -e -o pipefail
LC_ALL=C

export VERSION=1.0.0

#/ Help section:
usage(){
  echo "

## Bam2Bigwig.sh Version: ${VERSION}

=> Options:
   
------------------------------------------------------------------------------------------------------------------
-h | --help          : show this message                                                           {}
-b | --bams          : space-delimited string of input bam files in double quotes or *.bam         {}
-m | --mode          : <single> or <paired>, if paired will use -pc option in bedtools genomecov   {}
-a | --atacseq       : use +4/-5 shifted 5' end of reads for pileup                                {FALSE}
-e | --extend        : numeric value to extend reads to fragment size, see details                 {0}
-n | --normalize     : if set then normalizes bedGraphs using TMM from edgeR, see details          {FALSE}
-k | --useexistingcm : use this existing count matrix to calculate SFs from                        {FALSE}
-u | --useexistingsf : use existing scaling_factors.txt to grep SFs from                           {FALSE}
-p | --peaks         : peaks in BED format                                                         {}
-j | --njobs         : number of parallel (GNU parallel) jobs to create bedGraphs, see details.    {1}
-t | --threads       : number of threads for featureCounts (if --normalize)                        {1}
-q | --sortthreads   : number of threads for sorting the bedGraph                                  {1}
-w | --sortmem       : memory for sorting the bedGraph                                             {1G}
-s | --customsuffix  : a suffix to append to the output file name , e.g. outname<_suffix>.bigwig   {}
------------------------------------------------------------------------------------------------------------------

Details: * Option --extend is simply -fs in bedtools genomecov. If --atacseq is set then this can be
           used to smooth the coverage of the cutting sites for visualization. It will then use
           bedtools slop -b to extend the cutting sites to both directions, so if a 100bp window is
           desired one should set it to 50 to have 50bp extension to each direction.
           
         * Option --normalize works together with --peaks. The peak list will be used to create a count
           matrix based on all BAM files in --bams which is then passed to edgeR::calcNormFactors() to
           calculate normalization factors (library size + composition).
           The 4th column of the bedGraphs is then divided by the per-library normalization factor.
           
         *  Option --njobs => one job requires 2-5 threads plus the --sortthreads threads so be careful.

"
}; if [[ -z "$1" ]] || [[ $1 == -h ]] || [[ $1 == --help ]]; then usage; exit; fi

#------------------------------------------------------------------------------------------------------------------------

#/ Set args:
for arg in "$@"; do                         
 shift                                      
 case "$arg" in
   "--bams")          set -- "$@" "-b"    ;;   
   "--mode")          set -- "$@" "-m"    ;;   
   "--atacseq")       set -- "$@" "-a"    ;;   
   "--extend")        set -- "$@" "-e"    ;;   
   "--normalize")     set -- "$@" "-n"    ;;   
   "--useexistingsf") set -- "$@" "-u"    ;;   
   "--useexistingcm") set -- "$@" "-k"    ;;   
   "--peaks")         set -- "$@" "-p"    ;;   
   "--njobs")         set -- "$@" "-j"    ;;   
   "--threads")       set -- "$@" "-t"    ;;   
   "--sortthreads")   set -- "$@" "-q"    ;;   
   "--sortmem")       set -- "$@" "-w"    ;;   
   "--customsuffix")  set -- "$@" "-s"    ;;
   *)                 set -- "$@" "$arg"  ;;   
 esac
done

#/ Set defaults and check sanity:
bams=""
mode=""
atacseq="FALSE"
extend="0"
normalize="FALSE"
useexistingsf="FALSE"
useexistingcm="FALSE"
peaks=""
njobs="2"
threads="1"
sortthreads="1"
sortmem="1G"
customsuffix=""

#/ getopts and export:
while getopts b:m:e:p:t:j:q:w:k:s:anu OPT           
  do   
  case ${OPT} in
    b) bams="${OPTARG}"          ;;
    m) mode="${OPTARG}"          ;;
    e) extend="${OPTARG}"        ;;
    p) peaks="${OPTARG}"         ;;
    t) threads="${OPTARG}"       ;;
    j) njobs="${OPTARG}"         ;;
    a) atacseq="TRUE"            ;;
    n) normalize="TRUE"          ;;
    u) useexistingsf="TRUE"      ;;
    k) useexistingcm="${OPTARG}" ;;
    q) sortthreads="${OPTARG}"   ;;
    w) sortmem="${OPTARG}"       ;;
    s) customsuffix="${OPTARG}"  ;;
  esac
done	

if [[ "${bams}" == "" ]]; then
  echo 'Either --bam or --bigwig must be specified'
  exit 1
fi

if [[ "${atacseq}" == "TRUE" ]]; then mode="single"; fi

if [[ "${mode}" == "" ]] && [[ "${mode}" != "single" ]] && [[ "${mode}" != "paired" ]]; then
  echo 'Argument --mode must be specified and must be <single> or <paired>'
  exit 1
fi

if [[ "${threads}" == "0" ]]; then threads=1; fi

OPTS=(bams mode extend peaks threads njobs atacseq \
     normalize sortthreads sortmem useexistingsf useexistingcm customsuffix)
for i in ${OPTS[*]}; do export $i; done

#/ Print summary:

echo ''
echo '---------------------------------------------------------------------------------------------'
echo '[Version] ::: '"${VERSION}"
echo ''
echo '[Info] Running with these parameters:'
echo '       --bams          = '"${bams}"
echo '       --mode          = '"${mode}"
echo '       --extend        = '"${extend}"
echo '       --peaks         = '"${peaks}"
echo '       --threads       = '"${threads}"
echo '       --njobs         = '"${njobs}"
echo '       --atacseq       = '"${atacseq}"
echo '       --normalize     = '"${normalize}"
echo '       --useexistingsf = '"${useexistingsf}"
echo '       --useexistingcm = '"${useexistingcm}"
echo '       --sortthreads   = '"${sortthreads}"
echo '       --sortmem       = '"${sortmem}"
echo '       --customsuffix  = '"${customsuffix}"
echo '---------------------------------------------------------------------------------------------'
echo ''

#------------------------------------------------------------------------------------------------------------------------

#/ Function that checks if required tools are callable:
function PathCheck {
  
  if [[ $(command -v $1 | wc -l | xargs) == 0 ]]; then 
  echo ${1} >> missing_tools.txt
  fi
  
}; export -f PathCheck

Tools=(bedGraphToBigWig bedtools bgzip featureCounts parallel Rscript)

#/ if that argument is set check for tools, then exit:
if [[ -e missing_tools.txt ]]; then rm missing_tools.txt; fi

for i in $(echo ${Tools[*]}); do PathCheck $i; done

if [[ $(Rscript -e '"edgeR" %in% rownames(installed.packages())' | cut -d ' ' -f2) == "FALSE" ]]; then
    echo 'R package: edgeR' >> missing_tools.txt
fi    

if [[ -e missing_tools.txt ]]; then 
  echo '[Error] Tools missing in PATH, see missing_tools.txt'
  exit 1
fi 

#------------------------------------------------------------------------------------------------------------------------

#/ A function to create a count matrix based on the provided BED file, and then use edgeR to calculate
#/ scaling factors:

if [[ "${normalize}" == "TRUE" ]] && [[ "${useexistingsf}" == "FALSE" ]]; then 

cat <<EOF > calculateTMM.R
  
#/ Rscript to calculate TMM factors based on a count matrix:
  
pkg <- c("edgeR")
  
if(sum(unlist(lapply(pkg, function(x) requireNamespace(x, quietly = TRUE)))) != length(pkg)){
  stop("Not all required packages installed or cannot be loaded!")
} else invisible(lapply(pkg, require, character.only = TRUE))
  
#/ read data from stdin:
raw.counts <- read.table(file = file('stdin'), skip = 1, header = TRUE, sep="\t")
  
#/ stop if only one (or no) samples found, featureCounts has 6 mandatory columns,
#/ and 7th to ncol are the samples
if (ncol(raw.counts) < 8) stop("Only one sample present, so skipping calculation!")

raw.counts <- raw.counts[,7:ncol(raw.counts)]
  
#/ total lib size and nf to get effective libsize per million:
if (file.exists("scaling_factors.txt")) apd <- TRUE else apd <- FALSE
nf      <- calcNormFactors(object = raw.counts, method = "TMM")
libsize <- colSums(raw.counts)
  
write.table(x = data.frame(Sample = names(libsize), ScalingFactor = nf * libsize / 10^6),
            file = "scaling_factors.txt",
            append = apd, sep = "\t", col.names = !apd, row.names = FALSE, quote = FALSE)
            
EOF
  
fi  

function GetSizeFactors {
  
  if [[ "${useexistingcm}" == "FALSE" ]]; then
  
    #/ get SAF format from peak list:
    mawk 'OFS="\t" {print $1"_"$2"_"$3, $1, $2, $3, "."}' "${peaks}" > "${peaks}".saf
  
    #/ make count matrix:
    if [[ "${atacseq}" == "TRUE" ]]; then r2p='--read2pos 5'; else r2p=''; fi
    if [[ "${mode}" == "paired" ]]; then pe='-p'; else pe=''; fi
  
    eval featureCounts -a "${peaks}".saf -o "${peaks}".counts -F SAF -T "${threads}" "${r2p}" "${pe}" "${bams}"
  
    if [[ -e scaling_factors.txt ]]; then mv scaling_factors.txt scaling_factors_old.txt; fi
      
    #/ feed to edgeR to get size factors:
    cat "${peaks}".counts | Rscript calculateTMM.R || exit 1
    
  else
  
    #/ use existing count matrix:
    cat "${useexistingcm}" | Rscript calculateTMM.R || exit 1
    
  fi  
  
}; export -f GetSizeFactors

#------------------------------------------------------------------------------------------------------------------------
#/ Function to produce BigWig from Bam:
function Bam2Bw {
  
  #/ Input is one bam file per sample:
  singlebam="${1}"
  
  if [[ ! -e "${singlebam}" ]]; then 
    echo '[ERROR]' "${singlebam}" 'does not exist'
    exit 1
  fi  
  
  Basename=${singlebam%.bam}
  
  #/ Get chromsizes based on idxstats
  samtools view -H "${singlebam}" | grep 'SN:' | cut -f2,3 | awk '{gsub("SN:|LN:", "");print}' > "${singlebam}".chromsize
  
  #---------------------------------------------
  
  #/ grep the per-sample scaling factor:
  if [[ "${normalize}" == "TRUE" ]]; then
  
    #/ take the reciprocal since the scale factor is intended for division but genomecov multiplies:
    SF=$(bc <<< "scale=10; $(grep -w ${singlebam} scaling_factors.txt | cut -f2)^-1")
    
    if [[ "${SF}" == "" ]]; then
      echo "[Error]" "${singlebam}" "not found in scaling_factors.txt"
      exit 1
    fi
    
  else SF=1
  fi
  
  #---------------------------------------------
  
  #/ function to get cutsites from bam2bed output to stdout:
  getCutsite() { 
  
    if (( $(head -n1 $1 | awk -F'\t' 'NR==1{print NF}') < 6 ));
      then
      echo '[Error]' $1 "has fewer than six columns!"
      exit 1
    fi
    
    mawk 'OFS="\t" {if ($6 == "+") print $1,$2+4,$2+5,".",".",$6} 
                   {if ($6 == "-") print $1,$3-5,$3-4,".",".",$6}' $1
                   
  }
  
  #---------------------------------------------
  
  #/ chain the commands into a generic call:
  if [[ "${atacseq}" == "TRUE" ]]; then
  
    do_bamtobed='bedtools bamtobed -i "${singlebam}" | getCutsite /dev/stdin |'
    
    if [[ "${extend}" != 0 ]]; then 
    
      do_slop='bedtools slop -g "${singlebam}".chromsize -b "${extend}" -i - |'
      suffix='_extended'
      
    else do_slop=''; suffix=''
      
    fi  
    
    input='-i -'
    csize='-g "${singlebam}".chromsize'
    is_paired=''
    
  else
  
    do_bamtobed=''
    do_slop=''
    input='-ibam "${singlebam}"'
    csize=''
    
    if [[ "${extend}" != 0 ]]; then 
      do_extend='-fs "${extend}"'
      suffix='_extended'
    else 
      do_extend=''
      suffix=''
    fi  
    
    if [[ "${mode}" == "paired" ]]; then is_paired='-pc'; else is_paired=''; fi
    
  fi

  do_scale='-scale "${SF}"'
  
  #/ put it together:
  eval \
  "${do_bamtobed}" \
  "${do_slop}" \
  bedtools genomecov -bga "${input}" "${csize}" "${do_extend}" "${is_paired}" "${do_scale}" \
  | sort -k1,1 -k2,2n --parallel="${sortthreads}" -S "${sortmem}" > "${Basename}""${suffix}".bedGraph
  
  bedGraphToBigWig "${Basename}""${suffix}".bedGraph "${singlebam}".chromsize "${Basename}""${suffix}""${customsuffix}".bigwig \
  && rm "${Basename}""${suffix}".bedGraph
  
}; export -f Bam2Bw

#------------------------------------------------------------------------------------------------------------------------

#/ Run it:
if [[ $(echo "${bams}" | tr " " "\n" | awk NF | wc -l) == 1 && "${normalize}" == "TRUE" ]]; then 

  echo 'Only one sample -- skipping normalization part'
  normalize="FALSE"
  
fi

#/ Run normalization:
if [[ "${normalize}" == "TRUE" ]]; then GetSizeFactors; fi

#/ Check whether scaling_factors.txt exists:
if [[ "${useexistingsf}" == "TRUE" ]]; then
    if [[ ! -f "scaling_factors.txt" ]]; then
        echo '[Error] --useexistingsf is set but scaling_factors.txt does not exist'
        exit 1
    fi
fi 

#/ Run the main function:
echo ${bams} | tr " " "\n" | awk NF | sort -u | parallel -j "${njobs}" "Bam2Bw {}"