#!/bin/bash

#/ RNA-seq quantification with salmon

set -e -o pipefail
LC_ALL=C

export VERSION=1.0.0

usage(){
  
echo "

Version: ${VERSION}

=> RNA-seq quantification with salmon.

=> Naming conventions for input fastq files are:
   
   1) single-end fastq : Basename.fastq.gz
   2) paired-end fastq : Basename_1/2.fastq.gz

=> Options:

---------------------------------------------------------------------------------------------
-h | --help        : Show this message                               {}
-i | --idx         : the transcriptome index folder                  {}
-m | --mode        : single or paired-end data (single,paired)       {}
-n | --noLength    : turn off length correction for end-tagged libs  {FALSE}
-t | --threads     : number of threads per run                       {16}
-j | --njobs       : number of parallel jobs for salmon              {4}
-l | --libtype     : library type                                    {A}
-s | --fldMean     : mean insert size for single-end data            {250}    
-q | --fldSD       : standard deviation for --fldMean                {25}
-a | --additional  : any additional salmon arguments                 {}    
-c | --trim        : whether to trim adapters via cutadapt           {FALSE}
-d | --adapter     : the adapter sequence to trim, default is TruSeq {AGATCGGAAGAGC}
-y | --trimthreads : threads per job for cutadapt                    {2}
-x | --trimjobs    : GNU parallel jobs for cutadapt                  {10}
---------------------------------------------------------------------------------------------

"

}; if [[ -z "$1" ]] || [[ $1 == -h ]] || [[ $1 == --help ]]; then usage; exit 0; fi
	
#-----------------------------------------------------------------------------------------------------------------------------------------

for arg in "$@"; do                         
 shift                                      
 case "$arg" in
   "--idx")         set -- "$@" "-i"   ;;    
   "--mode")        set -- "$@" "-m"   ;;    
   "--noLength")    set -- "$@" "-n"   ;;   
   "--threads")     set -- "$@" "-t"   ;;   
   "--njobs")       set -- "$@" "-j"   ;;   
   "--libtype")     set -- "$@" "-l"   ;;   
   "--fldMean")     set -- "$@" "-s"   ;;   
   "--fldSD")       set -- "$@" "-q"   ;;   
   "--additional")  set -- "$@" "-a"   ;; 
   "--trim")        set -- "$@" "-c"   ;;
   "--adapter")     set -- "$@" "-d"   ;;
   "--trimthreads") set -- "$@" "-y"   ;; 
   "--trimjobs")    set -- "$@" "-x"   ;; 
   *)               set -- "$@" "$arg" ;;   
 esac
done
	
#/ Defaults:
idx=""
mode=""
noLength="FALSE"
threads="16"
njobs="4"
libtype="A"
fldMean="250"
fldSD="25"
additional=""
trim="FALSE"
adapter="AGATCGGAAGAGC"
trimthreads="2"
trimjobs="10"

while getopts i:m:t:j:l:s:q:a:d:y:x:nc OPT     
  do   
  case ${OPT} in                         
    i) idx="${OPTARG}"           ;;
    m) mode="${OPTARG}"          ;;                                     
    n) noLength="TRUE"           ;;     
    t) threads="${OPTARG}"       ;;
    j) njobs="${OPTARG}"         ;;
    l) libtype="${OPTARG}"       ;;
    s) fldMean="${OPTARG}"       ;;
    q) fldSD="${OPTARG}"         ;;
    a) additional="${OPTARG}"    ;;
    c) trim="TRUE"               ;;
    d) adapter="${OPTARG}"       ;;
    y) trimthreads="${OPTARG}"   ;;
    x) trimjobs="${OPTARG}"      ;;
  esac
done	

OPTIONS=(idx mode noLength threads njobs libtype fldMean fldSD additional trim adapter \
         runfastqc trimthreads trimjobs)
         
for i in $(echo ${OPTIONS[*]}); do export $i; done  

#-----------------------------------------------------------------------------------------------------------------------------------------

#/ Return summary:

echo ''
echo '---------------------------------------------------------------------------------------------'
echo '[Version] '"${VERSION}"
echo ''
echo '[Info] Running with these parameters:'
echo '       --idx         = '"${idx}"
echo '       --mode        = '"${mode}"
echo '       --noLength    = '"${noLength}"
echo '       --jobs        = '"${njobs}"
echo '       --threads     = '"${threads}"
echo '       --libtype     = '"${libtype}"
echo '       --fldMean     = '"${fldMean}"
echo '       --fldSD       = '"${fldSD}"
echo '       --additional  = '"${additional}"
echo '       --trim        = '"${trim}"
echo '       --adapter     = '"${adapter}"
echo '       --trimthreads = '"${trimthreads}"
echo '       --trimjobs    = '"${trimjobs}"
echo '---------------------------------------------------------------------------------------------'
echo ''

#-----------------------------------------------------------------------------------------------------------------------------------------
	
#/ Check if required tools are in PATH and/or callable:

if [[ -e missing_tools.txt ]]; then rm missing_tools.txt; fi

#/ Function that checks if required tools are callable:
function PathCheck {
  
  if [[ $(command -v $1 | wc -l | xargs) == 0 ]]; then 
    echo ${1} >> missing_tools.txt
  fi
    
}; export -f PathCheck

#/ All tools:
TOOLS=(salmon cutadapt)

#/ Loop through the tools and write to missing_tools.txt all those not in PATH / callable:
for i in $(echo ${TOOLS[*]}); do
  PathCheck $i; done
  
#/ If tools are missing write them to <missing_tools.txt>
if [[ -e missing_tools.txt ]] && [[ $(cat missing_tools.txt | wc -l | xargs) > 0 ]]; then
  echo '[Error] Tools missing in PATH -- see missing_tools.txt' && exit 1; fi

#-----------------------------------------------------------------------------------------------------------------------------------------

function TRIM {
  
  Basename="${1}"
  
  if [[ ! -e trimdir ]]; then mkdir trimdir; fi
  if [[ ! -e untrimdir ]]; then mkdir untrimdir; fi
  
  CUT='cutadapt --quiet -j "${trimthreads}" -m 31 --max-n 0.1 -a "${adapter}"'
  
  if [[ "${mode}" == "single" ]]; then 
    
    CMD_CUT=$(echo "${CUT}" -o trimdir/"${Basename}"".fastq.gz")
    
    eval "${CMD_CUT}" "${Basename}"".fastq.gz" && \
    mv ${Basename}*.gz ./untrimdir/ && \
    mv ./trimdir/${Basename}*.gz .
    
  fi  
  
  if [[ ${mode} == "paired" ]]; then 
  
    CMD_CUT=$(echo "${CUT}" -A "${adapter}" -o trimdir/"${Basename}"_1.fastq.gz -p trimdir/"${Basename}"_2.fastq.gz)
    
    eval "${CMD_CUT}" "${Basename}""_1.fastq.gz" "${Basename}""_2.fastq.gz" && \
    mv ${Basename}*.gz ./untrimdir/ && \
    mv ./trimdir/${Basename}*.gz .
    
  fi  
  
}; export -f TRIM

#-----------------------------------------------------------------------------------------------------------------------------------------

function SALMON {
  
  (>&2 echo '')
  (>&2 paste -d " " <(echo '[INFO]' 'Running salmon for' "${1}" 'in' ${mode} 'mode -- started on') <(date))
  
  #/ Basic command that is always executed:
  SalmonBasic="salmon quant \
                 -l ${libtype} -i ${idx} -p ${threads} --validateMappings --no-version-check"
  
  #---------------------------------------------------------------------------------------------------------------------------------------
  
  #/ Settings for paired-end mode:
  if [[ ${mode} == "paired" ]]; then 
    
    if [[ ! -e ${1}_1.fastq.gz || ! -e ${1}_2.fastq.gz ]]; then
      echo '[Error]: At least on of the input files is missing for' $1 && exit 1
    fi
    
    Files="-1 ${1}_1.fastq.gz -2 ${1}_2.fastq.gz"
    BiasFlags="--gcBias --seqBias"
    
    #/ If data are end-tagged then set the option, if not add seqBias to BiasFlags:
    if [[ "${noLength}" == "TRUE" ]]; then
      noLength="--noLengthCorrection"
      BiasFlags=""
    else noLength=""
    fi
    
  fi
  
  #/ Settings for single-end:
  if [[ ${mode} == "single" ]]; then 
    
    if [[ ! -e ${1}.fastq.gz ]]; then
      echo '[Error]: Input files is missing for' $1 && exit 1
    fi
    
    Files="-r ${1}.fastq.gz"
    BiasFlags="--seqBias"
    
    #/ If data are end-tagged then set the option, if not add seqBias to BiasFlags:
    if [[ "${NoLength}" == "TRUE" ]]; then
      NoLength="--noLengthCorrection"
      BiasFlags=""
    else noLength=""
    fi
    
  fi
  
  #---------------------------------------------------------------------------------------------------------------------------------------
  
  #/ Run it:
  eval "${SalmonBasic}" "${additional}" "${BiasFlags}" "${noLength}" -o ${1}_salmon "${Files}"
  
  #---------------------------------------------------------------------------------------------------------------------------------------
  
  (>&2 echo '')
  (>&2 paste -d " " <(echo '[INFO]' 'Running salmon for' "${1}" 'in' ${mode} 'mode -- ended on') <(date))
 
}; export -f SALMON

#-----------------------------------------------------------------------------------------------------------------------------------------

#/ Sanity checks:
if [[ ! -d "${idx}" ]]; then 
  echo '[ERROR] --idx directory does not seem to exist'; exit 1
fi  

array=(paired single)
if [[ ! "${array[@]}" =~ "${mode}" ]]; then
  echo "[Error] -m/--mode must be <single> or <paired> -- exiting" 
  exit 1
fi

#/ expected suffix:
if [[ ${mode} == "paired" ]]; then Suffix="_1.fastq.gz"; fi
if [[ ${mode} == "single" ]]; then Suffix=".fastq.gz"; fi

#/ optional trimming:
if [[ ${trim} == "TRUE" ]]; then
  
  ls *"${Suffix}" \
  | awk -F "${Suffix}" '{print $1}' \
  | parallel -j "${trimjobs}" "TRIM {} 2>> {}.log"

fi

#/ run salmon:
ls *"${Suffix}" \
  | awk -F "${Suffix}" '{print $1}' \
  | parallel -j ${njobs} "SALMON {} 2>> {}.log"
