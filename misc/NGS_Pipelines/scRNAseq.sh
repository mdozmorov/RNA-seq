#!/bin/bash

LC_ALL=C

usage(){
  echo '

------------------------------------------------------------------------------------------------------------
[Deprecated] Use our nextflow workflow: https://github.com/ATpoint/sc_preprocess
------------------------------------------------------------------------------------------------------------

Quantification of droplet-based scRNA_seq data using selective alignment via Alevin.

Assumes => CB/UMI in Read1 (Basename_1.fastq.gz) and
        => cDNA   in Read2 (Basename_2.fastq.gz)
        
------------------------------------------------------------------------------------------------------------
-h | --help       : Show this message                                                                    
-i | --idx        : the transcriptome index                              {}                                                                                                    
-c | --chemistry  : the kit chemistry (dropseq, chromium, chromiumV3)    {chromiumV3}
-g | --tgmap      : transcript2gene conversion table                     {}                                
-m | --mtgenes    : list with mitochondrial genes for CB whitelisting    {}                                    
-r | --rrnagenes  : list with rRNA genes for CB whitelisting             {}                                
-l | --libtype    : library type                                         {ISR}                             
-t | --threads    : number of threads per sample                         {36}
-j | --njobs      : number of parallel jobs                              {2}
-a | --additional : any additional parameters to pass to alevin          {}
------------------------------------------------------------------------------------------------------------

'
  
}; if [[ -z "$1" ]] || [[ $1 == -h ]] || [[ $1 == --help ]]; then usage; exit; fi

#-----------------------------------------------------------------------------------------------------------------------------------------

for arg in "$@"; do                         
shift                                      
case "$arg" in
"--idx")         set -- "$@" "-i"   ;;   # 
"--chemistry")   set -- "$@" "-c"   ;;   # 
"--tgmap")       set -- "$@" "-g"   ;;   #
"--mtgenes")     set -- "$@" "-m"   ;;   #
"--rrnagenes")   set -- "$@" "-r"   ;;   #
"--libtype")     set -- "$@" "-l"   ;;   #
"--threads")     set -- "$@" "-t"   ;;   #
"--njobs")       set -- "$@" "-j"   ;;   #
"--additional")  set -- "$@" "-a"   ;;   #
*)               set -- "$@" "$arg" ;;   # 
esac
done

#/ defaults:
idx=""
chemistry="chromiumV3"
tgmap=""
mtgenes=""
rrnagenes=""
libtype="ISR"
threads="36"
njobs="2"
additional=""

while getopts i:c:g:m:r:l:t:j:a: OPT      
do   
case ${OPT} in                      
i) idx="${OPTARG}"        ;;
c) chemistry="${OPTARG}"  ;;
g) tgmap="${OPTARG}"      ;;
m) mtgenes="${OPTARG}"    ;;
r) rrnagenes="${OPTARG}"  ;;
l) libtype="${OPTARG}"    ;;
t) threads="${OPTARG}"    ;;
j) njobs="${OPTARG}"      ;;
a) additional="${OPTARG}" ;;
esac
done	

#-----------------------------------------------------------------------------------------------------------------------------------------

## Check for args, set defaults if empty:
function Exit1 {
  echo '[ERROR]:' $1 'is empty or does not exist' && exit 1
}; export -f Exit1

#/ Check chemistry:
Chem=(dropseq chromium chromiumV3)
if [[ ! " ${Chem[@]} " =~ " ${chemistry} " ]]; then 
echo '[ERROR]: --chemistry must be one of {dropseq, chromium, chromiumV3}'
exit 1
fi

#/ export options for use within GNU parallel:
OPTIONS=(idx chemistry tgmap mtgenes rrnagenes libtype threads njobs additional)
for i in $(echo ${OPTIONS[*]}); do
export $i
done  

#-----------------------------------------------------------------------------------------------------------------------------------------

#/ Return summary:

echo ''
echo '----------------------------------------------------------------------------------------------------------------------------------'
echo '[Info] Running with these parameters:'
echo '       --idx        = '"${idx}"
echo '       --chemistry  = '"${chemistry}"
echo '       --tgmap      = '"${tgmap}"
echo '       --mtgenes    = '"${mtgenes}"
echo '       --rrnagenes  = '"${rrnagenes}"
echo '       --libtype    = '"${libtype}"
echo '       --threads    = '"${threads}"
echo '       --jobs       = '"${njobs}"
echo '       --additional = '"${additional}"
echo '----------------------------------------------------------------------------------------------------------------------------------'
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
TOOLS=(salmon)

#/ Loop through the tools and write to missing_tools.txt all those not in PATH / callable:
for i in $(echo ${TOOLS[*]}); do
PathCheck $i; done

#/ If tools are missing write them to <missing_tools.txt>
if [[ -e missing_tools.txt ]] && [[ $(cat missing_tools.txt | wc -l | xargs) > 0 ]]; then
echo '[ERROR] Tools missing in PATH -- see {missing_tools.txt}' && exit 1; fi

#-----------------------------------------------------------------------------------------------------------------------------------------

function ALEVIN {
  
  chem=--"${chemistry}"
  
  AlevinBasic='salmon alevin "${chem}" -l "${libtype}" -i "${idx}" -p "${threads}" \
               --tgMap "${tgmap}" --mrna "${mtgenes}" --rrna "${rrnagenes}" -o ${1}'
  
  (>&2 echo '')
  (>&2 paste -d " " <(echo '[INFO]' 'Running alevin for' "${1}" 'in' ${Mode} 'mode -- started on') <(date))
  
  eval "${AlevinBasic}" "${additional}" -1 "${1}"_1.fastq.gz -2 "${1}"_2.fastq.gz
  
  (>&2 paste -d " " <(echo '[INFO]' 'Running alevin for' "${1}" 'in' ${Mode} 'mode -- ended on') <(date))
  (>&2 echo '')
  
}; export -f ALEVIN

## Run it on all fastq files in $(pwd):
if [[ $(ls *_1.fastq.gz 2> /dev/null | awk NF | wc -l) == 0 ]]; then
echo '[ERROR] No input files found'; exit 1
fi  

ls *_1.fastq.gz \
| awk -F "_1.fastq.gz" '{print $1}' \
| parallel -j ${njobs} "ALEVIN {} 2> {}_alevin.log"
