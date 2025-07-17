#!/bin/bash

#/ DNA-seq processing pipeline

set -e -o pipefail
LC_ALL=C

export VERSION=1.0.0------------------------------------------------------------------------------------------------------------------------

#/ Display help section:
usage(){

  echo "

## DNAseq.sh Version: ${VERSION}
  
   Script performs alignment, filtering and basic QC for assays such as ATAC-seq and ChIP-seq.
   It runs on all input files in the working directory from where it is launched.
   
   => The expected naming conventions for input files (either fastq.gz or uBAM) are:
   1) single-end fastq : Basename.fastq.gz
   2) paired-end fastq : Basename_1/2.fastq.gz
   3)            ubam  : Basename_ubam.bam

   
=> Options:

------------------------------------------------------------------------------------------------------------------
-h | --help          : show this message                                                           {}
-g | --genome        : path to genome index prefix                                                 {}
-j | --jobnumber     : number of parallel jobs (GNU parallel)                                      {4}
-t | --threads       : number of threads per job for alignment                                     {16}          
-f | --format        : type of input files (fq_se/pe, bam_se/pe)                                   {}           
-a | --atacseq       : turn on ATAC-seq mode                                                       {FALSE}
-o | --onlyAln       : perform only the alignment and then exit                                    {FALSE}
-m | --chrM          : name of the mitochondrial chromosome (for ATAC-seq % reads mapping to it)   {chrM}
-l | --memSort       : memory per thread for BAM sorting by samtools in the form of e.g. "1G"      {1G}
-y | --threadsSort   : threads for sorting operations                                              {8}
-c | --chrRegex      : use this regex to select chromosomes to keep during filtering               {chr[1,9,X,Y]}
-k | --keepDup       : do not remove PCR duplicates from alignment                                 {FALSE}  
-q | --minMAPQ       : keep only alignments with that minimal MAPQ                                 {20}
-x | --checktools    : check whether required software is in PATH                                  {}
-d | --noalignment   : do not run alignment                                                        {FALSE}
-s | --nofrips       : do not run the QC assessment via peak calling and FRiP calculation          {FALSE}
-b | --layout        : whether se or pe, only relevant if --noalignment is set (se,pe)             {se}
-w | --fripqcJobs    : number of parallel jobs for the QC/FRiP part                                {36}
-u | --gflag         : for macs2 flag -g (hs for human, mm for mouse)                              {mm}
------------------------------------------------------------------------------------------------------------------
"
  
}; if [[ $# -eq 0 ]] || [[ $1 == -h ]] || [[ $1 == --help ]]; then usage; exit 0; fi
           
#------------------------------------------------------------------------------------------------------------------------
         
#/ Parse arguments
for arg in "$@"; do                         
  shift                                      
  case "$arg" in
    "--genome")      set -- "$@" "-g"   ;;   # 
    "--jobnumber")   set -- "$@" "-j"   ;;   #
    "--threads")     set -- "$@" "-t"   ;;   #
    "--format")      set -- "$@" "-f"   ;;   #
    "--atacseq")     set -- "$@" "-a"   ;;   #
    "--onlyaln")     set -- "$@" "-o"   ;;   #
    "--chrM")        set -- "$@" "-m"   ;;   #
    "--memSort")     set -- "$@" "-l"   ;;   #
    "--threadsSort") set -- "$@" "-y"   ;;   #
    "--chrRegex")    set -- "$@" "-c"   ;;   #
    "--keepDup")     set -- "$@" "-k"   ;;   #
    "--minMAPQ")     set -- "$@" "-q"   ;;   #
    "--checktools")  set -- "$@" "-x"   ;;   #
    "--noalignment") set -- "$@" "-d"   ;;   #
    "--nofrips")     set -- "$@" "-s"   ;;   #
    "--layout")      set -- "$@" "-b"   ;;   #
    "--fripqcJobs")  set -- "$@" "-w"   ;;   #
    "--gflag")       set -- "$@" "-u"   ;;   #
    *)               set -- "$@" "$arg" ;;   # 
  esac
done

#/ Default opts:
jobnumber=4
threads=16
atacseq="FALSE"
onlyaln="FALSE"
chrM="chrM"
memSort="1G"
threadsSort=8
chrRegex="chr[1-9,X,Y]"
keepDup=FALSE
minMAPQ=20
checktools="FALSE"
runfastqc="FALSE"
noalignment="FALSE"
nofrips="FALSE"
layout="se"
fripqcJobs=16
gflag="mm"


while getopts g:j:t:f:m:l:y:c:q:b:w:u:aokxds OPT           
  do   
  case ${OPT} in
    g) genome="${OPTARG}"      ;;
    j) jobnumber="${OPTARG}"   ;;             
    t) threads="${OPTARG}"     ;;
    f) format="${OPTARG}"      ;;
    a) atacseq="TRUE"          ;;
    o) onlyaln="TRUE"          ;;
    m) chrM="${OPTARG}"        ;;
    l) memSort="${OPTARG}"     ;;
    y) threadsSort="${OPTARG}" ;;
    c) chrRegex="${OPTARG}"    ;;
    k) keepDup="TRUE"          ;;
    q) minMAPQ="${OPTARG}"     ;;
    x) checktools="TRUE"       ;;
    d) noalignment="TRUE"      ;;
    s) nofrips="TRUE"          ;;
    b) layout="${OPTARG}"      ;;
    w) fripqcJobs="${OPTARG}"  ;;
    u) gflag="${OPTARG}"       ;;
  esac
done	

#/ Export all opts so they can be used inside functions:
OPTS=(genome jobnumber threads format atacseq onlyaln chrM \
      memSort threadsSort chrRegex keepDup minMAPQ checktools \
      noalignment nofrips layout fripqcJobs gflag)

for i in ${OPTS[*]}; do export $i; done

#------------------------------------------------------------------------------------------------------------------------

#/ Function that checks if required tools are callable:
function PathCheck {
  
  if [[ $(command -v $1 | wc -l | xargs) == 0 ]]; then 
  echo ${1} >> missing_tools.txt
  fi
  
}; export -f PathCheck

#/ List of required tools:

TOOLS=(bedGraphToBigWig bgzip bowtie2 cutadapt fastqc featureCounts macs2 mawk \
       parallel picard pigz samblaster samtools seqtk)

#/ if that argument is set check for tools, then exit:
if [[ -e missing_tools.txt ]]; then rm missing_tools.txt; fi

for i in $(echo ${TOOLS[*]}); do PathCheck $i; done

if [[ ${checktools} == "TRUE" ]] && [[ -e missing_tools.txt ]]; then 
  echo 'Missing tools listed in missing_tools.txt'
  exit 0
fi

if [[ ${checktools} == "TRUE" ]] && [[ ! -e missing_tools.txt ]]; then 
  echo 'All required tools found in PATH'
  exit 0
fi  

if [[ ${checktools} == "FALSE" ]] && [[ -e missing_tools.txt ]]; then 
  echo '[ERROR] Tools missing in PATH, see missing_tools.txt'
  exit 1
fi  

if [[ ${noalignment} == "FALSE" ]]; then

  if [[ "${genome}" == "" ]]; then echo '[ERROR] -g/--genome most be provided'; exit 1; fi

  if [[ "${genome}" != "" ]]; then 
    if (( $(ls "${genome}".1.bt2 2> /dev/null | awk NF | wc -l) == 0 )); then
      echo '[ERROR] -g/--genome path does not seem to direct to the index files!'
      exit 1
    fi
  fi  
fi  

if [[ $format == "" ]]; then 
  echo '[ERROR] -f/--format most be provided'; exit 1
fi  

if [[ ${checktools} == "FALSE" ]]; then

  echo ''
  echo '---------------------------------------------------------------------------------------------'
  echo '[Info] Running with these parameters:'
  echo '       -g/--genome       =' "${genome}"
  echo '       -j/--jobnumber    =' "${jobnumber}"
  echo '       -t/--threads      =' "${threads}"
  echo '       -f/--format       =' "${format}"  
  echo '       -a/--atacseq      =' "${atacseq}"
  echo '       -o/--onlyaln      =' "${onlyaln}"
  echo '       -m/--chrM         =' "${chrM}"
  echo '       -l/--memSort      =' "${memSort}"
  echo '       -y/--threadsSort  =' "${threadsSort}"
  echo '       -c/--chrRegex     =' "${chrRegex}"
  echo '       -k/--keepdup      =' "${keepDup}"
  echo '       -q/--minMAPQ      =' "${minMAPQ}"
  echo '       -d/--noalignment  =' "${noalignment}"
  echo '       -s/--nofrips      =' "${nofrips}"
  echo '       -b/--layout       =' "${layout}"
  echo '       -w/--fripqcJobs   =' "${fripqcJobs}"
  echo '       -u/--gflag        =' "${gflag}"
  echo '---------------------------------------------------------------------------------------------'
  echo ''

fi  

#------------------------------------------------------------------------------------------------------------------------

#/ Exit function if BAM file looks corrupted or is missing after a step:
function ExitBam {
  
  (>&2 echo '[ERROR]' "${1}" 'looks suspicious or is empty -- exiting') && exit 1
  
}; export -f ExitBam

#/ Check if BAM is corrupted or empty:
function BamCheck {
  
  Basename="${1%_*}"
  samtools quickcheck -q $1 && echo '' >/dev/null || ExitBam "${1}"
  
  ## Also check if file is not empty:
  if [[ $(samtools view "${1}" | head -n 1 | wc -l) < 1 ]]; then
    ExitBam $Basename
  fi  
  
}; export -f BamCheck  

#------------------------------------------------------------------------------------------------------------------------

#/ Function to get the % of reads mapped to chrM:
function mtDNA {
  
  Bam=$1
  
  mtReads=$(samtools idxstats ${Bam} | awk -v var="${chrM}" '{if($1==var) print $3}')
  totalReads=$(samtools idxstats ${Bam} | awk '{SUM += $3} END {print SUM}')
  
  echo $(bc <<< "scale=2;100*$mtReads/$totalReads")
  
}; export -f mtDNA

#------------------------------------------------------------------------------------------------------------------------

#/ Main alignment /filtering function:

function Fq2Bam {
  
  Basename="${1}"
  
  (>&2 echo '')
  (>&2 paste -d " " <(echo '[INFO]' 'Fq2Bam for' "${Basename}" 'in' ${SePe} 'mode started on') <(date))
  
  #/ single or paired:
  SePe=$(echo "${format}" | tr "_" "\t" < /dev/stdin | cut -f2)
  
  #(>&2 paste -d " " <(echo '[INFO]' 'Pipeline for' "${Basename}" 'in' ${SePe} 'mode started on') <(date))
  
  #/ default adapter types (TruSeq or Nextera)
  if [[ ${atacseq} == "TRUE" ]]; then ADAPTER="CTGTCTCTTATACACATCT"; else ADAPTER="AGATCGGAAGAGC"; fi  
  
  #----------------------------------------------------------------------------------------------------------------------
  
  #/ the basic cutadapt and bowtie2 commands:
  CUT='cutadapt --quiet -j 1 -m 18 --max-n 0.1 -a "${ADAPTER}"'
  B2FQ='samtools fastq -n "${Basename}"_ubam.bam'
  MERGEPE='seqtk mergepe "${Basename}"_1.fastq.gz "${Basename}"_2.fastq.gz'
  BOWTIE2='bowtie2 -q --very-sensitive --threads "${threads}" --rg-id "${Basename}" -x "${genome}"'
  B_SE='-U -'
  B_PE='-X 2000 --interleaved -'
  
  #/ if paired:
  if [[ ${format} == "fq_se" ]];  then CMD_CUT=$(echo "${CUT}" "${Basename}.fastq.gz"); fi
  if [[ ${format} == "fq_pe" ]];  then CMD_CUT=$(echo "${MERGEPE}" "|" "${CUT}" '-A "${ADAPTER}" --interleaved -'); fi
  if [[ ${format} == "bam_se" ]]; then CMD_CUT=$(echo "${B2FQ}" "|" "${CUT}" '-'); fi
  if [[ ${format} == "bam_pe" ]]; then CMD_CUT=$(echo "${B2FQ}" "|" "${CUT}" '-A "${ADAPTER}" --interleaved -'); fi
  
  #/ concat the commands, output would be SAM file to stdout:
  if [[ ${SePe} == "se" ]]; then CMD_RUN=$(echo "${CMD_CUT}" "|" "${BOWTIE2}" "${B_SE}"); fi
  if [[ ${SePe} == "pe" ]]; then CMD_RUN=$(echo "${CMD_CUT}" "|" "${BOWTIE2}" "${B_PE}"); fi
  
  #/ Alignment, dup marking and sorting plus index on-the-fly (samtools v1.11, force bai index)
  eval "${CMD_RUN}" \
  | samblaster --quiet --ignoreUnmated \
  | samtools sort 2> /dev/null -m "${memSort}" -T "${Basename}" \
    -@ "${threadsSort}" -l 5 -O BAM \
  | tee "${Basename}"_raw.bam \
  | samtools index - "${Basename}"_raw.bam.bai
  
  BamCheck "${Basename}"_raw.bam
  
  #/ percent mt:
  if [[ ${atacseq} == "TRUE" ]]; then 
    paste <(echo $Basename)  <(mtDNA "${Basename}"_raw.bam) >> mtDNA_percent.txt
  fi  
  
  #/ chromosome sizes
  if [[ ! -e tmp_chromSizes.txt ]]; then
    samtools idxstats "${Basename}"_raw.bam | grep -v '*' > tmp_chromSizes.txt
    fi
  
  ## Decide appropriate flag to keep only mapped reads depending on paired/single-end sequencing:
  if [[ $SePe == "se" ]]; then FlagKeep=0; fi
  if [[ $SePe == "pe" ]]; then FlagKeep=1; fi  
  
	#/ Now define filtering that happens for all DNA-seq:
	if [[ $keepDup == "TRUE" ]]; then FlagRemove=2308; else FlagRemove=3332; fi
	if [[ $chrRegex == "" ]]; then grepChr='.'; else grepChr="${chrRegex}";	fi 
	
  Filter='samtools idxstats "${Basename}"_raw.bam \
          | cut -f1 \
          | grep "${grepChr}" | grep -v "*" \
          | xargs \
          samtools view -O bam,level=6 -@ 2 -q "${minMAPQ}" -f "${FlagKeep}" -F "${FlagRemove}" "${Basename}"_raw.bam'
                   
  eval "${Filter}" \
  | tee "${Basename}"_dedup.bam \
  | samtools index - "${Basename}"_dedup.bam.bai
  
  BamCheck "${Basename}"_dedup.bam                 
  ls "${Basename}"_*.bam | parallel "samtools flagstat -@ 2 {} > {.}.flagstat"
  
  if [[ $onlyAln == "TRUE" ]]; then
    echo '--onlyAln is set, therefore now exiting'
    exit 0
  fi
  
  #/ If ATAC-seq mode produce bigwig of cutting sites:
  if [[ ${atacseq} == "TRUE" ]]; then
  
    bedtools bamtobed -i "${Basename}"_dedup.bam \
    | mawk 'OFS="\t" {if ($6 == "+") print $1, $2+4, $2+5, ".", ".", $6} {if ($6 == "-") print $1, $3-5, $3-4, ".", ".", $6}' \
    | sort -k1,1 -k2,2n -k3,3n -k6,6 -S"${memSort}" --parallel="${threadsSort}" \
    | tee >(bgzip -@ 6 > "${Basename}"_cutsites.bed.gz) \
    | bedtools genomecov -bg -i - -g tmp_chromSizes.txt > "${Basename}"_cutsites.bedGraph
    bedGraphToBigWig "${Basename}"_cutsites.bedGraph tmp_chromSizes.txt "${Basename}"_cutsites.bigwig && \
    rm "${Basename}"_cutsites.bedGraph
  
  fi
  
  #/ if paired-end then collect TLENs:
  if [[ $SePe == "pe" ]]; then
  
    picard CollectInsertSizeMetrics \
    --INPUT "${Basename}"_dedup.bam \
    --OUTPUT "${Basename}"_InsertSizes.txt \
    --Histogram_FILE "${Basename}"_InsertSizes.pdf \
    --QUIET true --VERBOSITY ERROR --VALIDATION_STRINGENCY LENIENT 2> /dev/null
    
  fi
  
  (>&2 paste -d " " <(echo '[INFO]' 'Fq2Bam for' "${Basename}" 'in' ${SePe} 'mode ended on') <(date))
  (>&2 echo '')
  
}; export -f Fq2Bam

#------------------------------------------------------------------------------------------------------------------------

## QC via peak calling and FRiP calculation:
function Peaks {
  
  Basename="${1}"
  
  (>&2 echo '')
  (>&2 paste -d " " <(echo '[INFO] Peaks/FRiPs for' $Basename 'started on') <(date))
  
  if [[ "${noalignment}" == "TRUE" ]]; then
    SePe2="${layout}"
  else SePe2=$(echo "${format}" | tr "_" "\t" < /dev/stdin | cut -f2)
  fi
  
  #/ check for input files, if present call peaks:
  if [[ "${atacseq}" == "TRUE" ]]; then
  
    if [[ $(ls "${Basename}"_cutsites.bed.gz 2> /dev/null | awk NF | wc -l | xargs) == 0 ]]; then
      echo "[Error]: No input files found for peak calling"
      exit 1
    fi
    
    macs2 callpeak \
      --tempdir ./ \
      -t "${Basename}"_cutsites.bed.gz \
      -n "${Basename}" -g ${gflag} \
      --extsize 100 --shift -50 --nomodel --keep-dup=all \
      -f BED -q 0.005 --min-length 150
      
  else 
    
    BamCheck "${Basename}"_dedup.bam
    if [[ ${SePe2} == "se" ]]; then BAMF="BAM"; isPE=""; else BAMF="BAMPE"; isPE="-p"; fi
    
    macs2 callpeak \
     --tempdir ./ \
     -t "${Basename}"_dedup.bam \
     -n "${Basename}" -g ${gflag} \
     --keep-dup=all -f ${BAMF}
    
  fi
  
  #/ Case no peaks or some kind of error
  if [[ ! -e ${Basename}_peaks.narrowPeak ]]; then 
    echo 'Error during peak calling -- no narrowPeak file found'
    paste <(echo $Basename) <(echo 'ERROR') >> FRiPs.txt
    exit 1
  else
    if (( $(wc -l "${Basename}"_peaks.narrowPeak | cut -d " " -f1) == 0 )); then
      echo 'Error during peak calling -- no peaks found in narrowPeak file'
      paste <(echo $Basename) <(echo 'ERROR') >> FRiPs.txt
      exit 1
    fi
  fi  
  
  #/ If ok make countmatrix:
  awk 'OFS="\t" {print $1"_"$2+1"_"$3, $1, $2+1, $3, "+"}' "${Basename}"_peaks.narrowPeak > "${Basename}"_peaks.saf
  
  featureCounts ${isPE} -a "${Basename}"_peaks.saf -T ${threadsSort} -F SAF \
    -o "${Basename}"_peaks.counts "${Basename}"_dedup.bam
  
  ASSIGNED=$(grep -w 'Assigned' "${Basename}"_peaks.counts.summary | cut -f2)
  UNASSIGNED=$(grep -w 'Unassigned_NoFeatures' "${Basename}"_peaks.counts.summary | cut -f2)
    
  paste <(echo "${Basename}") <(bc <<< "scale=6;${ASSIGNED}/(${ASSIGNED}+${UNASSIGNED})") >> FRiPs.txt
  
  (>&2 paste -d " " <(echo '[INFO] Peaks/FRiPs for' $Basename 'ended on') <(date))
  (>&2 echo '')
  
}; export -f Peaks

#------------------------------------------------------------------------------------------------------------------------

####//// Start alignment / filtering:  
if [[ ${noalignment} == "FALSE" ]]; then

  if [[ "${format}" == "fq_se" ]]; then 
    FileDefine=$(ls *.fastq.gz 2> /dev/null | awk -F ".fastq.gz" '{print $1 | "sort -u"}' | tr " " "\n")
  fi

  if [[ "${format}" == "fq_pe" ]]; then 
    FileDefine=$(ls *_1.fastq.gz 2> /dev/null | awk -F "_1.fastq.gz" '{print $1 | "sort -u"}' | tr " " "\n")
  fi
    
  if [[ "${format}" == "bam_se" ]] || [[ "${format}" == "bam_pe" ]] ; then 
    FileDefine=$(ls *_ubam.bam 2> /dev/null | awk -F "_ubam.bam" '{print $1 | "sort -u"}' | tr " " "\n")
  fi  
    
  if [[ $(echo $FileDefine | awk NF | wc -l) == 0 ]]; then 
    echo "[Error]: No input files found!"
    exit 1
  fi
  
  #/ Run pipeline:
  echo "${FileDefine}" \
  | parallel -j "${jobnumber}" "Fq2Bam {} 2> {}.log"
  
fi  

####//// Peaks/QC/FRiP:
if [[ ${nofrips} == "FALSE" ]]; then

  if [[ "${atacseq}" == "TRUE" ]]; then Suffix="_cutsites.bed.gz"; fi
  if [[ "${atacseq}" == "FALSE" ]]; then Suffix="_dedup.bam"; fi
  
  ls *"${Suffix}" \
  | awk -F "${Suffix}" '{print $1 | "sort -u"}' \
  | parallel -j "${fripqcJobs}" "Peaks {} 2>> {}.log"
  
fi
