#!/bin/bash
# bash strict mode
set -euo pipefail
IFS=$'\n\t'

#-------------------------------------------------------------------------------------------------

REF=$1
REF_DIR="$(dirname ${REF})"

# bail gracefully if de-novo failed
if [ ! -s ${REF} ]; then
    echo "De Novo assembly failed" > ${scratch}/guided.err
    exit 0
fi

READ_1=$2
WELL="$(basename -s ".ecc.fq.gz" ${READ_1})"
OUT_DIR="$(dirname ${READ_1})"

THREADS="threads=1"
MEM="-Xmx8g"

# send all temp files to the scratch directory
# trap ensures that scratch will be deleted no matter what
scratch="$(mktemp -d -t tmp.XXXXXXXXXX)"
function finish {
    mv -f ${scratch}/guided.err ${OUT_DIR}/${WELL}.guided.err
    rm -rf ${scratch}
}
trap finish EXIT

#-------------------------------------------------------------------------------------------------

# set maxindel to desired bp length of guaranteed indel call
bbmap.sh \
    ${MEM} \
    ${THREADS} \
    -eoom \
    in=${READ_1} \
    ref=${REF} \
    outm=stdout.sam \
    nodisk \
    overwrite \
    maxindel=16000 \
    2>> ${scratch}/guided.err \
    | samtools sort > ${scratch}/map.bam

samtools index ${scratch}/map.bam

echo ">>>map Complete" >> ${scratch}/guided.err

#-------------------------------------------------------------------------------------------------

# freebayes caller
freebayes \
    -f ${REF} \
    --pooled-continuous \
    --ploidy 1 \
    --haplotype-length 1 \
    --min-base-quality 20 \
    --min-alternate-fraction 0.2 \
    --min-alternate-count 1 \
    --min-mapping-quality 10 \
    ${scratch}/map.bam \
    | bcftools convert -Oz \
    > ${scratch}/freebayes.vcf.gz

bcftools index ${scratch}/freebayes.vcf.gz

echo ">>>freebayes Complete" >> ${scratch}/guided.err

#-------------------------------------------------------------------------------------------------

# generate consensus sequence
bcftools consensus -f ${REF} ${scratch}/freebayes.vcf.gz > ${scratch}/consensus.fasta 2> /dev/null

echo ">>>consensus Complete" >> ${scratch}/guided.err

#-------------------------------------------------------------------------------------------------

# move relevant files
if [ -e ${scratch}/map.bam ]; then
   mv -f ${scratch}/map.bam ${OUT_DIR}/${WELL}.map.bam
fi
if [ -e ${scratch}/map.bam.bai ]; then
   mv -f ${scratch}/map.bam.bai ${OUT_DIR}/${WELL}.map.bam.bai
fi
if [ -e ${scratch}/freebayes.vcf.gz ]; then
   mv -f ${scratch}/freebayes.vcf.gz ${OUT_DIR}/${WELL}.freebayes.vcf.gz
fi
if [ -e ${scratch}/consensus.fasta ]; then
   mv -f ${scratch}/consensus.fasta ${OUT_DIR}/${WELL}.consensus.fasta
fi

echo ">>>mv Complete" >> ${scratch}/guided.err

