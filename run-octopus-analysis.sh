#!/bin/bash

Help()
{
  echo 'Usage: [-h] [-a <absolute path>] [-o <absolute path>]'
  echo '       [-d <docker:image>] -q <absolute path>'
  echo
  echo 'Automatically run through the steps of the OCTOPUS plasmid sequencing pipeline:'
  echo '   1. src/link-data.py'
  echo '   2. src/generate-input.py'
  echo '   3. make'
  echo
  echo 'Paths MUST be absolute paths e.g. /path/to/octopus3, ~/path/to/fastas, $(pwd)'
  echo 'Required argument:'
  echo '   -q   mise[q] run root path: contains SampleSheet.csv and further path to Fastqs'
  echo
  echo 'Optional arguments:'
  echo '   -h   this help text and exit'
  echo '   -a   fast[a] path: contains expected fastas for the run'
  echo '        default is no fastas and performing de-novo alignment only'
  echo '   -o   [o]ctopus3 root path: contains Makefile, pipeline/, src/'
  echo '        default is the absolute path to current directory'
  echo '   -d   [d]ocker image name'
  echo '        default is octant/octopus3:latest'
  echo
}

# process flags
while getopts ha:o:d:q: flag
do
  case "${flag}" in
    h) Help
       exit ;;
    a) FASTAS_DIR="${OPTARG}" ;;
    o) OCTOPUS3_DIR="${OPTARG}" ;;
    d) DOCKER_IMAGE="${OPTARG}" ;;
    q) MISEQ_OUTPUTS="${OPTARG}" ;;
  esac
done

if [ "$#" -lt 2 ]
then
  echo 'Missing arguments! Use help flag -h to check proper usage.' >&2
  exit
fi

if [ -z "$MISEQ_OUTPUTS" ]
then
  echo '-q: mise[q] run root path missing! Use help flag -h to check proper usage.' >&2
  exit
fi

if [ -z "$DOCKER_IMAGE" ]
then
  DOCKER_IMAGE="octant/octopus3"
fi

if [ -z "$OCTOPUS3_DIR" ] || [ "$OCTOPUS3_DIR" == "." ]
then
  OCTOPUS3_DIR="$(pwd)"
fi

# background highlighting colors
RED='\033[41m'
GREEN='\033[42m'
DBLUE='\033[44m'
LBLUE='\033[104m'
NC='\033[0m'

# provide error if pipeline is interrupted after this point
set -e
trap 'echo -e "${RED}Exiting pipeline prematurely.${NC}"' EXIT

# begin pipeline by making working directory the octopus3 root directory
echo -e "${DBLUE}Automated OCTOPUS computational pipeline started!${NC}"
cd "$OCTOPUS3_DIR"

# link sample sheet and fastqs to pipeline
echo -e "${LBLUE}Linking data . . .${NC}"
python3 src/link-data.py "$MISEQ_OUTPUTS" -o pipeline/

# run full analysis if fastas provided and "input.fasta" generated, otherwise perform make denovo
if [ -z "$FASTAS_DIR" ]
then
  echo -e "${RED}No valid input fasta directory provided. Performing de-novo alignment only.${NC}"
  echo -e "${LBLUE}Running de-novo alignment . . .${NC}"
  MAKE_TARGET="spades-contigs.fasta"
else
  echo -e "${LBLUE}Generating "input.fasta" . . .${NC}"
  python3 src/generate-input.py "$FASTAS_DIR" -o "pipeline/$(basename "$MISEQ_OUTPUTS")"
  echo -e "${LBLUE}Running analysis pipeline . . .${NC}"
  MAKE_TARGET="aggregated-stats.tsv"
fi

echo -e "${LBLUE}Running docker image \"$DOCKER_IMAGE\"${NC}"
docker run --rm -it --user $(id -u) -v "$OCTOPUS3_DIR":/data/octopus3 -v "$MISEQ_OUTPUTS":"$MISEQ_OUTPUTS" -w /data/octopus3 $DOCKER_IMAGE bash -c "make pipeline/$(basename "$MISEQ_OUTPUTS")/$MAKE_TARGET"

# end
trap "echo -e \"${GREEN}OCTOPUS analysis pipeline completed successfully!${NC}\"" EXIT
echo -e "${DBLUE}Output \"$MAKE_TARGET\" for $(basename "$MISEQ_OUTPUTS") can be found in $OCTOPUS3_DIR/pipeline/$(basename "$MISEQ_OUTPUTS")${NC}"
