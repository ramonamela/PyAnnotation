#!/bin/bash

#
# BASH OPTIONS
#

set -e # Exit when command fails
set -u # Exit when undefined variable

#
# SCRIPT GLOBAL VARIABLES
#

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
BASE_DIR="${SCRIPT_DIR}/../"

#
# HELPER METHODS
#


#
# MAIN METHOD
#

main() {
  ## Cache files for GRCh38
  if [ ! -f "${HOME}/.vep/homo_sapiens/94_GRCh38/10/99000001-100000000.gz" ]; then
    mkdir -p "${HOME}/.vep/homo_sapiens/94_GRCh38"
    pushd "${HOME}/.vep"
      wget ftp://ftp.ensembl.org/pub/release-94/variation/VEP/homo_sapiens_vep_94_GRCh38.tar.gz
      tar xzf homo_sapiens_vep_94_GRCh38.tar.gz
      rm -r homo_sapiens_vep_94_GRCh38.tar.gz
      pushd ${BASE_DIR}/dependencies/ensembl-vep
      perl convert_cache.pl -species "Homo_sapiens" -version 94
      popd
    popd
  fi
  if [ ! -f "${HOME}/.vep/homo_sapiens/94_GRCh38/10/all_vars.gz" ]; then
    pushd ${BASE_DIR}/dependencies/ensembl-vep
    perl convert_cache.pl -species "Homo_sapiens" -version 94
    popd
    "${HOME}/.vep/homo_sapiens/94_GRCh38/10/all_vars.gz does not exist"
  fi
  if [ ! -f "${BASE_DIR}/data/cache/clinvar.vcf.gz" ]; then
    pushd "${BASE_DIR}/data/cache/"
    wget ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
    wget ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi
    wget ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.md5
    popd
  fi
  if [ ! -f "${BASE_DIR}/data/cache/dbNSFP4.0a.zip" ]; then
    pushd "${BASE_DIR}/data/cache/"
    pip install gdown
    popd
    gdown https://drive.google.com/uc?id=1BNLEdIc4CjCeOa7V7Z8n8P8RHqUaF5GZ
  fi
  if [ ! -f "${BASE_DIR}/data/cache/dbNSFP4.0a.txt.gz" ]; then
    path_to_executable=$(whereis htsfile | awk '{ print $2 }')
    if [ ! -x "$path_to_executable" ] ; then
      pushd /opt
      sudo wget https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2
      sudo tar -vxjf htslib-1.9.tar.bz2
      sudo rm htslib-1.9.tar.bz2
      pushd htslib-1.9
      sudo make
      echo 'export PATH=$PATH:/opt/htslib-1.9' >> ~/.bashrc
      echo 'export HTSLIB_DIR=/opt/htslib-1.9/' >> ~/.bashrc
      source ~/.bashrc
      popd
    fi
    pushd "${BASE_DIR}/data/cache/"
      unzip -u dbNSFP4.0a.zip
      zcat dbNSFP4.0a_variant.chr1 | head -n1 > dbNSFP4.0a.txt
      zcat dbNSFP4.0a_variant.chr* | grep -v "#" >> dbNSFP4.0a.txt
      rm dbNSFP4.0a_variant.chr*
      bgzip dbNSFP4.0a.txt
      tabix -s 1 -b 2 -e 2 dbNSFP4.0a.txt.gz
      rm dbNSFP4.0a.txt
    popd
  fi
  #${BASE_DIR}/dependencies/ensembl-vep/vep 
}

#
# ENTRY POINT
#

main "$@"
