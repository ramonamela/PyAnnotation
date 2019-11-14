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
  if [ ! -f "${HOME}/.vep/homo_sapiens_vep_94_GRCh38.tar.gz" ]; then
    mkdir -p "${HOME}/.vep"
    pushd "${HOME}/.vep"
      wget ftp://ftp.ensembl.org/pub/release-94/variation/VEP/homo_sapiens_vep_94_GRCh38.tar.gz
      pushd ${BASE_DIR}/dependencies/ensembl-vep
      perl convert_cache.pl -species "Homo_sapiens" -version 94
      popd
    popd
  fi
  #${BASE_DIR}/dependencies/ensembl-vep/vep 
}

#
# ENTRY POINT
#

main "$@"
