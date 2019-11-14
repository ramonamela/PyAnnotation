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

#
# HELPER METHODS
#

install_vep() {
  install_dir=$1
  pushd ${install_dir}
  rm -rf ensembl-vep
  git clone --branch release/94 https://github.com/Ensembl/ensembl-vep.git
  popd
}

#
# MAIN METHOD
#

main() {
  install_vep ${SCRIPT_DIR}
  #pushd ${SCRIPT_DIR}/../../
  #git submodule update --init
  #popd
}

#
# ENTRY POINT
#

main "$@"
