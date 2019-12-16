#!/bin/bash

#
# SCRIPT GLOBAL VARIABLES
#

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
BASE_DIR="${SCRIPT_DIR}/../"

python3 ${BASE_DIR}/pyannotation/pysnpeff.py "${BASE_DIR}/data/examples/small_example.vcf.gz" "${BASE_DIR}/data/examples/small_example_snpeff.vcf.gz" "${BASE_DIR}/config/simple.json"
#python3 ${BASE_DIR}/pyannotation/pysnpeff.py "${BASE_DIR}/data/examples/small_example_norm.vcf.gz" "${BASE_DIR}/data/examples/small_example_norm_snpeff.vcf.gz"
