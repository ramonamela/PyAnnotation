#!/bin/bash

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
BASE_DIR="${SCRIPT_DIR}/../"

docker_image_name=${1:-"pyannotationimage"}
echo "${BASE_DIR}/Dockerfile"

sudo docker build -f "${BASE_DIR}/Dockerfile" -t "${docker_image_name}" .

