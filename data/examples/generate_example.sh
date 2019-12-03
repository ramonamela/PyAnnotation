#!/bin/bash

vcf_name=${1:-example.vcf}

cat ${vcf_name} | bgzip -c > "${vcf_name}.gz"
tabix -p vcf "${vcf_name}.gz"
