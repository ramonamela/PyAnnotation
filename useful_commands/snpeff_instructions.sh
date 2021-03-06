#!/bin/bash

#
# SCRIPT GLOBAL VARIABLES
#

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
BASE_DIR="${SCRIPT_DIR}/../"

#
# FULL LOCAL RUN
#

#java -Xmx64g -jar /gpfs/projects/bsc05/apps/MN4/SNPEFF/SRC/snpEff/SnpSift.jar dbnsfp -f '1000Gp3_AC,1000Gp3_AF,1000Gp3_AFR_AC,1000Gp3_AFR_AF,1000Gp3_EUR_AC,1000Gp3_EUR_AF,1000Gp3_AMR_AC,1000Gp3_AMR_AF,1000Gp3_EAS_AC,1000Gp3_EAS_AF,1000Gp3_SAS_AC,1000Gp3_SAS_AF,ESP6500_AA_AC,ESP6500_AA_AF,ESP6500_EA_AC,ESP6500_EA_AF,ExAC_AC,ExAC_AF,ExAC_Adj_AC,ExAC_Adj_AF,ExAC_AFR_AC,ExAC_AFR_AF,ExAC_AMR_AC,ExAC_AMR_AF,ExAC_EAS_AC,ExAC_EAS_AF,ExAC_FIN_AC,ExAC_FIN_AF,ExAC_NFE_AC,ExAC_NFE_AF,ExAC_SAS_AC,ExAC_SAS_AF,ExAC_nonTCGA_AC,ExAC_nonTCGA_AF,ExAC_nonTCGA_Adj_AC,ExAC_nonTCGA_Adj_AF,ExAC_nonTCGA_AFR_AC,ExAC_nonTCGA_AFR_AF,ExAC_nonTCGA_AMR_AC,ExAC_nonTCGA_AMR_AF,ExAC_nonTCGA_EAS_AC,ExAC_nonTCGA_EAS_AF,ExAC_nonTCGA_FIN_AC,ExAC_nonTCGA_FIN_AF,ExAC_nonTCGA_NFE_AC,ExAC_nonTCGA_NFE_AF,ExAC_nonTCGA_SAS_AC,ExAC_nonTCGA_SAS_AF,ExAC_nonpsych_AC,ExAC_nonpsych_AF,ExAC_nonpsych_Adj_AC,ExAC_nonpsych_Adj_AF,ExAC_nonpsych_AFR_AC,ExAC_nonpsych_AFR_AF,ExAC_nonpsych_AMR_AC,ExAC_nonpsych_AMR_AF,ExAC_nonpsych_EAS_AC,ExAC_nonpsych_EAS_AF,ExAC_nonpsych_FIN_AC,ExAC_nonpsych_FIN_AF,ExAC_nonpsych_NFE_AC,ExAC_nonpsych_NFE_AF,ExAC_nonpsych_SAS_AC,ExAC_nonpsych_SAS_AF,gnomAD_exomes_AC,gnomAD_exomes_AN,gnomAD_exomes_AF,gnomAD_exomes_AFR_AC,gnomAD_exomes_AFR_AN,gnomAD_exomes_AFR_AF,gnomAD_exomes_AMR_AC,gnomAD_exomes_AMR_AN,gnomAD_exomes_AMR_AF,gnomAD_exomes_ASJ_AC,gnomAD_exomes_ASJ_AN,gnomAD_exomes_ASJ_AF,gnomAD_exomes_EAS_AC,gnomAD_exomes_EAS_AN,gnomAD_exomes_EAS_AF,gnomAD_exomes_FIN_AC,gnomAD_exomes_FIN_AN,gnomAD_exomes_FIN_AF,gnomAD_exomes_NFE_AC,gnomAD_exomes_NFE_AN,gnomAD_exomes_NFE_AF,gnomAD_exomes_SAS_AC,gnomAD_exomes_SAS_AN,gnomAD_exomes_SAS_AF,gnomAD_exomes_OTH_AC,gnomAD_exomes_OTH_AN,gnomAD_exomes_OTH_AF,gnomAD_genomes_AC,gnomAD_genomes_AN,gnomAD_genomes_AF,gnomAD_genomes_AFR_AC,gnomAD_genomes_AFR_AN,gnomAD_genomes_AFR_AF,gnomAD_genomes_AMR_AC,gnomAD_genomes_AMR_AN,gnomAD_genomes_AMR_AF,gnomAD_genomes_ASJ_AC,gnomAD_genomes_ASJ_AN,gnomAD_genomes_ASJ_AF,gnomAD_genomes_EAS_AC,gnomAD_genomes_EAS_AN,gnomAD_genomes_EAS_AF,gnomAD_genomes_FIN_AC,gnomAD_genomes_FIN_AN,gnomAD_genomes_FIN_AF,gnomAD_genomes_NFE_AC,gnomAD_genomes_NFE_AN,gnomAD_genomes_NFE_AF,gnomAD_genomes_OTH_AC,gnomAD_genomes_OTH_AN,gnomAD_genomes_OTH_AF,CADD_phred,CADD_raw,Polyphen2_HDIV_score,Polyphen2_HDIV_pred,Polyphen2_HVAR_score,Polyphen2_HVAR_pred,SIFT_score,SIFT_pred,MutationAssessor_score,MutationAssessor_pred,MutationTaster_score,MutationTaster_pred,FATHMM_score,FATHMM_pred,GERP++_NR,GERP++_RS,LRT_Omega,LRT_score,LRT_pred,MetaLR_score,MetaLR_pred,MetaSVM_score,MetaSVM_pred,Reliability_index,PROVEAN_score,PROVEAN_pred,VEST3_score,M-CAP_score,M-CAP_pred,fathmm-MKL_coding_score,fathmm-MKL_coding_pred,Eigen_coding_or_noncoding,Eigen-raw,Eigen-phred,Eigen-PC-raw,Eigen-PC-phred,integrated_fitCons_score,integrated_confidence_value,GM12878_fitCons_score,GM12878_confidence_value,H1-hESC_fitCons_score,H1-hESC_confidence_value,HUVEC_fitCons_score,HUVEC_confidence_value,GenoCanyon_score,DANN_score,MutPred_score,MutPred_Top5features,REVEL_score,phastCons100way_vertebrate,phastCons20way_mammalian,phyloP100way_vertebrate,phyloP20way_mammalian,SiPhy_29way_pi,SiPhy_29way_logOdds,cds_strand,codonpos,refcodon,codon_degeneracy,Ancestral_allele,Interpro_domain,aapos,Uniprot_aapos_Polyphen2,Uniprot_acc_Polyphen2,Uniprot_id_Polyphen2' -db /gpfs/projects/bsc05/apps/MN4/SNPEFF/SRC/snpEff/db/GRCh37/dbNSFP3.5a/dbNSFP3.5a_hg19.txt.gz -m -a $results_combine/snvs$trim/snvs.PASS.vcf | \
#java -Xmx64g -jar /gpfs/projects/bsc05/apps/MN4/SNPEFF/SRC/snpEff/SnpSift.jar annotate -info 'CNT,SOMATIC_STATUS,PRIMARY_SITE,PRIMARY_HISTOLOGY' /gpfs/projects/bsc05/apps/reference_data/COSMIC/Cosmic_v84.vcf /dev/stdin -a | \
#java -Xmx64g -jar /gpfs/projects/bsc05/apps/MN4/SNPEFF/SRC/snpEff/SnpSift.jar annotate -id -noInfo /gpfs/projects/bsc05/apps/MN4/SNPEFF/SRC/snpEff/db/GRCh37/00-All.vcf.gz /dev/stdin -a | \
#java -Xmx64g -jar /gpfs/projects/bsc05/apps/MN4/SNPEFF/SRC/snpEff/SnpSift.jar annotate -id -info 'CLNDISDB,CLNDISDBINCL,CLNDN,CLNSIG,CLNREVSTAT' /gpfs/projects/bsc05/apps/MN4/SNPEFF/SRC/snpEff/db/GRCh37/clinvar.vcf.gz /dev/stdin -a | \
#java -Xmx64g -jar /gpfs/projects/bsc05/apps/MN4/SNPEFF/SRC/snpEff/snpEff.jar GRCh37.p13.RefSeq > $results_combine/snvs$trim/snvs.PASS.ANNOT.GRCh37.p13.RefSeq.vcf 

#
# SUBSET OF COMMANDS THAT WE ARE INTERESTED IN
#

## Scores taken into account
#MutationAssessor_score
#MutationAssessor_pred
#SIFT_score
#SIFT_pred
#FATHMM_score
#FATHMM_pred

#java -Xmx64g -jar /gpfs/projects/bsc05/apps/MN4/SNPEFF/SRC/snpEff/SnpSift.jar dbnsfp -f 'MutationAssessor_score,MutationAssessor_pred,SIFT_score,SIFT_pred,FATHMM_score,FATHMM_pred' -db /gpfs/projects/bsc05/apps/MN4/SNPEFF/SRC/snpEff/db/GRCh37/dbNSFP3.5a/dbNSFP3.5a_hg19.txt.gz -m -a input.vcf | \
#java -Xmx64g -jar /gpfs/projects/bsc05/apps/MN4/SNPEFF/SRC/snpEff/SnpSift.jar annotate -id -info 'CLNDISDB,CLNDISDBINCL,CLNDN,CLNSIG,CLNREVSTAT' /gpfs/projects/bsc05/apps/MN4/SNPEFF/SRC/snpEff/db/GRCh37/clinvar.vcf.gz /dev/stdin -a | \
#java -Xmx64g -jar /gpfs/projects/bsc05/apps/MN4/SNPEFF/SRC/snpEff/snpEff.jar GRCh37.p13.RefSeq > $results_combine/snvs$trim/snvs.PASS.ANNOT.GRCh37.p13.RefSeq.vcf

## By now, we don't support this one since it's not annotated in the ICGC portal ##
##??##java -Xmx64g -jar /gpfs/projects/bsc05/apps/MN4/SNPEFF/SRC/snpEff/SnpSift.jar annotate -id -noInfo /gpfs/projects/bsc05/apps/MN4/SNPEFF/SRC/snpEff/db/GRCh37/00-All.vcf.gz /dev/stdin -a

#
# HELPER METHODS
#

available_punctuations() {
  cat ${BASE_DIR}/useful_commands/snpeff_instructions.sh | head -n3 | tail -n1 | awk '{ print $7 }' | tr "," "\n" | grep score
}

get_all_options() {
  option=${1}
  cat ${BASE_DIR}/useful_commands/snpeff_instructions.sh | head -n3 | tail -n1 | awk '{ print $7 }' | tr "," "\n" | grep ${option}
}

annotate_file() {
  input_file=${1}
  output_file=${2}
  compress=false
  if [[ ${output_file} == *vcf.gz ]]; then
    compress=true
    output_file=${output_file::-3}
  elif [[ ! ${output_file} == *vcf ]]; then
    false
  fi
  if [[ ${input_file} == *vcf.gz ]]; then
    zcat ${input_file} | \
    java -Xmx64g -jar ${BASE_DIR}/dependencies/snpEff/SnpSift.jar dbnsfp -f 'MutationAssessor_score,MutationAssessor_pred,SIFT_score,SIFT_pred,FATHMM_score,FATHMM_pred' -db ${BASE_DIR}/data/cache/dbNSFP4.0a.txt.gz -m /dev/stdin -a | \
    java -Xmx64g -jar ${BASE_DIR}/dependencies/snpEff/SnpSift.jar annotate -info 'CLNDISDB,CLNDISDBINCL,CLNDN,CLNSIG,CLNREVSTAT' ${BASE_DIR}/data/cache/clinvar.vcf.gz /dev/stdin -a | \
    java -Xmx64g -jar ${BASE_DIR}/dependencies/snpEff/snpEff.jar GRCh38.86 > ${output_file}
  elif [[ ${input_file} == *vcf ]]; then
    java -Xmx64g -jar ${BASE_DIR}/dependencies/snpEff/SnpSift.jar dbnsfp -f 'MutationAssessor_score,MutationAssessor_pred,SIFT_score,SIFT_pred,FATHMM_score,FATHMM_pred' -db ${BASE_DIR}/data/cache/dbNSFP4.0a.txt.gz -m -a ${input_file} | \
    java -Xmx64g -jar ${BASE_DIR}/dependencies/snpEff/SnpSift.jar annotate -id -info 'CLNDISDB,CLNDISDBINCL,CLNDN,CLNSIG,CLNREVSTAT' ${BASE_DIR}/data/cache/clinvar.vcf.gz /dev/stdin -a | \
    java -Xmx64g -jar ${BASE_DIR}/dependencies/snpEff/snpEff.jar GRCh38.86 > ${output_file}
  else
    false
  fi
  if [[ ${compress} = true ]]; then
    bgzip -f "${output_file}"
    tabix -f "${output_file}.gz"
  fi
}

#
# MAIN METHOD
#

main() {
  annotate_file ${BASE_DIR}/data/examples/example.vcf.gz ${BASE_DIR}/data/examples/example_snpeff.vcf.gz
}

#
# ENTRY POINT
#

main "$@"
