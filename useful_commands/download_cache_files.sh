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

## In order to have everything together, it is recommended to have a symlink of ~/.vep to ${BASE_DIR}/data/cache/vep_data
## and snpEff/data to ${BASE_DIR}/data/cache/snpEff_data
## This is done by default when installing vep and snpEff with this scripts
## More precisely, there are the commands used to do so:
## ln -s "${BASE_DIR}/data/cache/vep_data" "${HOME}/.vep"
## ln -s "${BASE_DIR}/data/cache/snpEff_data" "${install_dir}/data"

download_vep_files() {
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
  fi
	if [ ! -d "${HOME}/.vep/homo_sapiens/dna/" ]; then
		mkdir -p "${HOME}/.vep/homo_sapiens/dna/"
	fi
	if [ ! -f "${HOME}/.vep/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz" ]; then
		pushd "${HOME}/.vep/homo_sapiens/dna/"
		curl -O ftp://ftp.ensembl.org/pub/release-94/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
		gzip -d Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
		bgzip Homo_sapiens.GRCh38.dna.primary_assembly.fa
		rm -f Homo_sapiens.GRCh38.dna.primary_assembly.fa
		popd
	fi
}

download_snpeff_files() {
  ## Cache files for GRCh38 - It seems it's the version 86. Check for updates
  if [ ! -d "${BASE_DIR}/data/cache/snpEff_data" ]; then
    mkdir -p "${BASE_DIR}/data/cache/snpEff_data"
  fi
  ## Cache files for GRCh38
  if [ ! -f "${BASE_DIR}/data/cache/snpEff_data/snpEff_v4_3_GRCh38.86.zip" ]; then
    pushd "${BASE_DIR}/data/cache/snpEff_data/"
    wget https://sourceforge.net/projects/snpeff/files/databases/v4_3/snpEff_v4_3_GRCh38.86.zip
    popd
  fi
  if [ ! -d "${BASE_DIR}/data/cache/snpEff_data/GRCh38.86" ]; then
    pushd "${BASE_DIR}/data/cache/snpEff_data/"
    unzip -o snpEff_v4_3_GRCh38.86.zip
    rm snpEff_v4_3_GRCh38.86.zip
    mv ./data/GRCh38.86 .
    rm -r data
    popd
  fi
}

download_common_files() {
  ## CLINVAR ##
	if [ ! -d "${BASE_DIR}/data/cache/clinvar/clinvar.vcf.gz" ]; then
		mkdir -p "${BASE_DIR}/data/cache/clinvar"
	fi
  if [ ! -f "${BASE_DIR}/data/cache/clinvar/clinvar.vcf.gz" ]; then
    pushd "${BASE_DIR}/data/cache/clinvar"
    wget ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
    wget ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi
    wget ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.md5
    popd
  fi

  ## dbNSFP ##
	if [ ! -d "${BASE_DIR}/data/cache/dbNSFP" ]; then
		mkdir -p "${BASE_DIR}/data/cache/dbNSFP"
	fi
  if [ ! -f "${BASE_DIR}/data/cache/dbNSFP/dbNSFP4.0a.txt.gz" ]; then
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
  	if [ ! -f "${BASE_DIR}/data/cache/dbNSFP/dbNSFP4.0a.zip" ]; then
    	pushd "${BASE_DIR}/data/cache/dbNSFP"
    	pip3 install gdown
    	gdown https://drive.google.com/uc?id=1BNLEdIc4CjCeOa7V7Z8n8P8RHqUaF5GZ
	popd
  fi
    pushd "${BASE_DIR}/data/cache/dbNSFP"
    unzip -u dbNSFP4.0a.zip
    zcat dbNSFP4.0a_variant.chr1 | head -n1 > dbNSFP4.0a.txt
    zcat dbNSFP4.0a_variant.chr* | grep -v "#" >> dbNSFP4.0a.txt
    rm dbNSFP4.0a_variant.chr*
    bgzip dbNSFP4.0a.txt
    tabix -s 1 -b 2 -e 2 dbNSFP4.0a.txt.gz
    rm dbNSFP4.0a.txt
    popd
  fi
	## CiViC ##
	if [ ! -d "${BASE_DIR}/data/cache/civic" ]; then
		mkdir "${BASE_DIR}/data/cache/civic"
	fi
	if [ ! -f "${BASE_DIR}/data/cache/civic/nightly-GeneSummaries.tsv" ]; then
		pushd "${BASE_DIR}/data/cache/civic"
		wget https://civicdb.org/downloads/nightly/nightly-GeneSummaries.tsv
		popd
	fi
	if [ ! -f "${BASE_DIR}/data/cache/civic/nightly-VariantSummaries.tsv" ]; then
		pushd "${BASE_DIR}/data/cache/civic"
    wget https://civicdb.org/downloads/nightly/nightly-VariantSummaries.tsv
    popd
	fi
  if [ ! -f "${BASE_DIR}/data/cache/civic/nightly-ClinicalEvidenceSummaries.tsv" ]; then
    pushd "${BASE_DIR}/data/cache/civic"
    wget https://civicdb.org/downloads/nightly/nightly-ClinicalEvidenceSummaries.tsv
    popd
  fi
  if [ ! -f "${BASE_DIR}/data/cache/civic/nightly-AssertionSummaries.tsv" ]; then
    pushd "${BASE_DIR}/data/cache/civic"      
    wget https://civicdb.org/downloads/nightly/nightly-AssertionSummaries.tsv
    popd
  fi
	## dbSNP ##
	base_path_dbSNP="${BASE_DIR}/data/cache/dbSNP/"
	if [ ! -d "${base_path_dbSNP}" ]; then
		mkdir -p "${base_path_dbSNP}"
	fi
	files_to_download_dbsnp=( "All_20180418.vcf.gz" "All_20180418.vcf.gz.md5" "All_20180418.vcf.gz.tbi" )
	base_ftp_path="ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/"	
	for i in "${files_to_download_dbsnp[@]}"; do
		if [ ! -f "${base_path_dbSNP}${i}" ]; then
			pushd ${base_path_dbSNP}
			wget "${base_ftp_path}${i}"
			popd
		fi
	done
	## COSMIC ##
	base_path_cosmic="${BASE_DIR}/data/cache/cosmic/"
  if [ ! -d "${base_path_cosmic}" ]; then
    mkdir -p "${base_path_cosmic}"
  fi
	vcf_files_to_download_cosmic=( "CosmicNonCodingVariants.vcf.gz" "CosmicCodingMuts.vcf.gz" )
	files_to_download_cosmic=( "CosmicSample.tsv.gz" "CosmicCompleteTargetedScreensMutantExport.tsv.gz" "CosmicGenomeScreensMutantExport.tsv.gz" "CosmicMutantExport.tsv.gz" "CosmicMutantExportCensus.tsv.gz" "CosmicNCV.tsv.gz" "CosmicCompleteGeneExpression.tsv.gz" "CosmicResistanceMutations.tsv.gz" )
	files_to_download_cosmic=( "CosmicMutantExportCensus.tsv.gz" )
	cosmic_credentials=$(echo "${COSMIC_USER}:${COSMIC_PASS}" | base64)
	base_url="https://cancer.sanger.ac.uk/cosmic/file_download/GRCh38/cosmic/v90/"
	for i in "${files_to_download_cosmic[@]}"; do
    if [ ! -f "${base_path_cosmic}${i}" ]; then
			pushd ${base_path_cosmic}
			response=$(curl -H "Authorization: Basic ${cosmic_credentials}" ${base_url}${i})
			curl "${response:8:-2}" --output "${i}"
			popd
    fi
		if [ ! -f "${base_path_cosmic}${i}.tbi" ]; then
			pushd ${base_path_cosmic}
			python3 ${BASE_DIR}/useful_commands/transform_cosmic_mutant_export_census.py "${base_path_cosmic}${i}" "${base_path_cosmic}${i::-7}Indexable.tsv"
			bgzip "${base_path_cosmic}${i::-7}Indexable.tsv"
			mv "${base_path_cosmic}${i::-7}Indexable.tsv.gz" "${base_path_cosmic}${i}"
			tabix -s 1 -b 2 -e 2 "${base_path_cosmic}${i}"
			popd
		fi
  done
  for i in "${vcf_files_to_download_cosmic[@]}"; do
    if [ ! -f "${base_path_cosmic}${i}" ]; then
			pushd ${base_path_cosmic}
			response=$(curl -H "Authorization: Basic ${cosmic_credentials}" ${base_url}VCF/${i})
			curl "${response:8:-2}" --output "${i}"
			gunzip ${i}
			bgzip ${i::-3}
			tabix -p vcf ${i}
			popd
    fi
		if [ ! -f "${base_path_cosmic}${i}.tbi" ]; then
			pushd ${base_path_cosmic}
			gunzip ${i}
      bgzip ${i::-3}
      tabix -p vcf ${i}
      popd
		fi
  done
}

#
# MAIN METHOD
#

main() {
  download_common_files
  download_vep_files
  download_snpeff_files
}

#
# ENTRY POINT
#

main "$@"
