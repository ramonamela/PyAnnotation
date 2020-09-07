#!/bin/bash

#
# BASH OPTIONS
#

set -e # Exit when command fails
#set -u # Exit when undefined variable

#
# SCRIPT GLOBAL VARIABLES
#

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
BASE_DIR="${SCRIPT_DIR}/../"

#
# HELPER METHODS
#

install_suse_dependencies() {
  sudo zypper update
  sudo zypper -y --non-interactive --no-recommends install libmysqlclient-dev libgd-dev bioperl bcftools python3 python3-pip python3-dev unzip perl perl-base libipc-run-perl libxml-dom-perl libxml-dom-xpath-perl perl-DBD-mysql curl
}

install_ubuntu_dependencies() {
  sudo apt-get update
  sudo apt-get -y --no-install-recommends --no-install-recommends  install libmysqlclient-dev libgd-dev bioperl bcftools python3 python3-pip python3-dev unzip perl perl-base libipc-run-perl libxml-dom-perl libxml-dom-xpath-perl libdbd-mysql-perl curl
}

install_general_dependencies() {
  dist=$(cat /etc/os-release | grep ID_LIKE | tr "=" "\t" | awk '{ print $2 }' | tr -d "\"")
  case "${dist}" in
		debian)
			install_ubuntu_dependencies
			;;
    ubuntu)
      install_ubuntu_dependencies
      ;;
    suse)
      install_suse_dependencies
      ;;
    esac
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
		export PATH=$PATH:/opt/htslib-1.9
		export HTSLIB_DIR=/opt/htslib-1.9/
    source ~/.bashrc
    popd
  fi
	pip3 install setuptools wheel
	pip3 install pytabix nose
	pip3 install jsonschema pandas
}

install_vep() {
  install_dir=$1
  mkdir -p ${install_dir}
  pushd ${install_dir}
  cpan App::cpanminus
  cpanm Module::Build
  cpanm Exception
  cpanm Test::Pod
  cpanm Test::Pod::Coverage
  cpanm GD
  cpanm MIROD/XML-DOM-XPath-0.13.tar.gz
  cpanm Bio::DB::BioFetch
  cpanm Bio::Root::Root
  cpanm Bio::DB::BioFetch
  cpanm Bio::DB::WebDBSeqI
  cpanm Bio::DB::HTS
  cpanm Bio::DB::HTS::VCF
  cpanm Bio::DB::HTS::Tabix
  cpanm Bio::SeqFeature::Lite
  cpanm DBI
  #cpanm DBD::mysql
  cpanm Archive::Zip 
  echo "export PATH=\${PATH}:${install_dir}" >> ~/.bashrc
  echo "export PERL5LIB=${install_dir}/:${install_dir}/modules/:\${PERL5LIB}" >> ~/.bashrc
	export PERL5LIB=${install_dir}/:${install_dir}/modules/:${PERL5LIB}
	export PATH=${PATH}:${install_dir}
  source ~/.bashrc
  #perl INSTALL.pl --NO_HTSLIB --NO_UPDATE --AUTO p
  perl INSTALL.pl --NO_HTSLIB --NO_UPDATE --AUTO p --PLUGINS dbNSFP
  echo "${BASE_DIR}/data/cache/vep_data"
  mkdir -p "${BASE_DIR}/data/cache/vep_data"
  #mkdir -p "${HOME}/.vep/Plugins"
  #pushd "${HOME}/.vep/Plugins"
  #wget https://github.com/Ensembl/VEP_plugins/blob/postreleasefix/94/dbNSFP.pm
  #popd
  popd
  if [[ -L "${HOME}/.vep" && -d "${HOME}/.vep" ]]; then
    rm -r "${HOME}/.vep"
  fi
  if [[ ! -d "${HOME}/.vep" ]]; then
    echo "Symlinks"
    ln -s "${BASE_DIR}/data/cache/vep_data" "${HOME}/.vep"
  fi
  echo "End install VEP"
}

install_snpEff() {
  install_dir=${1}
  if [ ! -d "${install_dir}" ]; then
    pushd /tmp
    wget http://sourceforge.net/projects/snpeff/files/snpEff_v4_3t_core.zip
    unzip snpEff_v4_3t_core.zip
    mv /tmp/snpEff ${install_dir}
    rm snpEff_v4_3t_core.zip
    popd
    mkdir -p "${BASE_DIR}/data/cache/snpEff_data"
    ln -s "${BASE_DIR}/data/cache/snpEff_data" "${install_dir}/data"
  fi
}

#
# MAIN METHOD
#

main() {
  pushd ${SCRIPT_DIR}/../
  git submodule update --init
  popd
  echo "## INSTALL GENERAL DEPENDENCIES ##"
  install_general_dependencies
  echo "## INSTALL ENSEMBL-VEP ##"
  install_vep ${BASE_DIR}/dependencies/ensembl-vep
	echo "## INSTALL SNPEFF ##"
  install_snpEff ${BASE_DIR}/dependencies/snpEff
}

#
# ENTRY POINT
#

main "$@"
