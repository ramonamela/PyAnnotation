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
  sudo zypper -y --non-interactive --no-recommends install libmysqlclient-dev libgd-dev bioperl
}

install_ubuntu_dependencies() {
  sudo apt-get update
  sudo apt-get -y --no-install-recommends --no-install-recommends  install libmysqlclient-dev libgd-dev bioperl
}

install_general_dependencies() {
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
  dist=$(cat /etc/os-release | grep ID_LIKE | tr "=" "\t" | awk '{ print $2 }' | tr -d "\"")
  case "${dist}" in
    ubuntu)
      install_ubuntu_dependencies
      ;;
    suse)
      install_suse_dependencies
      ;;
    esac
}

install_vep() {
  install_dir=$1
  mkdir -p ${install_dir}
  pushd ${install_dir}
  cpan App::cpanminus
  cpanm Exception
  cpanm Test::Pod
  cpanm Test::Pod::Coverage
  cpanm GD
  cpanm MIROD/XML-DOM-XPath-0.13.tar.gz
  cpanm Bio::DB::BioFetch
  cpanm Bio::Root::Root
  cpanm Bio::DB::BioFetch
  cpanm Bio::DB::WebDBSeqI
  cpanm Bio::SeqFeature::Lite
  cpan Bio::DB:HTS
  cpanm DBI
  cpanm DBD::mysql
  cpanm Archive::Zip 
  echo "export PATH=\${PATH}:${install_dir}" >> ~/.bashrc
  echo "export PERL5LIB=\${PERL5LIB}:${install_dir}/dependencies/ensembl-vep/modules/:${install_dir}/dependencies/ensembl-vep/" >> ~/.bashrc
  source ~/.bashrc
  echo "n\nn\nn\nn" | perl INSTALL.pl
  popd
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
  install_snpEff ${BASE_DIR}/dependencies/snpEff
}

#
# ENTRY POINT
#

main "$@"
