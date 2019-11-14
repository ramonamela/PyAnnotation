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
echo ${SCRIPT_DIR}

#
# HELPER METHODS
#

install_suse_dependencies() {
  sudo zypper update
  sudo zypper -y --non-interactive --no-recommends install libmysqlclient-dev libgd-dev bioperl expect
}

install_ubuntu_dependencies() {
  sudo apt-get update
  sudo apt-get -y --no-install-recommends --no-install-recommends  install libmysqlclient-dev libgd-dev bioperl expect
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
  #spawn perl INSTALL.pl
  #expect "Do you wish to exit so you can get updates (y) or continue (n):"
  #send -- "n"
  #expect "Do you want to install any cache files (y/n)?"
  #send -- "n"
  #expect "Do you want to install any FASTA files (y/n)?"
  #send -- "n"
  #expect "Do you want to install any plugins (y/n)?"
  #send -- "n"  
  popd
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
  install_vep ${SCRIPT_DIR}/ensembl-vep
}

#
# ENTRY POINT
#

main "$@"
