#!/bin/bash
set -e

tmpdir=$(mktemp -d)
pushd "$tmpdir"
wget https://sourceforge.net/projects/snpeff/files/snpEff_v4_3t_core.zip
unzip snpEff_v4_3t_core.zip
mkdir -p "$HOME/.local/bin"
if [[ :$PATH: != *:"$HOME/.local/bin":* ]] ; then
  echo 'export $PATH:"$HOME/.local/bin"' >> "$HOME/.profile"
fi
mkdir -p "$HOME/.local/etc/"
scp -r snpEff "$HOME/.local/etc/"
chmod +x "$HOME/.local/etc/snpEff/snpEff.jar" "$HOME/.local/etc/snpEff/SnpSift.jar"
ln -s "$HOME/.local/etc/snpEff/snpEff.jar" "$HOME/.local/bin/SnpEff.jar"
ln -s "$HOME/.local/etc/snpEff/SnpSift.jar" "$HOME/.local/bin/SnpSift.jar"
popd
rm -r "$tmpdir"
