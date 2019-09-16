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
scp snpEff/snpEff.jar "$HOME/.local/bin/SnpEff.jar"
scp snpEff/SnpSift.jar "$HOME/.local/bin/SnpSift.jar"
chmod +x "$HOME/.local/bin/SnpEff.jar" "$HOME/.local/bin/SnpSift.jar"
popd
rm -r "$tmpdir"
