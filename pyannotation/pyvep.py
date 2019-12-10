#!/usr/bin/env python3
import sys
import os
from os.path import dirname
from os.path import abspath
from os.path import devnull
from subprocess import Popen, PIPE
from pyannotation.utils import bgzip_tabix_file
from pyannotation.utils import open_potential_gz
import pandas as pd


def annotate_vep(in_file, out_file, clinvar_fields=[], dbNSFP_fields=[], dbSNP_fields=[]):
    compress = False
    if out_file[-3:] == ".gz":
        out_file = out_file[:-3]
        compress = True
    module_path = dirname(dirname(abspath(__file__)))
    command = [module_path + "/dependencies/ensembl-vep/vep", "-i", in_file, "-o", out_file]
    # Reference genome
    command.extend("--species homo_sapiens --assembly GRCh38".split(" "))
    # Clinvar
    if len(clinvar_fields) > 0:
        command.extend(
            ["--custom", module_path + "/data/cache/clinvar/clinvar.vcf.gz,ClinVar_ID,vcf,exact,0," + ",".join(clinvar_fields)])

    # dbSNP
    if len(dbSNP_fields) > 0:
        command.extend(["--custom", module_path + "/data/cache/dbSNP/All_20180418.vcf.gz,dbSNP_ID,vcf,exact,0," + ",".join(["TOPMED"])])

    # dbNSFP
    if len(dbNSFP_fields) > 0:
        command.extend(
            ["--plugin", "dbNSFP," + module_path + "/data/cache/dbNSFP/dbNSFP4.0a.txt.gz," + ",".join(dbNSFP_fields)])
    # General options
    command.extend((
            "--cache --offline --force_overwrite --vcf_info_field ANN --vcf --hgvsg --hgvs --fasta "
            "" + "~/.vep/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz" + " --symbol").split(
        " "))
    null_out = open(devnull, 'w')
    p = Popen(" ".join(command), stdout=null_out, shell=True)
    p.wait()
    if compress:
        bgzip_tabix_file(out_file, out_file=out_file + ".gz")


if __name__ == "__main__":
    in_file = sys.argv[1]
    out_file = sys.argv[2]
    annotate_vep(in_file, out_file)
    #, clinvar_fields=["CLNDISDB", "CLNDISDBINCL", "CLNDN", "CLNSIG", "CLNREVSTAT"])
    """,dbNSFP_fields=["MutationAssessor_score", "MutationAssessor_pred", "SIFT_score", "SIFT_pred",
                                "FATHMM_score", "FATHMM_pred", "Polyphen2_HDIV_score",
                                "Polyphen2_HDIV_pred", "Polyphen2_HVAR_score",
                                "Polyphen2_HVAR_pred"])"""
