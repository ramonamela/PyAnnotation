#!/usr/bin/env python3
import sys
from pyannotation.utils import open_potential_gz
from os.path import dirname
from os.path import abspath
from subprocess import Popen, PIPE
from pyannotation.utils import bgzip_tabix_file
from pyannotation.utils import is_gz_file

def __annotate_stream_snpeff(input_stream, output_stream, clinvar_fields=[], dbNSFP_fields=[]):

    module_path = dirname(dirname(abspath(__file__)))

    # Reference genome
    command_string = "java -Xmx64g -jar " + module_path + "/dependencies/snpEff/snpEff.jar GRCh38.86"
    if len(dbNSFP_fields) > 0 or len(clinvar_fields) > 0:
        p = Popen(["bash", "-c", command_string], stdin=input_stream, stdout=PIPE)
        p_out = p.stdout
    else:
        p = Popen(["bash", "-c", command_string], stdin=input_stream, stdout=output_stream)

    # dbNSFP
    if len(dbNSFP_fields) > 0:
        command_string = "java -Xmx64g -jar " + module_path + "/dependencies/snpEff/SnpSift.jar dbnsfp -f '" + ",".join(dbNSFP_fields) + "' -db " + module_path + "/data/cache/dbNSFP/dbNSFP4.0a.txt.gz" + " -m -a /dev/stdin"
        if len(clinvar_fields) > 0:
            p = Popen(["bash", "-c", command_string], stdin=p_out, stdout=PIPE)
            p_out = p.stdout
        else:
            p = Popen(["bash", "-c", command_string], stdin=p_out, stdout=output_stream)

    # Clinvar
    if len(clinvar_fields) > 0:
        command_string = "java -Xmx64g -jar " + module_path + "/dependencies/snpEff/SnpSift.jar annotate -info '" + ",".join(clinvar_fields) + "' " + module_path + "/data/cache/clinvar/clinvar.vcf.gz /dev/stdin -a"
        p = Popen(["bash", "-c", command_string], stdin=p_out, stdout=output_stream)
    p.wait()

def annotate_snpeff(input_file, out_file, clinvar_fields=[], dbNSFP_fields=[]):
    compress = False
    if out_file[-3:] == ".gz":
        out_file = out_file[:-3]
        compress = True
    #with open_potential_gz(input_file) as input_stream, open(out_file, "w") as output_stream:
    with open(out_file, "w") as output_stream:
        if is_gz_file(input_file):
            p = Popen(["bash", "zcat", input_file], stdout=PIPE)
            input_stream = p.stdout
        else:
            p = Popen(["bash", "cat", input_file], stdout=PIPE)
            input_stream = p.stdout
        __annotate_stream_snpeff(input_stream, output_stream, clinvar_fields=clinvar_fields, dbNSFP_fields=dbNSFP_fields)
    if compress:
        bgzip_tabix_file(out_file, out_file=out_file + ".gz")

if __name__ == "__main__":
    in_file = sys.argv[1]
    out_file = sys.argv[2]
    annotate_snpeff(in_file, out_file, clinvar_fields=["CLNDISDB", "CLNDISDBINCL", "CLNDN", "CLNSIG", "CLNREVSTAT"],
                 dbNSFP_fields=["MutationAssessor_score", "MutationAssessor_pred", "SIFT_score", "SIFT_pred",
                                "FATHMM_score", "FATHMM_pred", "Polyphen2_HDIV_score", "Polyphen2_HDIV_rankscore",
                                "Polyphen2_HDIV_pred", "Polyphen2_HVAR_score", "Polyphen2_HVAR_rankscore",
                                "Polyphen2_HVAR_pred"])
