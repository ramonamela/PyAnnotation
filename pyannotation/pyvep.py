#!/usr/bin/env python3
import sys
import os
from os.path import dirname
from os.path import abspath
from os.path import devnull
from subprocess import Popen, PIPE
from pyannotation.utils import bgzip_tabix_file
from pyannotation.utils import is_gz_file
import io
from os import read, write
import select

def annotate_stream_vep(input_stream, output_stream, clinvar_fields=[], dbNSFP_fields=[], dbSNP_fields=[]):
    module_path = dirname(dirname(abspath(__file__)))
    command = [module_path + "/dependencies/ensembl-vep/vep", "-o", "STDOUT"]
    # Reference genome
    command.extend("--species homo_sapiens --assembly GRCh38".split(" "))
    # Clinvar
    if len(clinvar_fields) > 0:
        command.extend(
            ["--custom",
             module_path + "/data/cache/clinvar/clinvar.vcf.gz,ClinVar_ID,vcf,exact,0," + ",".join(clinvar_fields)])

    # dbSNP
    if len(dbSNP_fields) > 0:
        command.extend(["--custom",
                        module_path + "/data/cache/dbSNP/All_20180418.vcf.gz,dbSNP_ID,vcf,exact,0," + ",".join(
                            ["TOPMED"])])

    # dbNSFP
    if len(dbNSFP_fields) > 0:
        command.extend(
            ["--plugin", "dbNSFP," + module_path + "/data/cache/dbNSFP/dbNSFP4.0a.txt.gz," + ",".join(dbNSFP_fields)])
    # General options
    command.extend((
            "--cache --offline --force_overwrite --format vcf --vcf_info_field ANN --vcf --hgvsg --hgvs --fasta "
            "" + "~/.vep/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz" + " --symbol").split(
        " "))
    print(command)
    p = Popen(" ".join(command), stdin=input_stream, stdout=output_stream, shell=True)
    return p

def annotate_vep(in_file, out_file):
    ## in_is_file_like = all(hasattr(attr) for attr in ('seek', 'close', 'read', 'write'))
    compress = False
    if out_file[-3:] == ".gz":
        out_file = out_file[:-3]
        compress = True

    with open(out_file, "wb") as output_stream:
        if is_gz_file(in_file):
            p = Popen(["bash", "zcat", in_file], stdout=PIPE)
            input_stream = p.stdout
        else:
            p = Popen(["bash", "cat", in_file], stdout=PIPE)
            input_stream = p.stdout

        #p = annotate_stream_vep(input_stream, output_stream)

        #pread, pwrite = os.pipe()
        #read_pipe = os.fdopen(pread)
        #write_pipe = os.fdopen(pwrite)
        #os.set_blocking(out_pipe.fileno(), False)
        #os.set_blocking(out_pipe.fileno(), False)
        #other_pipe = os.fdopen(w)
        # p = annotate_stream_snpeff(input_stream, output_stream, clinvar_fields=clinvar_fields,
        p = annotate_stream_vep(input_stream, output_stream)
        """
        while True:
            print(p.poll())
            if not p.poll() is None:
                print("The process is closed")
            else:
                print("The process is not closed")
            #line = in_pipe.read()
            line = read(read_pipe.fileno(), 1024)
            print(line)
            if not line:
                print("Break")
                break
            else:
                print("Previous to write")
                output_stream.write(line)
        """
        #for line in in_pipe:
        #    print(line, end="")
        #    output_stream.write(line)
        #while True:
        #    print("hola")
        #    line = current_pipe.readline()
        #    print("adeu")
        #    if not line:
        #        break
        #    output_stream.write(line)
        #for line in other_pipe:
        #    print("linea")
        #    print(line)
        #    output_stream.write(line)
        output_stream.close()
        #other_pipe.close()
        #output_stream.write(other_pipe.read())
        # open('image.jpg', 'wb').write(p.body_file().read())
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
