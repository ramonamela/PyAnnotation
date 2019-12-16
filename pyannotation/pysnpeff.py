#!/usr/bin/env python3
import sys
from pyannotation.utils import open_potential_gz
from os.path import dirname
from os.path import abspath
from subprocess import Popen, PIPE
from pyannotation.utils import bgzip_tabix_file
from pyannotation.utils import is_gz_file
import io


def annotate_stream_snpeff(input_stream, output_stream, clinvar_fields=[], dbNSFP_fields=[], dbSNP_fields=[]):
    module_path = dirname(dirname(abspath(__file__)))

    # Reference genome
    command_string = "java -Xmx64g -jar " + module_path + "/dependencies/snpEff/snpEff.jar -t GRCh38.86"
    if len(dbNSFP_fields) > 0 or len(clinvar_fields) > 0 or len(dbNSFP_fields) > 0:
        p = Popen(["bash", "-c", command_string], stdin=input_stream, stdout=PIPE)
        p_out = p.stdout
    else:
        p = Popen(["bash", "-c", command_string], stdin=input_stream, stdout=output_stream)

    # dbNSFP
    if len(dbNSFP_fields) > 0:
        command_string = "java -Xmx64g -jar " + module_path + "/dependencies/snpEff/SnpSift.jar dbnsfp -f '" + ",".join(
            dbNSFP_fields) + "' -db " + module_path + "/data/cache/dbNSFP/dbNSFP4.0a.txt.gz" + " -m -a /dev/stdin"
        if len(clinvar_fields) > 0:
            p = Popen(["bash", "-c", command_string], stdin=p_out, stdout=PIPE)
            p_out = p.stdout
        else:
            p = Popen(["bash", "-c", command_string], stdin=p_out, stdout=output_stream)

    # dbSNP
    if len(dbSNP_fields) > 0:
        command_string = "java -Xmx64g -jar " + module_path + "/dependencies/snpEff/SnpSift.jar annotate -info '" + \
                         ",".join(
                             dbSNP_fields) + "' " + module_path + "/data/cache/dbSNP/All_20180418.vcf.gz /dev/stdin -a"
        if len(clinvar_fields) > 0:
            p = Popen(["bash", "-c", command_string], stdin=p_out, stdout=PIPE)
            p_out = p.stdout
        else:
            p = Popen(["bash", "-c", command_string], stdin=p_out, stdout=output_stream)

    # Clinvar
    if len(clinvar_fields) > 0:
        command_string = "java -Xmx64g -jar " + module_path + "/dependencies/snpEff/SnpSift.jar annotate -info '" + \
                         ",".join(
                             clinvar_fields) + "' " + module_path + "/data/cache/clinvar/clinvar.vcf.gz /dev/stdin -a"
        p = Popen(["bash", "-c", command_string], stdin=p_out, stdout=output_stream)
    return p


def annotate_snpeff(in_file, out_file, clinvar_fields=[], dbNSFP_fields=[]):
    compress = False
    if out_file[-3:] == ".gz":
        out_file = out_file[:-3]
        compress = True
    # with open_potential_gz(input_file) as input_stream, open(out_file, "w") as output_stream:
    with open(out_file, "w") as output_stream:
        if is_gz_file(in_file):
            p = Popen(["bash", "zcat", in_file], stdout=PIPE)
            input_stream = p.stdout
        else:
            p = Popen(["bash", "cat", in_file], stdout=PIPE)
            input_stream = p.stdout
        valid_stream = io.BufferedReader(io.BytesIO())
        #p = annotate_stream_snpeff(input_stream, output_stream, clinvar_fields=clinvar_fields,
        p = annotate_stream_snpeff(input_stream, valid_stream, clinvar_fields=clinvar_fields, dbNSFP_fields=dbNSFP_fields)
        with open('image.jpg', 'wb') as out_file:
            out_file.write(valid_stream.read())
        #open('image.jpg', 'wb').write(p.body_file().read())
        p.wait()
    if compress:
        bgzip_tabix_file(out_file, out_file=out_file + ".gz")


if __name__ == "__main__":
    in_file = sys.argv[1]
    out_file = sys.argv[2]
    annotate_snpeff(in_file, out_file)
