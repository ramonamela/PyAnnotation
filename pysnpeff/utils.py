#!/usr/bin/env python3
import json
from subprocess import Popen, PIPE
import sys

def _from_json_to_command(config):
    """Annotates the input_stream to the output_stream following the config
    directives.
    Parameters
    ----------
    input_stream : TextIOBase
        Input stream.
    output_stream : TextIOBase
        Output stream.
    config : dict
        Dictionary describing the annotation step to be performed
    """
    return "java -Xmx64g -jar /mnt/nfs/eucancan_vcf_vc/SnpSift/target/SnpSift-4.4.jar annotate -id -noInfo /mnt/nfs/eucancan_vcf_vc/databases/00-All.vcf.gz /dev/stdin -a"

def annotate(input_file, output_file, json_file):
    """Annotates the input_file to the output_file following the configuration
    present in json_file.
    Parameters
    ----------
    input_stream : TextIOBase
        Input stream.
    output_stream : TextIOBase
        Input stream.
    config : dict
        Dictionary containing the necessary information to perform the
        annotation
    """
    config_json = json.load(json_file)
    out_stream = open(input_file, "r")
    for i in range(config_json["steps"]):
        command = str.split(_from_json_to_command(config_json["steps"][i]))
        in_stream = out_stream
        if i == len(config_json["steps"] - 1):
            out_stream = open(output_file, "rw")
        else:
            out_stream = PIPE
        p = Popen(command, stdin=in_stream, stdout=out_stream)
        out_stream = p.stdout
    p.wait()

def main(in_file, out_file, json_file):
    annotate(in_file, out_file, json_file)

if __name__ == "__main__":
    in_file = sys.argv[1]
    out_file = sys.argv[2]
    json_file = sys.argv[3]
    main(in_file, out_file, json_file)