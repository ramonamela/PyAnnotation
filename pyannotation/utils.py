#!/usr/bin/env python3
import json
from subprocess import Popen, PIPE
import sys
from jsonschema import validate, validators, Draft7Validator

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
    return "java -Xmx64g -jar $(readlink -f $(which SnpSift.jar)) annotate -id -noInfo /home/ramela/Downloads/example_vcfs/databases/00-All.vcf.gz /dev/stdin -a"

def _validate_with_default(validator_class):
    validate_properties = validator_class.VALIDATORS["properties"]

    def set_defaults(validator, properties, instance, schema):
        for property, subschema in properties.items():
            if "default" in subschema:
                instance.setdefault(property, subschema["default"])

        for error in validate_properties(
            validator, properties, instance, schema,
        ):
            yield error

    return validators.extend(
        validator_class, {"properties" : set_defaults},
    )

def _check_json_file(input_json, json_schema):
    DefaultValidatingDraft7Validator = _validate_with_default(Draft7Validator)
    DefaultValidatingDraft7Validator(json_schema).validate(input_json)

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

    with open(json_file, 'r') as f:
        config_json = json.load(f)

    with open("/home/ramela/git/PySnpEff/pysnpeff/config/schema.json") as f:
        schema_json = json.load(f)

    _check_json_file(config_json, schema_json)

    p = Popen(["bash", "-c", "zcat " + input_file], stdout=PIPE)
    out_stream = p.stdout
    for i in range(len(config_json["steps"])):
        command = ["bash", "-c", _from_json_to_command(config_json["steps"][i])]
        in_stream = out_stream
        p = Popen(command, stdin=in_stream, stdout=PIPE)
        out_stream = p.stdout
    in_stream = out_stream
    out_stream = open(output_file, "w")
    p = Popen(["bash", "-c", "bgzip"], stdin=in_stream, stdout=out_stream)
    p.wait()

def main(in_file, out_file, json_file):
    return annotate(in_file, out_file, json_file)

if __name__ == "__main__":
    in_file = sys.argv[1]
    out_file = sys.argv[2]
    json_file = sys.argv[3]
    main(in_file, out_file, json_file)