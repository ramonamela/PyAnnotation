#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse
import textwrap
import sys
import json
from os.path import dirname
from os.path import abspath
from pyannotation.utils import open_potential_gz
from pyannotation.utils import annotate_stream_to_stream

if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=textwrap.dedent("""\
    snpEff and VEP can be used as effect predictors.
    Furthermore, dbNSFP, dbSNP, CIViC and COSMIC databases are available to annotate the predicted variants.
    """), epilog=textwrap.dedent("""\
    When the input file is not given, input is expected to be given through the standard input.
    
    When the output file is not given, the output is given through the standard output.
    
    If the json file is present, all the flags corresponding to the fields to annotate from the different databases 
    are ignored.
    """))
    parser.add_argument("-i", "--input_file", help="path to input vcf file")
    parser.add_argument("-o", "--output_file", help="path to output maf file")
    parser.add_argument("-j", "--json_file", help="path to the json file descriving the actions to be performed")
    parser.add_argument("-t", "--prediction_tool", choices=["vep", "snpEff"], default="", help="variant predictor tool")
    parser.add_argument("--dbSNP", help="comma separated fields from dbSNP added to the maf output", default="")
    parser.add_argument("--dbNSFP", help="comma separated fields from dbNSFP added to the maf output", default="")
    parser.add_argument("-dbCGS", "--dbCiVICGeneSummaries",
                        help="comma separated fields from the Gene Summaries CiVIC file added to the maf output",
                        default="")
    parser.add_argument("-dbCVS", "--dbCiVICVariantSummaries",
                        help="comma separated fields from the Variant Summaries CiVIC file added to the maf output",
                        default="")
    parser.add_argument("-dbCCE", "--dbCiVICClinicalEvidence",
                        help="comma separated fields from the Clinical Evidence CiVIC file added to the maf output",
                        default="")
    parser.add_argument("-dbCV", "--dbClinVar", help="comma separated fields from ClinVar added to the maf output",
                        default="")
    parser.add_argument("-dbCCM", "--dbCOSMICCodingMutations",
                        help="comma separated fields from the Coding Mutations COSMIC file added to the maf output",
                        default="")
    parser.add_argument("-dbNCV", "--dbCOSMICNonCodingVariants",
                        help="comma separated fields from the Non Coging Variants COSMIC file added to the maf output",
                        default="")
    parser.add_argument("-dbCMC", "--dbCOSMICMutantCensus",
                        help="comma separated fields from the Mutation Census COSMIC file added to the maf output",
                        default="")

    args = parser.parse_args()
    input_file = args.input_file
    if input_file is None:
        input_stream = sys.stdin
    else:
        input_stream = open_potential_gz(input_file)

    output_file = args.output_file
    if output_file is None:
        output_stream = sys.stdout
    else:
        output_stream = open(output_file, "w")

    module_path = dirname(dirname(abspath(__file__)))

    with open(module_path + "/PyAnnotation/config/schema.json") as f:
        schema_json = json.load(f)

    if args.json_file is None:
        ## Args handling -> build JSON so we can verify the format
        config_json = {}
        input_json = vars(args)
        config_json["prediction_tool"] = [] if input_json["prediction_tool"] == "" else input_json["prediction_tool"]

        array_input_parameters = ["dbSNP", "dbNSFP", "dbCiVICGeneSummaries",
                                  "dbCiVICVariantSummaries", "dbCiVICClinicalEvidence", "dbClinVar",
                                  "dbCOSMICCodingMutations", "dbCOSMICNonCodingVariants", "dbCOSMICMutantCensus"]

        config_json["annotated_fields"] = []

        for field in array_input_parameters:
            if input_json[field] == "" or input_json[field] is None:
                config_json["annotated_fields"].append({"database": field, "fields": []})
            elif input_json[field] == "all":
                config_json["annotated_fields"].append({"database": field, "fields": schema_json["definitions"][field + "_fields"]["enum"]})
            else:
                config_json["annotated_fields"].append({"database": field, "fields": input_json[field].split(",")})
    else:
        with open(args.json_file, 'r') as f:
            config_json = json.load(f)

    config_json.pop("input_field", None)
    config_json.pop("output_file", None)

    from pyannotation.utils import _check_json_file
    try:
        _check_json_file(config_json, schema_json)
    except Exception as e:
        print(e)
    annotate_stream_to_stream(input_stream, output_stream, config_json)
