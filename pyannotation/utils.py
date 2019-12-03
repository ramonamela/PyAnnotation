#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import json
from subprocess import Popen, PIPE
import sys
from jsonschema import validators, Draft7Validator
import gzip
import binascii
from itertools import islice
import multiprocessing
import pandas as pd
import os
from os.path import dirname
from os.path import abspath


def bgzip_tabix_file(in_file, out_file=None):
    with open(in_file, 'rb', 0) as f:
        bgzip_tabix_stream(f, in_file=in_file, out_file=out_file)


def bgzip_tabix_stream(in_stream, in_file=None, out_file=None):
    if in_file is None and out_file is None:
        raise Exception("Could not compute the output name in the function \"bgzip_tabix_stream\"")
    if out_file is None:
        out_stream = open(in_file + ".gz", 'w')
        tabix_file = in_file + ".gz"
    else:
        out_stream = open(out_file, "w")
        tabix_file = out_file
    p = Popen(["bash", "-c", "bgzip -c "], stdin=in_stream, stdout=out_stream)
    p.wait()
    p = Popen(["bash", "-c", "tabix -p vcf " + tabix_file], stdout=open(os.devnull, 'w'))
    p.wait()


def is_gz_file(filepath):
    with open(filepath, 'rb') as test_f:
        return binascii.hexlify(test_f.read(2)) == b'1f8b'


def open_potential_gz(input_vcf):
    if (is_gz_file(input_vcf)):
        return gzip.open(input_vcf, "r")
    else:
        return open(input_vcf, "r")


def from_info_field_to_dataframe(info_field):
    splitted_info = list(map(lambda x: x.split("="), info_field.split(";")))
    info_field_headers = list(map(lambda x: x[0], splitted_info))
    info_field_content = list(map(lambda x: x[1], splitted_info))
    return pd.DataFrame([info_field_content], columns=info_field_headers)


def compute_hgvs_id(vcf_dataframe, allele_separator="&"):
    vcf_header = vcf_dataframe.columns.tolist()
    format_index = vcf_header.index("FORMAT")
    # VEP CASE
    HGVSc_index = vcf_header.index("HGVSc") if "HGVSc" in vcf_header else None
    # SnpEff CASE
    HGVS_c_index = vcf_header.index("HGVS.c") if "HGVS.c" in vcf_header else None
    Feature_ID_index = vcf_header.index("Feature_ID") if "Feature_ID" in vcf_header else None
    if not HGVSc_index is None:
        ## VEP
        splitted_hgvs_colon = vcf_dataframe["HGVSc"].apply(lambda x: allele_separator.join(list(
            map(lambda y: "." if y == "." else ":".join([((y.split(":")[0]).split("."))[0], y.split(":")[1]]),
                x.split(allele_separator)))))
        vcf_dataframe["hgvsc_id"] = splitted_hgvs_colon
    elif not HGVS_c_index is None and not Feature_ID_index is None:
        ## SnpEff
        transcript = list(
            map(lambda y: list(map(lambda x: x.split(".")[0] if x[0:4] == "ENST" else "", y.split(allele_separator))),
                vcf_dataframe["Feature_ID"]))
        rna_position = list(map(lambda x: x.split(allele_separator), list(vcf_dataframe["HGVS.c"])))
        vcf_dataframe["hgvsc_id"] = list(map(lambda z: allele_separator.join(z), list(
            map(lambda x: list(map(lambda y: ":".join(list(y)), zip(x[0], x[1]))),
                list(zip(transcript, rna_position))))))
    vcf_header = vcf_header[:format_index] + ["hgvsc_id"] + vcf_header[format_index:]
    vcf_dataframe = vcf_dataframe[vcf_header]
    return vcf_dataframe


def split_info(vcf_dataframe):
    vcf_header = vcf_dataframe.columns.tolist()
    info_index = vcf_header.index("INFO")

    new_info_dataframe = pd.concat(list(map(from_info_field_to_dataframe, list(vcf_dataframe["INFO"]))),
                                   sort=False).reset_index(drop=True)
    new_info_headers = new_info_dataframe.columns.tolist()

    vcf_dataframe_without_info = vcf_dataframe.drop(columns="INFO")
    result_dataframe = pd.concat([vcf_dataframe_without_info, new_info_dataframe], axis=1)

    new_header_order = vcf_header[:info_index] + new_info_headers + vcf_header[info_index + 1:]
    result_dataframe = result_dataframe[new_header_order]
    return result_dataframe


def split_ann(vcf_dataframe, ann_header):
    vcf_header = vcf_dataframe.columns.tolist()
    ann_index = vcf_header.index("ANN")

    lists_with_ann_content = list(map(lambda z: [list(a) for a in zip(*z)], list(
        map(lambda x: list(map(lambda y: list(map(lambda w: "." if w == "" else w, y.split("|"))), x.split(","))),
            list(vcf_dataframe["ANN"])))))
    new_ann_dataframe = pd.DataFrame(list(map(lambda x: list(map(lambda y: ','.join(y), x)), lists_with_ann_content)),
                                     columns=ann_header)

    vcf_dataframe_without_ann = vcf_dataframe.drop(columns="ANN")
    result_dataframe = pd.concat([vcf_dataframe_without_ann, new_ann_dataframe], axis=1)

    new_header_order = vcf_header[:ann_index] + ann_header + vcf_header[ann_index + 1:]
    result_dataframe = result_dataframe[new_header_order]
    return result_dataframe


def get_multicolumn_by_id(dataframe, id_dataframe, id_column, columns_to_keep, split_separator="&", join_separator="&"):
    # Since we search for row id, we assume there is a single row at max
    def handle_field(element):
        current_row = dataframe.loc[dataframe[id_column] == element][columns_to_keep].values.tolist()
        return current_row[0] if len(current_row) == 1 else ["."] * len(columns_to_keep)

    col_values = list(
        map(lambda row: list(
            map(lambda x: split_separator.join(x), list(zip(*list(map(handle_field, row.split(join_separator))))))),
            id_dataframe))
    multicolumn_dataframe = pd.DataFrame(col_values, columns=columns_to_keep)
    return multicolumn_dataframe


def add_civic_variant_summaries_fields(vcf_dataframe, path_to_tsv, field_list):
    with open(path_to_tsv) as f:
        header = f.readline().rstrip("\n").split("\t")
        variant_summaries_fields_list = list(map(lambda y: (y[:len(header) - 1] + [",".join(y[len(header) - 1:])]),
                                                 list(map(lambda x: x.rstrip("\n").split("\t"), list(f)))))
    pandas_variant_summaries = pd.DataFrame(variant_summaries_fields_list, columns=header)
    pandas_variant_summaries["hgvsc_id"] = list(
        map(lambda x: next(
            (t.split(":")[0].split(".")[0] + ":" + t.split(":")[-1] for t in x.split(",") if (t[0:4] == 'ENST')), None),
            pandas_variant_summaries["hgvs_expressions"]))

    fields_to_merge = [x for x in field_list if x[0] in pandas_variant_summaries.columns]
    field_new_names = [x[1] for x in fields_to_merge]
    fields_to_merge.extend([("hgvsc_id", "hgvsc_id")])
    field_old_names = [x[0] for x in fields_to_merge]
    dataframe_to_merge = pandas_variant_summaries[field_old_names]
    multicolumn_dataframe = get_multicolumn_by_id(dataframe_to_merge, vcf_dataframe["hgvsc_id"], "hgvsc_id",
                                                  field_old_names, split_separator="&", join_separator="&")
    del multicolumn_dataframe["hgvsc_id"]
    multicolumn_dataframe.columns = field_new_names
    vcf_header = list(vcf_dataframe.columns)
    format_index = vcf_header.index("FORMAT")
    vcf_dataframe = pd.concat([vcf_dataframe, multicolumn_dataframe], axis=1, sort=False)
    vcf_dataframe = vcf_dataframe[vcf_header[:format_index] + field_new_names + vcf_header[format_index:]]
    return vcf_dataframe


def add_columns_from_several_ids(input_dataframe, auxiliar_dataframe, field_list, id_fields, input_separator="&"):
    ## ID fields (input, auxiliar)
    auxiliar_fields_to_select = [x[0] for x in field_list if x[0] in list(auxiliar_dataframe.columns)]
    new_auxiliar_fields_to_select = [x[1] for x in field_list if x[0] in list(auxiliar_dataframe.columns)]
    origin_fields_id = [x[0] for x in id_fields]

    def f(x):
        current_candidates = list(zip(*list(map(lambda x: x.split(input_separator), x))))
        solution_dictionary = {}
        for field in auxiliar_fields_to_select:
            solution_dictionary[field] = []
        for candidate in current_candidates:
            base_mask = (auxiliar_dataframe[id_fields[0][1]] == candidate[0])
            for i in range(1, len(id_fields)):
                base_mask = base_mask & (auxiliar_dataframe[id_fields[i][1]] == candidate[i])
            dataframe_to_return = auxiliar_dataframe[base_mask]
            for field in auxiliar_fields_to_select:
                if dataframe_to_return.shape[0] > 0:
                    solution_dictionary[field].append("|".join(list(dataframe_to_return[field])))
                else:
                    solution_dictionary[field].append(".")
        list_of_values = []
        for key in solution_dictionary.keys():
            list_of_values.append("&".join(solution_dictionary[key]))
        values = list_of_values
        return "\t".join(values)

    result = list(map(lambda x: x.split("\t"), input_dataframe[origin_fields_id].apply(f, axis=1)))
    new_values_dataframe = pd.DataFrame(result, columns=new_auxiliar_fields_to_select)
    input_dataframe = pd.concat([input_dataframe, new_values_dataframe], axis=1, sort=False)
    return input_dataframe


def add_civic_variant_clinical_evidence_fields(vcf_dataframe, path_to_tsv, field_list):
    variant_clinical_evidence = pd.read_csv(path_to_tsv, sep="\t")
    id_vector = [("gene", "gene"), ("variant", "variant")]
    auxiliar_fields = field_list + id_vector
    columns_old_names = [x[0] for x in auxiliar_fields if x[0] in list(variant_clinical_evidence.columns)]
    columns_new_names = [x[1] for x in auxiliar_fields if x[0] in list(variant_clinical_evidence.columns)]
    variant_clinical_evidence = variant_clinical_evidence[columns_old_names]
    vcf_dataframe = add_columns_from_several_ids(vcf_dataframe, variant_clinical_evidence, field_list, id_vector)
    vcf_dataframe[columns_old_names].columns = columns_new_names
    return vcf_dataframe


def collapse_dataframe_column(input_dataframe, single_cols=[], multiple_cols=[], external_sep=",",
                              old_internal_sep="&", new_internal_sep="&"):
    for pattern in single_cols:
        if pattern in list(input_dataframe.columns):
            input_dataframe[pattern][pd.notnull(input_dataframe[pattern])] = input_dataframe[pattern][
                pd.notnull(input_dataframe[pattern])].apply(
                lambda x: x.split(external_sep)[0].replace(old_internal_sep, new_internal_sep))
    for pattern in multiple_cols:
        if pattern in list(input_dataframe.columns):
            input_dataframe[pattern] = input_dataframe[pattern].apply(
                lambda x: x.replace(external_sep, new_internal_sep) if x.split(external_sep)[0] == "." else
                x.split(external_sep)[0].replace(old_internal_sep, new_internal_sep))
    return input_dataframe


def clean_dataframe_column(input_dataframe, cols=[], old_allele_sep=",", new_allele_sep=",",
                           old_internal_sep="&", new_internal_sep="&"):
    for pattern in cols:
        if pattern in list(input_dataframe.columns):
            input_dataframe[pattern][pd.notnull(input_dataframe[pattern])] = input_dataframe[pattern][
                pd.notnull(input_dataframe[pattern])].apply(
                lambda x: new_allele_sep.join(
                    list(map(lambda y: y.replace(old_internal_sep, new_internal_sep), x.split(old_allele_sep)))))
    return input_dataframe


def split_dataframe_column(input_dataframe, cols=[], old_internal_sep=",", new_internal_sep="&",
                           reference_column="Allele", reference_separator="&"):
    column_length = list(map(lambda x: len(x.split(reference_separator)), list(input_dataframe[reference_column])))
    for pattern in cols:
        if pattern in list(input_dataframe.columns):
            input_dataframe[pattern] = list(map(
                lambda x: new_internal_sep.join(["."] * x[0]) if x[1] == "." else x[1].replace(old_internal_sep,
                                                                                               new_internal_sep),
                list(zip(column_length, input_dataframe[pattern]))))
    return input_dataframe


def add_civic_gene_description_fields(input_dataframe, gene_summaries_tsv, gene_summaries_civic_fields):
    with open(gene_summaries_tsv) as f:
        header = f.readline().rstrip("\n").split("\t")
        gene_summaries_fields_list = list(map(lambda y: (y[:len(header) - 1] + [",".join(y[len(header) - 1:])]),
                                              list(map(lambda x: x.rstrip("\n").split("\t"), list(f)))))
    pandas_gene_summaries = pd.DataFrame(gene_summaries_fields_list, columns=header)
    id_vector = [("gene", "name")]
    new_columns_old_names = [x[0] for x in gene_summaries_civic_fields if x[0] in list(pandas_gene_summaries.columns)]
    new_columns_new_names = [x[0] for x in gene_summaries_civic_fields if x[1] in list(pandas_gene_summaries.columns)]
    auxiliar_fields = gene_summaries_civic_fields + id_vector
    columns_old_names = [x[0] for x in auxiliar_fields if x[0] in list(pandas_gene_summaries.columns)] + [x[1] for x in
                                                                                                          id_vector if
                                                                                                          x[1] in list(
                                                                                                              pandas_gene_summaries.columns)]
    variant_clinical_evidence = pandas_gene_summaries[columns_old_names]
    vcf_dataframe = add_columns_from_several_ids(input_dataframe, variant_clinical_evidence,
                                                 gene_summaries_civic_fields, id_vector)
    return vcf_dataframe


def treat_chunk(chunk, vcf_header, ann_header,
                input_additional_multifields=None,
                vep_fields=None,
                snpeff_fields=None,
                clinvar_prefix=None,
                clinvar_fields=None,
                dbNSFP_fields=None):
    if input_additional_multifields is None:
        input_additional_multifields = ["platformnames", "datasetnames", "callsetnames", "filt"]
    if vep_fields is None:
        vep_fields = ["Allele", "Consequence", "IMPACT", "SYMBOL", "Gene", "Feature_type", "Feature", "BIOTYPE" "EXON",
                      "INTRON", "HGVSc", "HGVSp", "cDNA_position", "CDS_position", "Protein_position", "Amino_acids",
                      "Codons", "Existing_variation", "DISTANCE", "STRAND", "FLAGS", "SYMBOL_SOURCE", "HGNC_ID",
                      "SOURCE", "HGVS_OFFSET", "HGVSg"]
    if snpeff_fields is None:
        snpeff_fields = ["Allele", "Annotation", "Annotation_Impact", "Gene_Name", "Gene_ID", "Feature_Type",
                         "Feature_ID", "Transcript_BioType", "Rank", "HGVS.c", "HGVS.p", "cDNA.pos/cDNA.length",
                         "CDS.pos/CDS.length", "AA.pos/AA.length", "Distance", "ERRORS/WARNINGS/INFO"]
    if clinvar_fields is None:
        clinvar_fields = ["CLNDISDB", "CLNDISDBINCL", "CLNDN", "CLNSIG", "CLNREVSTAT"]
    if dbNSFP_fields is None:
        dbNSFP_fields = ["MutationAssessor_score", "MutationAssessor_pred", "SIFT_score", "SIFT_pred",
                         "FATHMM_score", "FATHMM_pred", "Polyphen2_HDIV_score", "Polyphen2_HDIV_rankscore",
                         "Polyphen2_HDIV_pred", "Polyphen2_HVAR_score", "Polyphen2_HVAR_rankscore",
                         "Polyphen2_HVAR_pred"]
    if clinvar_prefix is None:
        clinvar_prefix = "ClinVar_ID"

    format_index = vcf_header.index("FORMAT")
    vcf_fields_header = vcf_header[:format_index]

    info_index = vcf_header.index("INFO")
    vcf_fields_header.pop(info_index)

    input_vcf_dataframe = pd.DataFrame(list(map(lambda x: x.decode("utf-8")[:-1].replace(" ", "").split("\t"), chunk)),
                                       columns=vcf_header)
    vcf_dataframe_with_info = split_info(input_vcf_dataframe)
    vcf_dataframe_with_info_ann = split_ann(vcf_dataframe_with_info, ann_header)
    vcf_dataframe_with_info_corrected = clean_dataframe_column(vcf_dataframe_with_info_ann,
                                                               cols=input_additional_multifields,
                                                               old_allele_sep="|", new_allele_sep="|",
                                                               old_internal_sep=",", new_internal_sep="|")
    vcf_dataframe_with_info_corrected = clean_dataframe_column(vcf_dataframe_with_info_corrected,
                                                               cols=vep_fields + snpeff_fields,
                                                               old_allele_sep=",", new_allele_sep="&",
                                                               old_internal_sep="&", new_internal_sep="|")
    clinvar_vep_result = [clinvar_prefix] + [clinvar_prefix + "_" + x for x in clinvar_fields]
    vcf_dataframe_with_vep_corrected = collapse_dataframe_column(vcf_dataframe_with_info_corrected,
                                                                 single_cols=clinvar_vep_result,
                                                                 multiple_cols=dbNSFP_fields, external_sep=",",
                                                                 old_internal_sep="&", new_internal_sep="&")
    vcf_dataframe_with_vep_corrected = clean_dataframe_column(vcf_dataframe_with_vep_corrected,
                                                              cols=clinvar_vep_result,
                                                              old_allele_sep=",", new_allele_sep="&",
                                                              old_internal_sep="&", new_internal_sep="|")
    dbNSFP_snpsift_result = ["dbNSFP_" + x for x in dbNSFP_fields]
    vcf_dataframe_with_splitted_cols = split_dataframe_column(vcf_dataframe_with_vep_corrected,
                                                              cols=dbNSFP_snpsift_result, old_internal_sep=",",
                                                              new_internal_sep="&", reference_column="Allele",
                                                              reference_separator="&")
    vcf_dataframe_with_hgvs_id = compute_hgvs_id(vcf_dataframe_with_splitted_cols)
    module_path = dirname(dirname(abspath(__file__)))
    variant_summaries_fields_tsv = module_path + "/data/cache/civic/nightly-VariantSummaries.tsv"
    vcf_dataframe_with_civic_variant_summaries = add_civic_variant_summaries_fields(vcf_dataframe_with_hgvs_id,
                                                                                    variant_summaries_fields_tsv,
                                                                                    [("random", "random"),
                                                                                     ("gene", "gene"),
                                                                                     ("variant", "variant"), (
                                                                                         "civic_variant_evidence_score",
                                                                                         "civic_variant_evidence_score")])
    variant_clinical_evidence_tsv = module_path + "/data/cache/civic/nightly-ClinicalEvidenceSummaries.tsv"
    variant_clinical_evidence_civic_fields = [("random", "random"),
                                              ("drugs", "drugs"),
                                              ("variant_summary", "variant_summary"),
                                              ("disease", "disease"),
                                              ("evidence_type", "evidence_type"),
                                              ("evidence_direction", "evidence_direction"),
                                              ("evidence_level", "evidence_level"),
                                              ("clinical_significance", "clinical_significance"),
                                              ("evidence_statement", "evidence_statement"),
                                              ("citation_id", "citation_id")]
    vcf_dataframe_with_civic_clinical_evidence = add_civic_variant_clinical_evidence_fields(
        vcf_dataframe_with_civic_variant_summaries,
        variant_clinical_evidence_tsv,
        variant_clinical_evidence_civic_fields)
    gene_summaries_tsv = module_path + "/data/cache/civic/nightly-GeneSummaries.tsv"
    gene_summaries_civic_fields = [("random", "random"), ("description", "gene_description")]
    vcf_dataframe_with_civic_gene_description = add_civic_gene_description_fields(
        vcf_dataframe_with_civic_clinical_evidence, gene_summaries_tsv, gene_summaries_civic_fields)
    return vcf_dataframe_with_civic_gene_description


def normalize_vcf_to_stream(path_to_vcf):
    p = Popen(["bash", "-c", "bcftools norm -m -both " + path_to_vcf], stdout=PIPE)
    return p.stdout


def from_vcf_to_maf(input_vcf, output_maf, chunk_size=10):
    """Converts vcf to maf.
    Parameters
    ----------
    input_vcf: path to the input vcf (can be gzipped) but the format is not checked
    output_maf: path where the MAF should be stored
    """

    with normalize_vcf_to_stream(input_vcf) as file:

        for line in file:
            line_string = line.decode("utf-8")
            if line_string.startswith("##INFO=<ID=ANN"):

                def detect_description(part):
                    splitted_part = part.split("=")
                    if splitted_part[0] == "Description":
                        return splitted_part[1].split(":")[1].replace(" ", "")
                    return None

                ann_header = list(filter(None, map(detect_description, list(
                    line_string[7:-1].replace("\"", "").replace("<", "").replace(">", "").split(",")))))[0].replace("'",
                                                                                                                    "").split(
                    "|")
            if line_string.startswith("#CHROM"):
                break

        vcf_header = line.decode("utf-8")[1:-1].split("\t")
        format_index = vcf_header.index("FORMAT")
        individual_headers = vcf_header[format_index:]

        mp_pool = multiprocessing.Pool()
        results = []

        for chunk in iter(lambda: list(islice(file, chunk_size)), []):
            results.append(mp_pool.apply_async(treat_chunk, (chunk, vcf_header, ann_header)))

        results = list(map(lambda x: x.get(), results))

        big_dataframe = pd.concat(results, sort=False)
        all_headers = list(big_dataframe.columns)
        common_fields = [x for x in all_headers if not x in individual_headers]
        new_order = common_fields + individual_headers
        big_dataframe = big_dataframe[new_order]
        big_dataframe.to_csv(output_maf, sep='\t', na_rep='.', index=False)


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
    return "java -Xmx64g -jar $(readlink -f $(which SnpSift.jar)) annotate -id -noInfo " \
           "/home/ramela/Downloads/example_vcfs/databases/00-All.vcf.gz /dev/stdin -a"


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
        validator_class, {"properties": set_defaults},
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
    null_out = open(os.devnull, 'w')
    p.wait()
    p = Popen(["bash", "-c", "tabix", "-f", output_file], stdout=null_out)
    p.wait()


def main(in_file, out_file, json_file):
    return annotate(in_file, out_file, json_file)


if __name__ == "__main__":
    in_file = sys.argv[1]
    out_file = sys.argv[2]
    # json_file = sys.argv[3]
    from_vcf_to_maf(in_file, out_file)
    # main(in_file, out_file, json_file)
