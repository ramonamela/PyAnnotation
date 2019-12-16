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
import time
import io
import numpy as np


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
    info_field_content = list(map(lambda x: x[1] if len(x) > 1 else "true", splitted_info))
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
    else:
        return vcf_dataframe
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


def split_ann(vcf_dataframe, ann_header, allele_separator=","):
    vcf_header = vcf_dataframe.columns.tolist()
    ann_index = vcf_header.index("ANN")

    lists_with_ann_content = list(map(lambda z: [list(a) for a in zip(*z)], list(
        map(lambda x: list(map(lambda y: list(map(lambda w: "." if w == "" else w, y.split("|"))), x.split(","))),
            list(vcf_dataframe["ANN"])))))
    new_ann_dataframe = pd.DataFrame(
        list(map(lambda x: list(map(lambda y: allele_separator.join(y), x)), lists_with_ann_content)),
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
                    dataframe_to_return = dataframe_to_return.fillna(value=".")
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
    for pattern in set(single_cols):
        if pattern in list(input_dataframe.columns):
            input_dataframe[pattern][pd.notnull(input_dataframe[pattern])] = input_dataframe[pattern][
                pd.notnull(input_dataframe[pattern])].apply(
                lambda x: x.split(external_sep)[0].replace(old_internal_sep, new_internal_sep))
    for pattern in set(multiple_cols):
        if pattern in list(input_dataframe.columns):
            input_dataframe[pattern] = input_dataframe[pattern].apply(
                lambda x: x.replace(external_sep, new_internal_sep) if x.split(external_sep)[0] == "." else
                x.split(external_sep)[0].replace(old_internal_sep, new_internal_sep))
    return input_dataframe


def clean_dataframe_column(input_dataframe, cols=[], old_allele_sep=",", new_allele_sep=",",
                           old_internal_sep="&", new_internal_sep="&"):
    for pattern in set(cols):
        if pattern in list(input_dataframe.columns):
            input_dataframe[pattern][pd.notnull(input_dataframe[pattern])] = input_dataframe[pattern][
                pd.notnull(input_dataframe[pattern])].apply(
                lambda x: new_allele_sep.join(
                    list(map(lambda y: y.replace(old_internal_sep, new_internal_sep), x.split(old_allele_sep)))))
    return input_dataframe


def split_dataframe_column(input_dataframe, cols=[], old_internal_sep=",", new_internal_sep="&",
                           reference_column="Allele", reference_separator="&"):
    column_length = list(map(lambda x: len(x.split(reference_separator)), list(input_dataframe[reference_column])))
    for pattern in set(cols):
        if pattern in list(input_dataframe.columns):
            input_dataframe[pattern] = list(map(
                lambda x: new_internal_sep.join(["."] * x[0]) if (str(x[1]) == 'nan') or (x[1] == ".") else x[
                    1].replace(old_internal_sep,
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


def add_dbsnp_fields(input_dataframe, path_to_dbsnp, dbsnp_multialt_fields=[], dbsnp_unique_fields=[],
                     dbsnp_boolean_fields=[]):
    current_values = input_dataframe[["CHROM", "POS", "REF", "ALT", "Allele"]]
    import tabix
    dbsnp = tabix.open(path_to_dbsnp)

    def get_columns(x):
        try:
            current_dataframe = pd.DataFrame(
                dbsnp.query(str(x["CHROM"]).replace("chr", ""), int(x["POS"]) - 1, int(x["POS"])),
                columns=["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"])
            current_row = current_dataframe[
                (current_dataframe["POS"] == x["POS"]) & (current_dataframe["REF"] == x["REF"])]
            current_alleles = [current_row["REF"].iloc[0]] + current_row["ALT"].iloc[0].split(",")
            dbsnp_info_fields = list(map(lambda y: y.split("="), current_row["INFO"].iloc[0].split(";")))
            dbsnp_info_names = [x[0] for x in dbsnp_info_fields]
            dbsnp_info_values = [x[1] if len(x) > 1 else "1" for x in dbsnp_info_fields]
            dbsnp_info_dataframe = pd.DataFrame([dbsnp_info_values], columns=dbsnp_info_names)
            result_cols = []
            for column in dbsnp_multialt_fields:
                try:
                    result_cols.append("|".join(list(map(lambda x: "=".join(x), list(
                        zip(current_alleles, dbsnp_info_dataframe[column].iloc[0].split(",")))))))
                except:
                    result_cols.append(".")
            for column in dbsnp_unique_fields:
                try:
                    result_cols.append(dbsnp_info_dataframe[column].iloc[0])
                except:
                    result_cols.append(".")
            for column in dbsnp_boolean_fields:
                try:
                    result_cols.append(dbsnp_info_dataframe[column].iloc[0])
                except:
                    result_cols.append("0")
            return result_cols
        except:
            return ["."] * (len(dbsnp_multialt_fields) + len(dbsnp_unique_fields) + len(dbsnp_boolean_fields))

    input_dataframe = pd.concat([input_dataframe, pd.DataFrame(list(current_values.apply(get_columns, axis=1).values),
                                                               columns=(
                                                                       dbsnp_multialt_fields +
                                                                       dbsnp_unique_fields +
                                                                       dbsnp_boolean_fields))],
                                axis=1)
    return input_dataframe


def add_dbnsfp_fields(input_dataframe, path_to_dbnsfp, dbnsfp_fields=[], variant_separator="&"):
    with open_potential_gz(path_to_dbnsfp) as f:
        header = f.readline().decode("utf-8").lstrip("#").rstrip("\n").rstrip("\r").split("\t")
    extended_dbnsfp_fields = dbnsfp_fields + ["Ensembl_transcriptid", "HGVSc_snpEff", "HGVSc_VEP"]

    # def pick_head(h):
    #    print(h + " in " + str(header.index(h)))
    # for dbsnp_field in extended_dbnsfp_fields:
    #    pick_head(dbsnp_field)

    import tabix
    dbnsfp = tabix.open(path_to_dbnsfp)

    def get_columns(x):
        # x is a single row of the input dataframe
        try:
            current_dbnsfp_dataframe = pd.DataFrame(
                dbnsfp.query(str(x["CHROM"]).replace("chr", ""), int(x["POS"]) - 1, int(x["POS"])),
                columns=header)
            current_row = current_dbnsfp_dataframe[
                (current_dbnsfp_dataframe["pos(1-based)"] == x["POS"]) & (
                        current_dbnsfp_dataframe["ref"] == x["REF"]) & (
                        current_dbnsfp_dataframe["alt"] == x["ALT"])]
            if current_row.shape[0] == 0:
                raise Exception("Current row is not in the database")
        except:
            returned_values = []
            for field in dbnsfp_fields:
                returned_values.append(variant_separator.join(["."] * len(x["hgvsc_id"].split(variant_separator))))
            return pd.Series(returned_values)
        current_dataframe = pd.DataFrame(
            _transpose(list(map(lambda x: x.split(";"), list(current_row[extended_dbnsfp_fields].values)[0]))),
            columns=extended_dbnsfp_fields)
        current_dataframe["vep_ids"] = list(map(lambda y: ":".join(y), list(
            zip(current_dataframe["Ensembl_transcriptid"].tolist(), current_dataframe["HGVSc_VEP"].tolist()))))
        current_dataframe["snpeff_ids"] = list(map(lambda y: ":".join(y), list(
            zip(current_dataframe["Ensembl_transcriptid"].tolist(), current_dataframe["HGVSc_snpEff"].tolist()))))

        returned_values = []
        for field in dbnsfp_fields:
            returned_values_current_field = []
            for value_to_check in x["hgvsc_id"].split("&"):
                try:
                    returned_values_current_field.append(current_dataframe[
                                                             (current_dataframe["snpeff_ids"] == value_to_check) | (
                                                                     current_dataframe[
                                                                         "vep_ids"] == value_to_check)][field].iloc[
                                                             0])
                except:
                    returned_values_current_field.append(".")
            returned_values.append(variant_separator.join(returned_values_current_field))
        return pd.Series(returned_values)

    input_dataframe[dbnsfp_fields] = pd.DataFrame(
        input_dataframe[["CHROM", "POS", "REF", "ALT", "Allele", "hgvsc_id"]].apply(get_columns, axis=1))

    # pd.DataFrame(list(current_values.apply(get_columns, axis=1).values), columns=dbnsfp_fields)

    # input_dataframe = pd.concat([input_dataframe, pd.DataFrame(list(current_values.apply(get_columns, axis=1).values),
    #                                                           columns=dbnsfp_fields)],
    #                            axis=1)
    return input_dataframe


def add_clinvar_fields(input_dataframe, path_to_clinvar, clinvar_fields=[], variant_separator="&"):
    if len(clinvar_fields) == 0:
        return input_dataframe

    import tabix
    clinvar = tabix.open(path_to_clinvar)

    def get_columns(x):
        # x is a single row of the input dataframe
        try:
            header = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
            current_clinvar_dataframe = pd.DataFrame(
                clinvar.query(str(x["CHROM"]).replace("chr", ""), int(x["POS"]) - 1, int(x["POS"])), columns=header)
            current_row = current_clinvar_dataframe[
                (current_clinvar_dataframe["POS"] == x["POS"]) & (current_clinvar_dataframe["ALT"] == x["ALT"])]
            if current_row.shape[0] == 0:
                raise Exception("Current row is not in the database")
        except:
            return pd.Series(["."] * len(clinvar_fields))

        returned_values = []
        info_field = current_row["INFO"].iloc[0]
        clinvar_dictionary = dict(list(map(lambda y: (y.split("=")[0], y.split("=")[1]), info_field.split(";"))))
        for field in clinvar_fields:
            if field == "CLINVAR_ID":
                returned_values.append(current_row["ID"].iloc[0])
            else:
                try:
                    returned_values.append(clinvar_dictionary[field])
                except:
                    returned_values.append(".")
        return pd.Series(returned_values)

    input_dataframe[clinvar_fields] = pd.DataFrame(
        input_dataframe[["CHROM", "POS", "ALT"]].apply(get_columns, axis=1))

    return input_dataframe


def add_cosmic_vcf_fields(input_dataframe, path_to_cosmic, cosmic_fields=[], variant_separator="&",
                          internal_separator="|"):
    if len(cosmic_fields) == 0:
        return input_dataframe

    vcf_header = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]

    cosmic_input_fields = [x[0] for x in cosmic_fields]
    cosmic_output_fields = [x[1] for x in cosmic_fields]

    import tabix
    cosmic = tabix.open(path_to_cosmic)

    def get_columns(x):
        # x is a single row of the input dataframe
        try:
            current_cosmic_dataframe = pd.DataFrame(
                cosmic.query(str(x["CHROM"]).replace("chr", ""), int(x["POS"]) - 1, int(x["POS"])), columns=vcf_header)
            current_row = current_cosmic_dataframe[
                (current_cosmic_dataframe["POS"] == x["POS"]) & (current_cosmic_dataframe["ALT"] == x["ALT"])]
            if current_row.shape[0] == 0:
                raise Exception("Current row is not in the database")
        except:
            return pd.Series(["."] * len(cosmic_input_fields))

        returned_values = []
        info_field_list = list(current_row["INFO"])
        cosmic_dictionaries = list(map(lambda w: dict(list(
            map(lambda y: (y.split("=")[0], y.split("=")[1]) if len(y.split("=")) == 2 else (y.split("=")[0], "true"),
                w.split(";")))), info_field_list))
        for field in cosmic_input_fields:
            def get_field(current_dictionary):
                if field == "COSMIC_ID":
                    return current_row["ID"].iloc[0]
                else:
                    try:
                        return current_dictionary[field]
                    except:
                        return "."

            returned_values.append(internal_separator.join(list(map(get_field, cosmic_dictionaries))))
        return pd.Series(returned_values)

    input_dataframe[cosmic_output_fields] = pd.DataFrame(
        input_dataframe[["CHROM", "POS", "ALT"]].apply(get_columns, axis=1))
    return input_dataframe


def add_cosmic_indexed_fields(input_dataframe, path_to_cosmic, cosmic_fields=[], cosmic_header=None,
                              variant_separator="&", internal_separator="|"):
    if len(cosmic_fields) == 0:
        return input_dataframe
    if cosmic_header is None:
        cosmic_header = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]

    cosmic_input_fields = [x[0] for x in cosmic_fields]
    cosmic_output_fields = [x[1] for x in cosmic_fields]

    import tabix
    cosmic = tabix.open(path_to_cosmic)

    def get_columns(x):
        # x is a single row of the input dataframe
        try:
            current_cosmic_dataframe = pd.DataFrame(
                cosmic.query(str(x["CHROM"]).replace("chr", ""), int(x["POS"]) - 1, int(x["POS"])),
                columns=cosmic_header)
            current_row = current_cosmic_dataframe[
                (current_cosmic_dataframe["POS"] == x["POS"]) & (current_cosmic_dataframe["ALT"] == x["ALT"])]
            if current_row.shape[0] == 0:
                raise Exception("Current row is not in the database")
        except:
            return pd.Series(["."] * len(cosmic_input_fields))

        returned_values = []
        for field in cosmic_input_fields:
            returned_values.append(internal_separator.join(list(current_cosmic_dataframe[field])))
        return pd.Series(returned_values)

    input_dataframe[cosmic_output_fields] = pd.DataFrame(
        input_dataframe[["CHROM", "POS", "ALT"]].apply(get_columns, axis=1))
    return input_dataframe


def _transpose(l):
    return list(map(list, zip(*l)))


def treat_chunk(chunk, vcf_header, ann_header,
                input_additional_multifields=None,
                vep_fields=None,
                snpeff_fields=None,
                clinvar_fields=None,
                dbNSFP_fields=None,
                dbSNP_multiallelic_fields=None,
                dbSNP_unique_fields=None,
                dbSNP_boolean_fields=None,
                cosmic_coding_mut_fields=None,
                cosmic_non_coding_vars_fields=None,
                cosmic_mutant_census_fields=None,
                renamings=None):
    if input_additional_multifields is None:
        input_additional_multifields = ["platformnames", "datasetnames", "callsetnames", "callable", "filt",
                                        "datasetsmissingcall"]
    if vep_fields is None:
        vep_fields = ["Allele", "Consequence", "IMPACT", "SYMBOL", "Gene", "Feature_type", "Feature", "BIOTYPE", "EXON",
                      "INTRON", "HGVSc", "HGVSp", "cDNA_position", "CDS_position", "Protein_position", "Amino_acids",
                      "Codons", "Existing_variation", "DISTANCE", "STRAND", "FLAGS", "SYMBOL_SOURCE", "HGNC_ID",
                      "SOURCE", "HGVS_OFFSET", "HGVSg"]
    if snpeff_fields is None:
        snpeff_fields = ["Allele", "Annotation", "Annotation_Impact", "Gene_Name", "Gene_ID", "Feature_Type",
                         "Feature_ID", "Transcript_BioType", "Rank", "HGVS.c", "HGVS.p", "cDNA.pos/cDNA.length",
                         "CDS.pos/CDS.length", "AA.pos/AA.length", "Distance", "ERRORS/WARNINGS/INFO"]
    if clinvar_fields is None:
        clinvar_fields = ["CLINVAR_ID", "CLNDISDB", "CLNDISDBINCL", "CLNDN", "CLNSIG", "CLNREVSTAT"]
    if dbNSFP_fields is None:
        """
        dbNSFP_fields = ["MutationAssessor_score", "MutationAssessor_pred", "SIFT_score", "SIFT_pred",
                         "FATHMM_score", "FATHMM_pred", "Polyphen2_HDIV_score", "Polyphen2_HDIV_rankscore",
                         "Polyphen2_HDIV_pred", "Polyphen2_HVAR_score", "Polyphen2_HVAR_rankscore",
                         "Polyphen2_HVAR_pred"]
        """
        dbNSFP_fields = ["MutationAssessor_score", "MutationAssessor_pred", "SIFT_score", "SIFT_pred", "FATHMM_score",
                         "FATHMM_pred", "Polyphen2_HDIV_score", "Polyphen2_HDIV_pred", "Polyphen2_HVAR_score",
                         "Polyphen2_HVAR_pred"]
    if dbSNP_multiallelic_fields is None:
        dbSNP_multiallelic_fields = ["TOPMED"]
    if dbSNP_unique_fields is None:
        dbSNP_unique_fields = ["SAO", "RS"]
    if dbSNP_boolean_fields is None:
        dbSNP_boolean_fields = ["ASP"]
    if cosmic_coding_mut_fields is None:
        cosmic_coding_mut_fields = ["COSMIC_ID", "LEGACY_ID", "GENE", "CDS"]
    if cosmic_non_coding_vars_fields is None:
        cosmic_non_coding_vars_fields = ["COSMIC_ID", "LEGACY_ID", "GENE", "SNP"]
    if cosmic_mutant_census_fields is None:
        cosmic_mutant_census_fields = ["Primary histology"]
    if renamings is None:
        renamings = {"cosmic_coding_mut_fields": {"GENE": "COSMIC_CODING_MUT_GENE", "COSMIC_ID": "COSMIC_ID_MUT",
                                                  "LEGACY_ID": "LEGACY_ID_MUT"},
                     "cosmic_non_coding_vars_fields": {"GENE": "COSMIC_NON_CODING_VARS_GENE",
                                                       "COSMIC_ID": "COSMIC_ID_NON_MUT",
                                                       "LEGACY_ID": "LEGACY_ID_NON_MUT"}}

    format_index = vcf_header.index("FORMAT")
    vcf_fields_header = vcf_header[:format_index]

    info_index = vcf_header.index("INFO")
    vcf_fields_header.pop(info_index)

    input_vcf_dataframe = pd.DataFrame(list(map(lambda x: x.decode("utf-8")[:-1].replace(" ", "").split("\t"), chunk)),
                                       columns=vcf_header)
    vcf_dataframe_with_info = split_info(input_vcf_dataframe)

    if not ann_header is None:
        vcf_dataframe_with_info_ann = split_ann(vcf_dataframe_with_info, ann_header)
    else:
        return vcf_dataframe_with_info

    vcf_dataframe_with_info_corrected = clean_dataframe_column(vcf_dataframe_with_info_ann,
                                                               cols=input_additional_multifields,
                                                               old_allele_sep="|", new_allele_sep="|",
                                                               old_internal_sep=",", new_internal_sep="|")
    vcf_dataframe_with_info_corrected = clean_dataframe_column(vcf_dataframe_with_info_corrected,
                                                               cols=vep_fields + snpeff_fields,
                                                               old_allele_sep=",", new_allele_sep="&",
                                                               old_internal_sep="&", new_internal_sep="|")
    # vcf_dataframe_with_snpsift_corrected = collapse_dataframe_column(vcf_dataframe_with_info_corrected,
    #                                                             single_cols=clinvar_fields,
    #                                                             multiple_cols=[], external_sep="&",
    #                                                             old_internal_sep=",", new_internal_sep="|")
    # clinvar_vep_result = [clinvar_prefix] + [clinvar_prefix + "_" + x for x in clinvar_fields]
    # vcf_dataframe_with_vep_corrected = collapse_dataframe_column(vcf_dataframe_with_snpsift_corrected,
    #                                                             single_cols=clinvar_vep_result,
    #                                                             multiple_cols=dbNSFP_fields, external_sep=",",
    #                                                             old_internal_sep="&", new_internal_sep="&")
    # vcf_dataframe_with_vep_corrected = clean_dataframe_column(vcf_dataframe_with_vep_corrected,
    #                                                          cols=clinvar_vep_result,
    #                                                          old_allele_sep=",", new_allele_sep="&",
    #                                                          old_internal_sep="&", new_internal_sep="|")
    # dbNSFP_snpsift_result = ["dbNSFP_" + x for x in dbNSFP_fields]
    # vcf_dataframe_with_splitted_cols = split_dataframe_column(vcf_dataframe_with_vep_corrected,
    #                                                          cols=dbNSFP_snpsift_result, old_internal_sep=",",
    #                                                          new_internal_sep="&", reference_column="Allele",
    #                                                          reference_separator="&")
    vcf_dataframe_with_hgvs_id = compute_hgvs_id(vcf_dataframe_with_info_corrected)

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

    path_to_dbsnp = module_path + "/data/cache/dbSNP/All_20180418.vcf.gz"
    vcf_dataframe_with_dbsnp = add_dbsnp_fields(vcf_dataframe_with_civic_gene_description, path_to_dbsnp,
                                                dbsnp_multialt_fields=dbSNP_multiallelic_fields,
                                                dbsnp_unique_fields=dbSNP_unique_fields,
                                                dbsnp_boolean_fields=dbSNP_boolean_fields)

    path_to_dbNSFP = module_path + "/data/cache/dbNSFP/dbNSFP4.0a.txt.gz"
    vcf_dataframe_with_dbnsfp = add_dbnsfp_fields(vcf_dataframe_with_dbsnp, path_to_dbNSFP, dbnsfp_fields=dbNSFP_fields)

    path_to_clinvar = module_path + "/data/cache/clinvar/clinvar.vcf.gz"
    vcf_dataframe_with_clinvar = add_clinvar_fields(vcf_dataframe_with_dbnsfp, path_to_clinvar,
                                                    clinvar_fields=clinvar_fields)

    if "cosmic_coding_mut_fields" in renamings.keys():
        cosmic_coding_mut_fields = [(x, x) if (not x in renamings["cosmic_coding_mut_fields"].keys()) else (
        x, renamings["cosmic_coding_mut_fields"][x]) for x in cosmic_coding_mut_fields]
    else:
        cosmic_coding_mut_fields = [(x, x) for x in cosmic_coding_mut_fields]
    path_to_cosmic_coding_mut = module_path + "/data/cache/cosmic/CosmicCodingMuts.vcf.gz"
    vcf_dataframe_with_cosmic_coding_mut = add_cosmic_vcf_fields(vcf_dataframe_with_clinvar,
                                                                 path_to_cosmic_coding_mut,
                                                                 cosmic_fields=cosmic_coding_mut_fields)
    if "cosmic_non_coding_vars_fields" in renamings.keys():
        cosmic_non_coding_vars_fields = [(x, x) if (not x in renamings["cosmic_non_coding_vars_fields"].keys()) else (
        x, renamings["cosmic_non_coding_vars_fields"][x]) for x in cosmic_non_coding_vars_fields]
    else:
        cosmic_non_coding_vars_fields = [(x, x) for x in cosmic_non_coding_vars_fields]
    path_to_cosmic_non_coding_vars = module_path + "/data/cache/cosmic/CosmicNonCodingVariants.vcf.gz"
    vcf_dataframe_with_cosmic_non_coding_mut = add_cosmic_vcf_fields(vcf_dataframe_with_cosmic_coding_mut,
                                                                     path_to_cosmic_non_coding_vars,
                                                                     cosmic_fields=cosmic_non_coding_vars_fields)

    if "cosmic_mutant_census_fields" in renamings.keys():
        cosmic_mutant_census_fields = [(x, x) if (not x in renamings["cosmic_mutant_census_fields"].keys()) else (
        x, renamings["cosmic_mutant_census_fields"][x]) for x in cosmic_mutant_census_fields]
    else:
        cosmic_mutant_census_fields = [(x, x) for x in cosmic_mutant_census_fields]
    path_to_cosmic_mutant_census = module_path + "/data/cache/cosmic/CosmicMutantExportCensus.tsv.gz"
    with open_potential_gz(path_to_cosmic_mutant_census) as mutant_cens:
        mutant_header = mutant_cens.readline().decode("utf-8").lstrip("#").rstrip("\n").split("\t")
    vcf_dataframe_with_cosmic_mutant_census = add_cosmic_indexed_fields(vcf_dataframe_with_cosmic_non_coding_mut,
                                                                        path_to_cosmic_mutant_census,
                                                                        cosmic_fields=cosmic_mutant_census_fields,
                                                                        cosmic_header=mutant_header)
    return vcf_dataframe_with_cosmic_mutant_census


def normalize_stream_to_stream(input_stream):
    p = Popen(["bash", "-c", "bcftools norm -m -both"], stdin=input_stream, stdout=PIPE)
    return p.stdout


def normalize_vcf_to_stream(path_to_vcf):
    p = Popen(["bash", "-c", "bcftools norm -m -both " + path_to_vcf], stdout=PIPE)
    return p.stdout


def from_vcf_stream_to_maf_stream(input_stream, output_stream, chunk_size=20):
    with normalize_stream_to_stream(input_stream) as input_stream:
        # with open_potential_gz(input_vcf) as file:

        ann_header = None
        for line in input_stream:
            line_string = line.decode("utf-8")
            print(line)
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

        for chunk in iter(lambda: list(islice(input_stream, chunk_size)), []):
            results.append(mp_pool.apply_async(treat_chunk, (chunk, vcf_header, ann_header)))

        results = list(map(lambda x: x.get(), results))
        big_dataframe = pd.concat(results, sort=False)
        all_headers = list(big_dataframe.columns)
        common_fields = [x for x in all_headers if not x in individual_headers]
        new_order = common_fields + individual_headers
        big_dataframe = big_dataframe[new_order]
        big_dataframe.to_csv(output_stream, sep='\t', na_rep='.', index=False)


def from_vcf_to_maf(input_vcf, output_maf, chunk_size=20):
    """Converts vcf to maf.
    Parameters
    ----------
    input_vcf: String
        Path to the input vcf (can be gzipped) but the format is not checked
    output_maf: String
        Path where the MAF should be stored
    """
    with open_potential_gz(input_vcf) as input_stream, open(output_maf, "w") as output_stream:
        from_vcf_stream_to_maf_stream(input_stream, output_stream)


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


def annotate_stream_to_stream(input_stream, output_stream, json_file):
    """Annotates the input_stream to the output_stream following the configuration
    present in json_file.
    Parameters
    ----------
    input_stream : TextIOBase
        Input stream.
    output_stream : TextIOBase
        Output stream.
    config : dict
        Dictionary containing the necessary information to perform the annotation
    """

    module_path = dirname(dirname(abspath(__file__)))

    with open(json_file, 'r') as f:
        config_json = json.load(f)

    with open(module_path + "/config/schema.json") as f:
        schema_json = json.load(f)

    _check_json_file(config_json, schema_json)

    if len(config_json["prediction_tool"]) > 0:
        if "vep" in config_json["prediction_tool"]:
            from pyannotation.pyvep import annotate_stream_vep
            p = annotate_stream_vep(input_stream, PIPE)
            from_vcf_stream_to_maf_stream(p.stdout, output_stream)
        elif "snpeff" in config_json["prediction_tool"]:
            from pyannotation.pysnpeff import annotate_stream_snpeff
            p = annotate_stream_snpeff(input_stream, PIPE)
            from_vcf_stream_to_maf_stream(p.stdout, output_stream)
    else:
        from_vcf_stream_to_maf_stream(input_stream, output_stream)


def annotate(input_file, output_file, json_file):
    """Annotates the input_file to the output_file following the configuration
    present in json_file.
    Parameters
    ----------
    input_stream : String
        Path to the input file.
    output_stream : String
        Path to the output file.
    config : dict
        Dictionary containing the necessary information to perform the annotation
    """

    module_path = dirname(dirname(abspath(__file__)))

    with open(json_file, 'r') as f:
        config_json = json.load(f)

    with open(module_path + "/config/schema.json") as f:
        schema_json = json.load(f)

    _check_json_file(config_json, schema_json)

    with open(output_file, "w") as output_stream:
        if is_gz_file(input_file):
            p = Popen(["bash", "zcat", input_file], stdout=PIPE)
            input_stream = p.stdout
        else:
            p = Popen(["bash", "cat", input_file], stdout=PIPE)
            input_stream = p.stdout
        annotate_stream_to_stream(input_stream, output_stream, json_file)


def main(input_file, output_file, json_file):
    return annotate(input_file, output_file, json_file)


if __name__ == "__main__":
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    json_file = sys.argv[3]
    start_time = time.time()
    main(input_file, output_file, json_file)
    print("Ellaspsed time: " + str(time.time() - start_time))
