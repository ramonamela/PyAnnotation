#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
import sys
from pyannotation.utils import open_potential_gz
from multiprocessing import Pool
import time

def parallelize_dataframe(df, func, n_cores=8):
    df_split = np.array_split(df, n_cores)
    pool = Pool(n_cores)
    df = pd.concat(pool.map(func, df_split))
    pool.close()
    pool.join()
    return df

def transform_dataframe(input_dataframe):
    ## Clean rows without an assignated position
    rows_to_modify_mask = (~input_dataframe["Mutation genome position"].isnull()) & (~input_dataframe["Mutation CDS"].isnull()) & list(map(lambda x: ">" in x, input_dataframe["Mutation CDS"]))
    input_dataframe = input_dataframe.reindex(columns=["#CHROM", "POS", "ALT"] + list(input_dataframe.columns))
    ## TRY TO DO THIS ONE BY ONE INSTEAD THAT WITH APPLY AND PROFILE IT IN CASE IS FASTER
    def generate_coordinates(x):
        first_split = x["Mutation genome position"].split(":")
        splitted_string = x["Mutation CDS"].split(">")
        alt = splitted_string[1]
        return pd.Series([first_split[0], first_split[1].split("-")[0], alt])
    input_dataframe[["#CHROM", "POS", "ALT"]] = input_dataframe[rows_to_modify_mask].apply(generate_coordinates, axis=1)
    return input_dataframe[rows_to_modify_mask]

def transform_tsv(input_file, output_file):
    beg_time = time.time()
    with open_potential_gz(input_file) as in_file:
        input_dataframe = pd.read_csv(in_file, sep="\t")
    read_time = time.time()
    result_dataframe = parallelize_dataframe(input_dataframe, transform_dataframe)
    create_cols_time = time.time()
    result_dataframe["POS"] = pd.to_numeric(result_dataframe["POS"])
    result_dataframe = result_dataframe.sort_values(by=["#CHROM", "POS"])
    sort_time = time.time()
    result_dataframe.to_csv(output_file, sep="\t", na_rep=".", index=False)
    write_time = time.time()

if __name__ == "__main__":
    if not len(sys.argv) == 3:
        print("Incorrect amount of input parameters")
    else:
        input_file = sys.argv[1]
        output_file = sys.argv[2]
        transform_tsv(input_file, output_file)