# PyAnnotation

The goal of this package is to annotate `.vcf` and `.vcf.gz` files either from the command line or from a Python script. 

**Under development**

## Table of contents

1. [Introduction](#intro)
2. [Installation](#install)
3. [Command line sintax](#comline)
4. [JSON sintax](#jsonsintax)
5. [Authors](#authors)

## 1. Introduction <a name="intro"></a>

First of all, the effect of each variant is predicted with either `snpEff` or `VEP`. Afterwards, this information is completed with fields coming from one or several of the following databases:

* [ClinVar](ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/)
* [dbNSFP](https://sites.google.com/site/jpopgen/dbNSFP)
* [dbSNP](ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/)
* [CiVIC gene summaries](https://civicdb.org/downloads/nightly/nightly-GeneSummaries.tsv)
* [CiVIC variant summaries](https://civicdb.org/downloads/nightly/nightly-VariantSummaries.tsv)
* [CiVIC clinical evidence](https://civicdb.org/downloads/nightly/nightly-ClinicalEvidenceSummaries.tsv)
* [COSMIC coding mutations](https://cancer.sanger.ac.uk/cosmic/download/CosmicCodingMuts.vcf.gz)
* [COSMIC non coding mutations](https://cancer.sanger.ac.uk/cosmic/download/CosmicNonCodingVariants.vcf.gz)
* [COSMIC mutant census](https://cancer.sanger.ac.uk/cosmic/download/CosmicMutantExportCensus.tsv.gz)

In order to define the desired fields from each database and the effect predictor to use, it is possible to supply them in two different ways. Through the command line parameters or with a JSON file.

## 2. Installation <a name="install"></a>

In order to get the package, it is enough to follow this steps:

1. Clone this repository in a local folder.
```
git clone https://github.com/ramonamela/PyAnnotation.git
```

2. Install the dependencies.
```
./PyAnnotation/dependencies/install_dependencies.sh
```

3. Download the databases.
```
./PyAnnotation/useful_commands/download_cache_files.sh
```
If the databases are already present in the system, it is possible to create soft link to the files containing them. It is important to realise that, for the most part of the sources, the files must be indexed with `tabix`. This is done automatically with the previous command.

4. Add the package folder to the PYTHONPATH.
```
echo "export PYTHONPATH=$(pwd)/PyAnnotation:\$PYTHONPATH" >> ~/.bashrc; source ~/.bashrc
```

## 3. Command line sintax <a name="comline"></a>

It is possible to annotate the files from the command line with the following command:

```
python3 /path/to/module/pyannotation.py -i /path/to/input.vcf -o /path/to/output.vcf -j /path/to/json.json
```

With the package there are two predefined JSON files:

* /path/to/module/config/default.json

Recommended JSON file with the considered most important fields from the different databases.

* /path/to/module/config/simple.json

Minimal example with some fields of each one of the different source databases.

On the other hand, it is possible to specify the required fields in the command line instead than giving a JSON file:

```
python3 /path/to/module/pyannotation.py -i /path/to/input.vcf -o /path/to/output.vcf -j /path/to/json.json
```

If the flag `-i` is not specified, the input is expected through the standard input. In a similar way, when the flag `-o` is not specified, the output is given in the standard output. This fact makes possible using this package in a piped workflow.

It is possible to see a full description of all the available options with the following command:

```
python3 /path/to/module/pyannotation.py -h
```

## 4. JSON sintax <a name="jsonsintax"></a>

In this section, a description of the fields found in a correct JSON file are described. The ```simple.json``` file is taken as reference. 

There are three main fields in the JSON file

### 1. Prediction tool
List with all the effect prefictors to use. Currently, `snpeff` and `vep` are accepted and only one value can be given. This field is implemented as a list since it is planned to extend the package to merge the files coming from several prediction tools. 
```
"prediction_tool": [
    "vep"
]
```

### 2. Annotated fields
List with all the desired annotations from the different databases. The format to specify a certain annotation is the following one:
```
{
    "database": "database_name",
    "fields": [ list_of_fields ]
}
```
The possible values for `database` are dbClinVar, dbCiVICGeneSummaries, dbCiVICVariantSummaries, dbCiVICClinicalEvidence, dbNSFP, dbSNP, dbCOSMICCodingMutations, dbCOSMICNonCodingVariants and dbCOSMICMutantCensus. Afterwards, the fields are checked to verify that the desired parameters are indeed in the input files.

In addition, there are some databases that need previous annotations in order to relate each variant with the database. This is the case for the fields present in dbCiVICClinicalEvidence. In order to get the fields from this database, it is mandatory to annotate `gene` and `variant` from dbCiVICVariantSummaries. 

### 3. Renamings
There is also the possibility to rename some field in order to avoid collisions between different databases. This field is called `all_renamings` in the JSON file and it is a list of the renamings to be done in each database. The format for each one of the sources is the following one:
```
{
    "database": "database_name",
    "database_renamings": [
        list_of_renamings
    ]
}
```
In addition, each renaming is defined as follows:
```
{
    "original_name": "original_field_name"
    "new_name":      "renamed_field_name"
}
```
Hence, is an example of a valid renaming field:
```
  "all_renamings": [
    {
      "database": "dbClinVar",
      "database_renamings": [
        {
          "original_name": "AF_ESP",
          "new_name": "AF_ESP_ClinVar"
        },
        {
          "original_name": "CLNDISDB"
          "new_name": "CLNDISDB_ClinVar"
        }
      ]
    },
    {
      "database": "dbCiVICGeneSummaries",
      "database_renamings": [
        {
          "original_name": "description",
          "new_name": "civic_gene_description"
        }
      ]
    }
  ]
```

## 5. Authors <a name="authors"></a>

* **Ramon Amela Milian** - [ORCID](https://orcid.org/0000-0001-5943-5824) [github](https://github.com/ramonamela)