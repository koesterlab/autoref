import pandas as pd
from snakemake.remote import HTTP

http = HTTP.RemoteProvider()

configfile: "config.yaml"

references = pd.read_table(config["references"])

rule all:
    input:
        expand("{ref.species}/{ref.source}/{ref.build}/{ref.type}", ref=references.itertuples())


def get_url(wildcards):
    category = "Sequence"
    if wildcards.type in ["smRNA", "Variation", "Genes"]:
        category = "Annotation"
    if wildcards.type == "bundle":
        pattern = ("s3://ngi-igenomes/igenomes/{species}/"
                   "{source}/{build}/")
    else:
        pattern = ("s3://ngi-igenomes/igenomes/{species}/"
                   "{source}/{build}/{category}/{type}/")
    return pattern.format(category=category, **wildcards)


rule download:
    input:
        manifest=http.remote("https://raw.githubusercontent.com/ewels/AWS-iGenomes/master/ngi-igenomes_file_manifest.txt", keep_local=True)
    output:
        directory("{species}/{source}/{build}/{type}")
    params:
        url=get_url
    wildcard_constraints:
        species="[^/]+",
        build="[^/]+",
        type="[^/]+",
        source="[^/]+"
    conda:
        "envs/aws.yaml"
    script:
        "scripts/download.py"
