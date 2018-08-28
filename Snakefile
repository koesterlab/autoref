import pandas as pd

configfile: "config.yaml"

references = pd.read_table(config["references"])
print(list(references.itertuples()))

rule all:
    input:
        expand("{ref.species}/{ref.build}/{ref.type}", ref=references.itertuples())


def get_category(wildcards):
    if wildcards.type in ["smRNA", "Variation", "GTF", "BED"]:
        return "Annotation"
    else:
        return "Sequence"


rule download:
    output:
        directory("{species}/{build}/{type}")
    params:
        category=get_category,
        source=config["source"]
    wildcard_constraints:
        species="[^/]+",
        build="[^/]+",
        type="[^/]+"
    conda:
        "envs/aws.yaml"
    shell:
        "aws s3 --no-sign-request --region eu-west-1 sync s3://ngi-igenomes/igenomes/{wildcards.species}/{params.source}/{wildcards.build}/{params.category}/{wildcards.type}/ {output}"
