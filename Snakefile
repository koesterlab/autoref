import pandas as pd
from snakemake.remote import HTTP, FTP

http = HTTP.RemoteProvider()
ftp = FTP.RemoteProvider()

configfile: "config.yaml"

references = pd.read_table(config["references"])

rule all:
    input:
        expand("{ref.species}/{ref.source}/{ref.build}/{ref.type}", ref=references.itertuples()),
        "ncbi-pipelines-genomes/homo_sapiens"        


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


rule ncbi_pipelines_genomes:
    """This is the best for human according to Heng Li:
       https://lh3.github.io/2017/11/13/which-human-reference-genome-to-use. 
       In addition we convert to ENSEMBL naming
       for compatibility with the rest of the ensembl annotation and in particular SnpEff."""
    output:
        directory("ncbi-pipelines-genomes/homo_sapiens")
    conda:
        "envs/bwa.yaml"
    shell:
        """
        mkdir {output};
        cd {output}
        wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz
        gzip -d GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz
        sed -i 's/chr//' GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna
        samtools faidx GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna
        bwa index GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna
        """

