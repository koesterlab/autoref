from snakemake.shell import shell
import os

with open(snakemake.input.manifest[0]) as manifest:
    if not any(url.startswith(snakemake.params.url) for url in manifest):
        raise RuntimeError("File {} not present it igenomes. Please check the manifest: {}".format(snakemake.params.url, snakemake.input.manifest))

shell("aws s3 --no-sign-request --region eu-west-1 sync "
      "{snakemake.params.url} {snakemake.output}")

if not os.listdir(snakemake.output[0]):
    raise RuntimeError("No file downloaded.")
