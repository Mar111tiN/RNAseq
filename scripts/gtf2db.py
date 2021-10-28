import gffutils

db = gffutils.create_db(
    snakemake.params.gtf,
    dbfn=snakemake.output.db,
    force=True,
    keep_order=True,
    merge_strategy="merge",
    sort_attribute_values=True,
    disable_infer_genes=True,
    disable_infer_transcripts=True,
)