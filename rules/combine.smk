rule combine_rsem:
    input:
        expand("results/RSEM/{unit.sample_name}-{unit.unit_name}.isoforms.results",unit=units.itertuples())
    output:
        expected_count = "results/combined/RSEM_expected_counts.tsv",
        TPM = "results/combined/RSEM_TPM.tsv",
        FPKM = "results/combined/RSEM_FPKM.tsv",
    threads: 4
    run:
        data_dict = {}
        cols = ["transcript_id", "gene_id", "length"]
        for data in ['expected_count', 'TPM', 'FPKM']:
            show_output(f"Combining {data} to {output[data]}")
            # load the df of shared rows into the data dict
            data_dict[data] = pd.read_csv(input[0], sep="\t").loc[:, cols]
            # add the 
            data_cols = cols + [data]
            for f in input:
                # retrieve the sample name
                sample_name = f.split("/")[-1].split("-")[0]
                show_output(f"Reading {f}")
                sample_df = pd.read_csv(f, sep="\t").loc[:, data_cols].rename({data:sample_name}, axis=1)
                # print(sample_df.columns)
                data_dict[data] = data_dict[data].merge(sample_df, how="outer")
            data_dict[data].to_csv(output[data], sep="\t", index=False)
            show_output(f"Combined {data} data written to {output[data]}", color="success")


rule combine_STARcounts:
    input:
        expand(
            "results/counts/{unit.sample_name}-{unit.unit_name}_ReadsPerGene.out.tab",
            unit=units.itertuples(),
        ),
    output:
        "results/combined/STAR_counts.tsv",
    log:
        "logs/count-matrix.log",
    params:
        samples=units["sample_name"].tolist()
    script:
        "../scripts/count-matrix.py"
