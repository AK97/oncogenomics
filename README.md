# Cancer Genomics
Multistep oncogenomic bioinformatics project. Aggregates and assembles data (mutational, transcriptomic, pathway, genomic) to perform some analysis on cancer genotypes and transcriptional consequences.

## Current Status
- Built clean SNV matrices. Genomic trinucleotide mutations are linked to individual tumor samples and projected into COSMIC 96-channel space.
- Coverage spans the four cancers above using merged TCGA MAFs.
- RNA-seq counts (STAR) pulled for the filtered samples and normalized.
- Signature exposures are estimated per sample and formatted for expression-derived pathway signals.
- Computed pathway activity scores (ssGSEA) from normalized RNA-seq and joined to mutation exposures.
- Compiled signature-pathway associations within and across cancer types.

## Functional Pipeline / Module Breakdown
1. **MAF ingestion**: `load_and_merge_maf_files()` reads downloaded TCGA MAFs, merges to `merged_data.parquet`.
2. **Trinucleotide context**: `import_reference_genome()` + `add_trinuc_context()` annotate SNVs with GRCh38 context.
3. **COSMIC channels**: `create_cosmic_channels()` and `create_96_matrix()` map SNVs into the 96-channel matrix; Potentially noisy samples are dropped in `filter_snv_counts()`.
4. **Signature exposures**: `import_cosmic_signatures()` aligns COSMIC SBS signatures, and `run_nnls()` learns per-sample fractional exposures.
5. **RNA-seq**: `get_rna_seq_for_samples()` queries GDC for matching STAR counts; `build_counts_matrix()` and `normalize_and_filter_rna_seq()` yields a filtered expression matrix.
6. **Pathway scoring**: `import_pathway_map()` + `score_ssgsea()` compute ssGSEA pathway activity per sample and standardize for modeling.

## Inputs
MAF files, GRCh38 FASTA, COSMIC signatures, RNA-seq quant files, pathway maps.

## Running
Run `python main.py` to execute the full pipeline. Need input data files (large, not uploaded here). Outputs are cached to parquet where possible.

## Ongoing Steps
- Construct a regression model that links signature exposures to pathway activity.
- Build visualizations: signature profiles, exposure heatmaps, pathway activity and exposure scatterplots.
