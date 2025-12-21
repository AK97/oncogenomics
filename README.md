# Cancer Genomics

Work in progress exploring mutational signatures across four TCGA cancers (BRCA, GBM, KIRC, LUAD) and preparing data for pathway-level analysis.

## Current status
- Clean SNV matrices built: genomic trinucleotide mutations are linked to individual tumor samples and projected into COSMIC 96-channel space.
- Coverage spans the four cancers above using merged TCGA MAFs; low-mutation samples are filtered out for downstream stability.
- RNA-seq counts (STAR, gene-level) pulled for the retained samples and normalized to log1p CPM to support pathway association work.
- Next step: associate mutation exposures with pathway activity scores.

## Pipeline outline
1. **MAF ingestion**: `load_and_merge_maf_files()` reads downloaded TCGA MAFs (`data/raw/`, manifest `gdc_manifest.txt`), merges to `merged_data.parquet`.
2. **Trinucleotide context**: `import_reference_genome()` + `add_trinuc_context()` annotate SNVs with GRCh38 context (`data/data_with_trinuc.parquet`).
3. **COSMIC channels**: `create_cosmic_channels()` and `create_96_matrix()` map SNVs into the 96-channel matrix; samples with <100 SNVs are dropped via `filter_snv_counts()`.
4. **Signature exposures**: `import_cosmic_signatures()` aligns COSMIC v3.5 SBS signatures, and `run_nnls()` learns per-sample fractional exposures.
5. **RNA-seq**: `get_rna_seq_for_samples()` queries GDC for matching STAR counts; `build_counts_matrix()` and `normalize_and_filter_rna_seq()` yield a filtered expression matrix (~543 samples x 17k genes).

Run `python main.py` to execute the full pipeline. Outputs are cached to parquet where possible (see `data/`).

## Inputs
- **TCGA MAFs**: place downloaded folders under `data/raw/`.
- **Reference genome**: GRCh38 FASTA + index in `data/`.
- **COSMIC SBS v3.5**: `data/COSMIC_v3.5_SBS_GRCh38.txt` (parquet cache written automatically).

## Next steps
- Compute pathway activity scores (e.g., via GSVA/ssGSEA) from normalized RNA-seq and join to mutation exposures.
- Evaluate signature-pathway associations within and across cancer types.
- Harden reproducibility: pin data versions in `config.yaml`, add CLI entrypoint/tests.
