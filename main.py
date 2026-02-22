import pandas as pd

from src.merge_mafs import load_and_merge_maf_files
from src.trinucleotide_context import import_reference_genome, add_trinuc_context, create_cosmic_channels, create_96_matrix, sample_proj_snp_key
from src.cosmic import filter_snv_counts, import_cosmic_signatures, run_nnls, build_signature_channel_matrix
from src.rna_seq import get_rna_seq_for_samples, build_counts_matrix, normalize_and_filter_rna_seq
from src.pathway_scoring import import_pathway_map, score_ssgsea
from src.modeling import clr, standardize_pathways
from src.regression import create_linreg_model

tumor_maf_data = load_and_merge_maf_files()

ref_genome = import_reference_genome()

with_trinucleotide_context = add_trinuc_context(tumor_maf_data, ref_genome)

cosmic_canonicals = create_cosmic_channels(with_trinucleotide_context)

# sample_key = sample_proj_snp_key(with_trinucleotide_context)
mat96 = create_96_matrix(cosmic_canonicals)
# samples x 96 matrix of counts, with metadata in attrs: n_snv, cancer_type

filtered_snv_counts = filter_snv_counts(mat96)
# at threshold 100, 612 samples remain; 612 x 96

# store metadata as maps. use 16 dig sample ids
n_snv_map = {k[:16]: v for k, v in filtered_snv_counts.attrs["n_snv"].items()}
cancer_type_map = {k[:16]: v for k, v in filtered_snv_counts.attrs["cancer_type"].items()}

cosmic_signatures = import_cosmic_signatures()

learned_exposures = run_nnls(filtered_snv_counts, cosmic_signatures)
learned_exposures_original = learned_exposures.copy()
# samples x fractional exposures; 612 x 97

signature_channel_matrix = build_signature_channel_matrix(cosmic_signatures, learned_exposures)
# 97 signatures x 96 channels of counts

rna_seq_data_mapping = get_rna_seq_for_samples(learned_exposures)

# write_manifest(rna_seq_data_mapping)
# downloaded files are in data/raw_rna_seq/

# build matrix of samples x genes
counts_matrix = build_counts_matrix(rna_seq_data_mapping)
# 543 x 60660

normalized_counts = normalize_and_filter_rna_seq(counts_matrix)
# 543 x 17635

# pathway enrichment analysis
pathway_map = import_pathway_map()
ssgsea_scores = score_ssgsea(normalized_counts, pathway_map)


# Now we have
# - mutational signatures (learned_exposures)
# - pathway activity (ssgsea_scores)
# - sample metadata (cancer_type, n_snv)

# splice sample ids to first 16 chars, to match
learned_exposures.index = learned_exposures.index.astype(str).str[:16]
# average scores to remove what are now replicates
learned_exposures_deduped = learned_exposures.groupby(level=0, sort=False).mean()

# overlapping samples - those with mutational signatures and pathway data
common_samples = learned_exposures_deduped.index.intersection(ssgsea_scores.index)

# 612 with mut sigs, 543 with pathway data; 543 overlap

# use centered log transform to avoid multilinearity in SBS data # but not for frac exp visualization
learned_exposures_clr = clr(learned_exposures_deduped.loc[common_samples])

# print(learned_exposures.head())
sp = standardize_pathways(ssgsea_scores) # (543 x 50)

# Data export for mutational signature landscape visualization
from src.output_data import output_mut_sig_landscape_data
output_mut_sig_landscape_data(
    (learned_exposures_deduped, n_snv_map, cancer_type_map),
    sp
)

# Data export for SBS QC check
from src.output_data import output_sbs_qc_check_data
output_sbs_qc_check_data(
    filtered_snv_counts,
    learned_exposures_original,
    signature_channel_matrix,
    filtered_snv_counts.attrs["n_snv"],
    filtered_snv_counts.attrs["cancer_type"]
)


# # linear regression
# model = create_linreg_model(
#     (learned_exposures_clr, n_snv_map, cancer_type_map),
#     sp
# )
# print("Model Coefficients:")
# print(model.coef_)
# print("\nModel Intercept:")
# print(model.intercept_)
# print("\nModel Parameters (alpha):")
# print(model.alpha) # hyperparam
