import pandas as pd

from src.merge_mafs import load_and_merge_maf_files
from src.trinucleotide_context import import_reference_genome, add_trinuc_context, create_cosmic_channels, create_96_matrix, sample_proj_snp_key
from src.channels import CHANNELS_96
from src.cosmic import filter_snv_counts, import_cosmic_signatures, run_nnls
from src.rna_seq import get_rna_seq_for_samples, build_counts_matrix, normalize_and_filter_rna_seq
from src.pathway_scoring import import_pathway_map, score_ssgsea
from src.modeling import clr, standardize_pathways

tumor_maf_data = load_and_merge_maf_files()

ref_genome = import_reference_genome()

with_trinucleotide_context = add_trinuc_context(tumor_maf_data, ref_genome)

cosmic_canonicals = create_cosmic_channels(with_trinucleotide_context)

sample_key = sample_proj_snp_key(with_trinucleotide_context)
mat96 = create_96_matrix(cosmic_canonicals)

# store metadata as maps. use 16 dig sample ids
n_snv_map = {k[:16]: v for k, v in mat96.attrs["n_snv"].items()}
cancer_type_map = {k[:16]: v for k, v in mat96.attrs["cancer_type"].items()}

filtered_snv_counts = filter_snv_counts(mat96)
# at threshold 100, 612 samples remain

cosmic_signatures = import_cosmic_signatures()

learned_exposures = run_nnls(filtered_snv_counts, cosmic_signatures)

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

'''
Now we have
- mutational signatures (learned_exposures)
- pathway activity (ssgsea_scores)
- sample metadata (project_name, n_snv)
'''

# splice sample ids to first 16 chars, to match
learned_exposures.index = learned_exposures.index.astype(str).str[:16]
# average scores to remove what are now replicates
learned_exposures = learned_exposures.groupby(level=0, sort=False).mean()

# overlapping samples - those with mutational signatures and pathway data
common_samples = learned_exposures.index.intersection(ssgsea_scores.index)

# 612 with mut sigs, 543 with pathway data; 543 overlap

# use centered log transform to avoid multilinearity in SBS data
learned_exposures = clr(learned_exposures)

# print(learned_exposures.head())
sp = standardize_pathways(ssgsea_scores)

# combine data thus far into one matrix (543 x 149)
aggregated_data = pd.concat(
    [
        learned_exposures,
        sp,
        pd.Series(n_snv_map, name="n_snv"),
        pd.Series(cancer_type_map, name="cancer_type"),
    ],
    axis=1,
    join="inner",
)
print(aggregated_data.shape)