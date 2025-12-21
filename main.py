from src.merge_mafs import load_and_merge_maf_files
from src.trinucleotide_context import import_reference_genome, add_trinuc_context, create_cosmic_channels, create_96_matrix, sample_proj_snp_key
from src.channels import CHANNELS_96
from src.cosmic import filter_snv_counts, import_cosmic_signatures, run_nnls
from src.rna_seq import get_rna_seq_for_samples, build_counts_matrix, normalize_and_filter_rna_seq


data = load_and_merge_maf_files()

ref_genome = import_reference_genome()

with_trinucleotide_context = add_trinuc_context(data, ref_genome)

cosmic_canonicals = create_cosmic_channels(with_trinucleotide_context)

sample_key = sample_proj_snp_key(with_trinucleotide_context)
mat96 = create_96_matrix(cosmic_canonicals)

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

