from pathlib import Path
import numpy as np
import pandas as pd
import pyfaidx
from src.channels import CHANNELS_96

def import_reference_genome() -> pyfaidx.Fasta:
    """
    Import reference genome from FASTA file using pyfaidx.
    Reference genome taken from GDC: GRCh38.
    Only autosomes are used downstream.
    """
    fasta_file = "data/GRCh38.d1.vd1.fa"
    # as_raw returns plain strings which keeps downstream indexing fast.
    return pyfaidx.Fasta(fasta_file, as_raw = True)


def _trinuc_for_chrom(chrom_df: pd.DataFrame, ref_seq) -> pd.DataFrame:
    """
    Vectorised trinucleotide lookup for a single chromosome.
    """
    # as_raw=True means FastaRecord doesn't expose .seq; slice to get plain string
    seq = ref_seq[:]
    positions = chrom_df["Start_Position"].to_numpy(dtype = np.int64)

    # Convert 1-based positions to 0-based indices for Python string access
    prev_bases = [seq[p - 2] for p in positions]
    ref_bases = [seq[p - 1] for p in positions]
    next_bases = [seq[p] for p in positions]

    result = pd.DataFrame(
        {
            "prev_base": prev_bases,
            "ref_base": ref_bases,
            "next_base": next_bases,
        },
        index = chrom_df.index,
    )
    result["ref_match"] = result["ref_base"].to_numpy() == chrom_df[
        "Reference_Allele"
    ].to_numpy()
    return result


def add_trinuc_context(snv_data: pd.DataFrame, ref_genome: pyfaidx.Fasta) -> pd.DataFrame:
    """
    Add trinucleotide context to SNV data.
    For each Start_Position, find nucleotide at, before, and after the position
    in the reference genome. Trinucleotide context is a string of 3 nucleotides,
    with the SNV nucleotide in the middle.
    Note: SNV positions are 1-based, but string indexing is 0-based.
    """

    parquet_filepath = "data/data_with_trinuc.parquet"
    if Path(parquet_filepath).exists():
        return pd.read_parquet(parquet_filepath)

    contexts = []
    for chrom, chrom_df in snv_data.groupby("Chromosome", sort = False):
        # if chrom not in ref_genome:
        #     raise ValueError(f"Chromosome {chrom} not found in reference genome")
        contexts.append(_trinuc_for_chrom(chrom_df, ref_genome[chrom]))

    context_df = pd.concat(contexts).sort_index()
    with_trinuc = snv_data.join(context_df)

    # Drop all rows where the reference allele does not match the reference genome
    with_trinuc = with_trinuc[with_trinuc["ref_match"]]

    # Save to parquet for future use
    with_trinuc.to_parquet(parquet_filepath)

    return with_trinuc

def create_cosmic_channels(snv_data: pd.DataFrame) -> pd.DataFrame:
    """
    Get mutation for each row in COSMIC-compatible format (pyrimidine context).
    """

    COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}

    # COSMIC format uses pyrimidine context (C or T as reference base)
    # Reverse complement if reference base is a purine (A or G)
    reverse_complement_needed = snv_data['ref_base'].isin(['A', 'G'])

    left  = snv_data["prev_base"].to_numpy()
    ref   = snv_data["ref_base"].to_numpy()
    right = snv_data["next_base"].to_numpy()
    alt   = snv_data["Tumor_Seq_Allele2"].to_numpy()

    complement_vec = np.vectorize(COMPLEMENT.get)

    # Avoid modifying the original arrays
    left_c  = left.copy()
    ref_c   = ref.copy()
    right_c = right.copy()
    alt_c   = alt.copy()

    # Apply reverse complement transformation for purine contexts
    idx = np.where(reverse_complement_needed.to_numpy())[0]
    if len(idx) > 0:
        left_c[idx]  = complement_vec(right[idx])
        ref_c[idx]   = complement_vec(ref[idx])
        right_c[idx] = complement_vec(left[idx])
        alt_c[idx]   = complement_vec(alt[idx])

    snv_data["left_c"] = left_c
    snv_data["ref_c"] = ref_c # should only be C or T
    snv_data["right_c"] = right_c
    snv_data["alt_c"] = alt_c
    snv_data["sub_c"] = snv_data["ref_c"] + ">" + snv_data["alt_c"]
    snv_data["channel96"] = snv_data["left_c"] + "[" + snv_data["ref_c"] + ">" + snv_data["alt_c"] + "]" + snv_data["right_c"]

    return snv_data

def sample_proj_snp_key(snv_data: pd.DataFrame) -> pd.DataFrame:
    '''
    Return df with columns: Tumor_Sample_Barcode, project_name, N_snv
    '''
    return (
        snv_data
        .groupby(["Tumor_Sample_Barcode", "project_name"])
        .size()
        .rename("N_snv")
        .reset_index()
    )

def create_96_matrix(snv_data: pd.DataFrame) -> np.ndarray:
    """
    Create 96-channel matrix from SNV data.
    Each row corresponds to a single SNV, and each column corresponds to a channel.
    Channels are ordered as per source of truth file (CHANNELS_96).
    """

    # Long-form counts
    counts_long = (
        snv_data.groupby(["Tumor_Sample_Barcode", "channel96"])
        .size()
        .rename("count")
        .reset_index()
    )

    # samples x channels
    counts_wide = (
        counts_long.pivot(
            index = "Tumor_Sample_Barcode", columns = "channel96", values="count"
        )
        .fillna(0)
        .astype(np.int32) # no need to be floats
    )

    # Freeze ordering
    counts_wide = counts_wide.reindex(columns = CHANNELS_96, fill_value = 0)

    # Count SNVs per sample, and store as metadata
    counts_wide.attrs["n_snv"] = counts_wide.sum(axis=1).to_dict()

    # Also store cancer type as metadata
    cancer_type_dict = (
        snv_data.groupby("Tumor_Sample_Barcode", sort=False)["project_name"]
        .first()
        .to_dict()
    )
    counts_wide.attrs["cancer_type"] = cancer_type_dict

    return counts_wide