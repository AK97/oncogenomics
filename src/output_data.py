import numpy as np
import pandas as pd

def output_mut_sig_landscape_data(inputs: tuple, responses: pd.DataFrame):
    sbs = inputs[0]
    n_snv = inputs[1]
    c_type = inputs[2]

    predictors = pd.concat(
        [
            sbs,
            np.log10(pd.Series(n_snv, name="n_snv")),
            pd.Series(c_type, name="cancer_type")
        ],
        axis = 1,
        join = "inner"
    )

    common_samples = predictors.index.intersection(responses.index)
    predictors = predictors.loc[common_samples]
    responses = responses.loc[common_samples]

    data_export = pd.concat(
        (predictors, responses),
        axis = 1
    )
    data_export.to_csv('output/aggr_data_mut_sig_landscape.csv')

def output_sbs_qc_check_data(signatures, exposures, signature_channel_matrix, snv_counts, cancer_type_map):
    '''
    Reconstruct the original mutation counts from the learned exposures and signature channel matrix, and compute cosine similarity and reconstruction error for each sample.
    Writes CSV with sample, cosine similarity, and reconstruction error.
    Writes CSV with cosine similarity, snv count, and cancer type metadata for each sample.
    '''

    # Reconstruct original mutation counts
    reconstructed = exposures @ signature_channel_matrix

    # Compute cosine similarity
    dot_products = (signatures * reconstructed).sum(axis=1)
    norms_signatures = np.linalg.norm(signatures.values, axis=1)
    norms_reconstructed = np.linalg.norm(reconstructed.values, axis=1)
    cosine_similarities = dot_products / (norms_signatures * norms_reconstructed)

    # Compute reconstruction error
    reconstruction_errors = np.linalg.norm((signatures - reconstructed).values, axis=1)

    # Create summary dataframe
    summary_df = pd.DataFrame({
        "sample": signatures.index,
        "cosine_similarity": cosine_similarities,
        "reconstruction_error": reconstruction_errors
    })

    summary_df.to_csv('output/sbs_qc_check_summary.csv', index=False)

    # Write additional CSV for cosine similarity, SNV count, and cancer type
    # log10 snv count for better visualization
    snv_counts = np.log10(snv_counts)
    snv_count_series = pd.Series(snv_counts, name="n_snv")
    cancer_type_series = pd.Series(cancer_type_map, name="cancer_type")

    snv_count_df = summary_df[["sample", "cosine_similarity"]].merge(
        snv_count_series.rename_axis("sample").reset_index(),
        on="sample",
        how="left"
    ).merge(
        cancer_type_series.rename_axis("sample").reset_index(),
        on="sample",
        how="left"
    )

    snv_count_df.to_csv('output/sbs_qc_check_cos_snv.csv', index=False)

    