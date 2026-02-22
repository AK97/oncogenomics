from matplotlib import pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd
import os

def _resolve_output_path(output_path):
    if os.path.exists(output_path):
        base, ext = os.path.splitext(output_path)
        counter = 1
        candidate = f"{base}_{counter}{ext}"
        while os.path.exists(candidate):
            counter += 1
            candidate = f"{base}_{counter}{ext}"
        output_path = candidate
    return output_path


def plot_mut_sig_landscape(data, output_path):
    sbs_cols = [col for col in data.columns if "SBS" in col]
    if not sbs_cols:
        raise ValueError("No SBS columns found in the input data.")

    grouped_means = data.groupby("cancer_type", sort=False)[sbs_cols].mean()
    top_sbs_set = set()
    TOP_N = 8
    for _, row in grouped_means.iterrows():
        top_sbs_set.update(row.sort_values(ascending=False).head(TOP_N).index.tolist())

    selected_sbs = (
        grouped_means[list(top_sbs_set)]
        .mean()
        .sort_values(ascending=False)
        .index.tolist()
    )
    other_sbs = [col for col in sbs_cols if col not in top_sbs_set]
    grouped = grouped_means[selected_sbs].copy()
    grouped["Other"] = grouped_means[other_sbs].sum(axis=1) if other_sbs else 0.0

    colors = sns.color_palette("tab20", n_colors=len(selected_sbs))
    color_map = {sbs: color for sbs, color in zip(selected_sbs, colors)}
    color_map["Other"] = (0.6, 0.6, 0.6)

    fig, ax = plt.subplots(figsize=(max(8, len(grouped) * 0.6), 6))
    x = np.arange(len(grouped))
    bar_width = 0.85
    for i, (cancer_type, row) in enumerate(grouped.iterrows()):
        sorted_sbs = row.sort_values(ascending=False).index.tolist()
        bottom = 0.0
        for sbs in sorted_sbs:
            value = row[sbs]
            if value == 0:
                continue
            ax.bar(
                x[i],
                value,
                width=bar_width,
                bottom=bottom,
                color=color_map[sbs],
                edgecolor="none",
            )
            bottom += value

    ax.set_xticks(x, grouped.index, rotation=45, ha="right")
    handles = [
        plt.Rectangle((0, 0), 1, 1, color=color_map[sbs])
        for sbs in (["Other"] + selected_sbs)
    ]
    ax.legend(handles, ["Other"] + selected_sbs, title="SBS", bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)
    ax.set_xlabel("Cancer Type")
    ax.set_ylabel("Mean SBS Fraction Exposure")
    ax.set_title("Top 10 SBS Signatures by Cancer Type")
    plt.tight_layout()

    output_path = _resolve_output_path(output_path)
    fig.savefig(output_path, dpi=200, bbox_inches="tight", pad_inches=0.2)
    plt.close(fig)


def plot_sbs_qc_check(
    summary_csv_path="output/sbs_qc_check_summary.csv",
    cos_snv_csv_path="output/sbs_qc_check_cos_snv.csv",
    hist_output_path="figures/sbs_qc_cosine_histogram.png",
    scatter_output_path="figures/sbs_qc_cosine_vs_log_n_snv.png",
):
    common_figsize = (8, 5)

    summary_df = pd.read_csv(summary_csv_path)
    cos_snv_df = pd.read_csv(cos_snv_csv_path)

    required_summary_cols = {"cosine_similarity"}
    required_cos_snv_cols = {"cosine_similarity", "n_snv", "cancer_type"}
    if not required_summary_cols.issubset(summary_df.columns):
        missing = required_summary_cols.difference(summary_df.columns)
        raise ValueError(f"Missing columns in {summary_csv_path}: {sorted(missing)}")
    if not required_cos_snv_cols.issubset(cos_snv_df.columns):
        missing = required_cos_snv_cols.difference(cos_snv_df.columns)
        raise ValueError(f"Missing columns in {cos_snv_csv_path}: {sorted(missing)}")

    fig_hist, ax_hist = plt.subplots(figsize=common_figsize)
    sns.histplot(
        data=summary_df,
        x="cosine_similarity",
        bins=np.linspace(0.0, 1.0, 41),
        binrange=(0.0, 1.0),
        edgecolor="white",
        ax=ax_hist,
    )
    ax_hist.set_xlim(0.5, 1.0)
    ax_hist.set_title("Distribution of Cosine Similarity")
    ax_hist.set_xlabel("Cosine Similarity")
    ax_hist.set_ylabel("Sample Count")
    fig_hist.tight_layout()
    hist_output_path = _resolve_output_path(hist_output_path)
    fig_hist.savefig(hist_output_path, dpi=200)
    plt.close(fig_hist)

    fig_scatter, ax_scatter = plt.subplots(figsize=common_figsize)
    cancer_types = cos_snv_df["cancer_type"].dropna().unique().tolist()
    base_palette = ["#808080", "#000000", "#FF8C00", "#0000FF"]  # gray, black, orange, blue
    palette = {
        ctype: base_palette[i % len(base_palette)]
        for i, ctype in enumerate(cancer_types)
    }
    sns.scatterplot(
        data=cos_snv_df,
        x="n_snv",
        y="cosine_similarity",
        hue="cancer_type",
        palette=palette,
        alpha=0.8,
        s=45,
        ax=ax_scatter,
    )
    ax_scatter.set_title("Cosine Similarity vs log10(SNV Count)")
    ax_scatter.set_xlabel("log10(SNV Count)")
    ax_scatter.set_ylabel("Cosine Similarity")
    ax_scatter.legend(title="Cancer Type", loc="best", frameon=False)
    fig_scatter.tight_layout()
    scatter_output_path = _resolve_output_path(scatter_output_path)
    fig_scatter.savefig(scatter_output_path, dpi=200)
    plt.close(fig_scatter)

    return {
        "histogram_path": hist_output_path,
        "scatter_path": scatter_output_path,
    }





# Develop figure for mutational signature landscape
# aggr_data_mut_sig_landscape = pd.read_csv('output/aggr_data_mut_sig_landscape.csv', index_col=0)
# plot_mut_sig_landscape(aggr_data_mut_sig_landscape, 'figures/MutationalSignatureLandscape.png')

# Develop figure for SBS QC check
# plot_sbs_qc_check()
