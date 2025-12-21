import requests
import pandas as pd
import json
from pathlib import Path
from numpy import log1p

def get_rna_seq_for_samples(samples:pd.DataFrame) -> pd.DataFrame:
    '''
    Pulls RNA-seq data from GDC.
    Pulls for only the samples we are working with (filtered, etc)
    '''

    if Path('data/rna_seq_data_for_manifest.parquet').exists():
        return pd.read_parquet('data/rna_seq_data_for_manifest.parquet')

    URL = 'https://api.gdc.cancer.gov/files'

    sample_ids_16 = samples.index.str[:16].tolist()

    def build_filters(submitter_ids):
        return {
            "op": "and",
            "content": [
                {
                    "op": "in",
                    "content": {
                        "field": "cases.samples.submitter_id",
                        "value": submitter_ids
                    }
                },
                {
                    "op": "in",
                    "content": {
                        "field": "data_category",
                        "value": ["Transcriptome Profiling"]
                    }
                },
                {
                    "op": "in",
                    "content": {
                        "field": "data_type",
                        "value": ["Gene Expression Quantification"]
                    }
                },
                {
                    "op": "in",
                    "content": {
                        "field": "analysis.workflow_type",
                        "value": ["STAR - Counts"]
                    }
                }
            ]
        }

    ## Craziness of http requesting errors... ##
    session = requests.Session()
    retries = requests.packages.urllib3.util.retry.Retry(
        total=5,
        backoff_factor=1,
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=["GET"]
    )
    adapter = requests.adapters.HTTPAdapter(max_retries=retries)
    session.mount("https://", adapter)
    session.mount("http://", adapter)

    hits = []
    chunk_size = 50
    for i in range(0, len(sample_ids_16), chunk_size):
        chunk = sample_ids_16[i:i+chunk_size]
        payload = {
            "filters": json.dumps(build_filters(chunk)),
            "format": "JSON",
            "size": len(chunk),
            "fields": ",".join([
                "file_id",
                "file_name",
                "cases.samples.submitter_id",
                "cases.project.project_id"
            ])
        }

        res = session.get(
            url = "https://api.gdc.cancer.gov/files",
            params = payload
        )
        res.raise_for_status()
        hits.extend(res.json()["data"].get("hits", []))
    
    rows = []
    for h in hits:
        file_id = h["file_id"]
        file_name = h.get("file_name")
        # Extract project + sample barcode(s)
        for case in h.get("cases", []):
            project_id = case.get("project", {}).get("project_id")
            for s in case.get("samples", []):
                sample16 = s.get("submitter_id")
                rows.append({
                    "file_id": file_id,
                    "file_name": file_name,
                    "project_id": project_id,
                    "sample16": sample16,
                })

    api_results = pd.DataFrame(rows)
    api_results.to_parquet('data/rna_seq_data_for_manifest.parquet')

    return api_results

def write_manifest(api_results) -> None:

    manifest = pd.DataFrame({
        'id': api_results['file_id'],
    })
    
    manifest.to_csv('data/rna_seq_manifest.tsv', sep='\t', index=False)


def build_counts_matrix(uuid_to_sample:pd.DataFrame) -> pd.DataFrame:
    '''
    counts = samples x genes (ints)  
    rna_meta = per-sample metadata
    '''
    directory = Path('data/raw_rna_seq')

    # Keep only the columns we need and ensure unique mapping
    uuid_to_sample = uuid_to_sample[['file_id', 'sample16']].drop_duplicates()
    file_to_sample = dict(zip(uuid_to_sample['file_id'], uuid_to_sample['sample16']))

    sample_to_counts = {}

    for folder in directory.iterdir():
        if not folder.is_dir():
            continue

        file_id = folder.name
        sample16 = file_to_sample.get(file_id)
        if sample16 is None:
            # Skip folders that don't map back to a sample
            continue

        tsv_files = list(folder.glob("*.tsv"))
        if not tsv_files:
            continue

        tsv_path = tsv_files[0]
        df = pd.read_csv(
            tsv_path,
            sep='\t',
            comment='#',
            usecols=['gene_name', 'unstranded']
        )

        # Drop rows without a gene name (e.g., summary rows)
        df = df.dropna(subset=['gene_name'])

        counts_series = pd.to_numeric(df['unstranded'], errors='coerce')
        counts_series.index = df['gene_name']

        sample_to_counts[sample16] = counts_series

    if not sample_to_counts:
        return pd.DataFrame()

    counts_matrix = pd.DataFrame(sample_to_counts).T
    counts_matrix.index.name = 'sample16'

    return counts_matrix

def normalize_and_filter_rna_seq(counts_matrix:pd.DataFrame) -> pd.DataFrame:
    '''
    Counts -> CPM (Counts per million)  
    log1p to stabilize variance  
    Normalized = CPM + log1p
    '''

    cpm = counts_matrix.div(counts_matrix.sum(axis=1), axis=0) * 1e6
    normalized_counts = log1p(cpm)

    # Filter lowly expressed genes
    to_keep = (normalized_counts >= log1p(1)).mean(axis=0) > 0.2

    return normalized_counts.loc[:, to_keep]

