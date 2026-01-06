from gseapy import ssgsea
import csv
import pandas as pd
from pathlib import Path

def import_pathway_map() -> dict:
    hallmarks_filepath = 'data/Hallmarks.txt'

    pathway_mapping = {}
    with open(hallmarks_filepath) as file:
        reader = csv.reader(file, delimiter = '\t')
        for row in reader:
            pathway = row[0][9:] # remove "HALLMARK_"
            genes = [gene for gene in row[1:] if gene]
            pathway_mapping[pathway] = genes

    return pathway_mapping

def score_ssgsea(expression_data, pathway_map) -> pd.DataFrame:

    path = "data/ssgsea_results.parquet"
    if Path(path).exists():
        return pd.read_parquet(path)

    ssgsea_scores = ssgsea(expression_data.T, pathway_map)
    results = ssgsea_scores.res2d.pivot(index = "Name", columns = "Term", values = "ES")
    results.index.name = "sample"
    results.columns.name = "pathway"

    results.to_parquet(path)

    return results
