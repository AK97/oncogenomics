import pandas as pd
from pathlib import Path
import glob
import requests

def get_project_ids(sample_ids: pd.Series) -> pd.Series:
    '''
    Get cancer name from GDC API using Tumor_Sample_Barcode (sample id)
    Returns a Series aligned to sample ids.
    '''
    url = "https://api.gdc.cancer.gov/cases"
    CHUNK_SIZE = 500 # so as not to overload the server

    # sample_id -> case_id (TCGA-XX-YYYY)
    case_ids = (
        sample_ids.astype(str)
        .str.split("-")
        .str[:3]
        .str.join("-")
    )

    uniq_cases = pd.Series(case_ids.dropna().unique())

    case_to_pname = {}
    for i in range(0, len(uniq_cases), CHUNK_SIZE):
        chunk = uniq_cases.iloc[i:i + CHUNK_SIZE].tolist()

        payload = {
            "filters": {
                "op": "in",
                "content": {
                    "field": "submitter_id",
                    "value": chunk
                },
            },
            "format": "JSON",
            "size": len(chunk),
            "fields": "submitter_id, project.name",
            "expand": "project"
        }

        r = requests.post(url, json = payload, timeout = 120)

        hits = r.json()["data"]["hits"]
       
        for hit in hits:
            case_to_pname[hit["submitter_id"]] = hit["project"]["name"]

        
    to_return = case_ids.map(case_to_pname).rename("project_name")
    return to_return

def load_and_merge_maf_files() -> pd.DataFrame:

    PATH = 'data/raw/'
    DESIRED_COLS = [
        "Hugo_Symbol",
        "Tumor_Sample_Barcode",
        "Chromosome",
        "Start_Position",
        "Reference_Allele",
        "Tumor_Seq_Allele2",
        "Variant_Type",
        "Variant_Classification",
    ]
    
    # Check if we already have a parquet file

    parquet_filepath = "data/merged_data.parquet"
    if Path(parquet_filepath).exists():
        return pd.read_parquet(parquet_filepath)

    all_dfs = []

    for i, folder in enumerate(Path(PATH).iterdir()):
        if folder.is_dir():
            maf_file = next(folder.glob("*.maf.gz"), None)
            if maf_file:
                df = pd.read_csv(
                    maf_file, 
                    sep = '\t', 
                    comment = '#', 
                    usecols = DESIRED_COLS,
                    low_memory = False
                )
                all_dfs.append(df)

            # User feedback
            if (i + 1) % 200 == 0:
                print(f"Loaded {i + 1} files...")

    print("Merging data...")

    combined_df = pd.concat(all_dfs, ignore_index=True)

    # Add project name column
    project_series = get_project_ids(combined_df["Tumor_Sample_Barcode"])
    combined_df["project_name"] = project_series

    # Finally, keep only single-nucleotide variant type rows (drop DEL/INS/ONP/TNP)
    # NOTE: MAF files store SNP/DEL/INS under Variant_Type, not Variant_Classification.
    combined_df = combined_df[combined_df["Variant_Type"] == "SNP"]

    # Additional filtering to ensure that we only have SNPs with single nucleotide alleles
    # Also, exclude X/Y chromosomes
    mask = (
        (combined_df["Reference_Allele"].str.len() == 1) &
        (combined_df["Tumor_Seq_Allele2"].str.len() == 1) &
        (combined_df["Reference_Allele"].isin(list("ACGT"))) &
        (combined_df["Tumor_Seq_Allele2"].isin(list("ACGT"))) &
        (~combined_df["Chromosome"].isin(["chrX", "chrY"]))
    )

    combined_df = combined_df.loc[mask].copy()

    # Save as parquet
    combined_df.to_parquet(parquet_filepath)

    return combined_df
    
