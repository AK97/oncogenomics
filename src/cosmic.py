import pandas as pd
import numpy as np
from pathlib import Path
from src.channels import CHANNELS_96
from scipy.optimize import nnls

def filter_snv_counts(mat96:pd.DataFrame, min_count:int = 100) -> pd.DataFrame:
    
    n_snv = pd.Series(mat96.attrs["n_snv"], name = "n_snv")

    desired = n_snv[n_snv >= min_count].index
    mat96_filtered = mat96.loc[desired]
    return mat96_filtered

def import_cosmic_signatures() -> pd.DataFrame:
    COSMIC_DATA_PATH = 'data/COSMIC_v3.5_SBS_GRCh38.txt'

    # check for parquet file
    if Path('data/COSMIC_v3.5_SBS_GRCh38.parquet').is_file():
        return pd.read_parquet('data/COSMIC_v3.5_SBS_GRCh38.parquet')
    
    cosmic_data = pd.read_csv(COSMIC_DATA_PATH, sep='\t')

    # Align to CHANNELS_96
    realigned = cosmic_data.set_index("Type").loc[CHANNELS_96].reset_index()

    realigned.to_parquet('data/COSMIC_v3.5_SBS_GRCh38.parquet')

    return realigned

def run_nnls(mutation_matrix: pd.Series, signatures: pd.DataFrame) -> pd.Series:
    # non-negative least squares optimization
    # returns fractional exposures

    M_mat = mutation_matrix.to_numpy(dtype=float)
    S_mat = signatures.iloc[:, 1:].to_numpy(dtype=float)

    # print(M_mat.shape, S_mat.shape)

    n_samples, n_channels = M_mat.shape
    _, n_signatures = S_mat.shape

    exposures = np.zeros((n_samples, n_signatures), dtype=float)

    for i in range(n_samples):
        y = M_mat[i]

        e, _ = nnls(S_mat, y)
        exposures[i] = e

    Exp = pd.DataFrame(
        exposures,
        index = mutation_matrix.index,
        columns = signatures.columns[1:]
    )

    fractional_exposures = Exp.div(Exp.sum(axis=1), axis=0)
    return fractional_exposures