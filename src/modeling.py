import numpy as np
import pandas as pd


def clr(sbs_data: pd.DataFrame) -> pd.DataFrame:
    EPSILON = 1e-8
    log_values = np.log(sbs_data + EPSILON)
    row_log_means = log_values.mean(axis=1)
    return log_values.sub(row_log_means, axis=0)

def standardize_pathways(pathway_data: pd.DataFrame) -> pd.DataFrame:
    mean = pathway_data.mean(axis=0)
    std = pathway_data.std(axis=0)
    return (pathway_data - mean) / std

