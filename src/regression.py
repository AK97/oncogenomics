import numpy as np
import pandas as pd
from sklearn.linear_model import Ridge
from sklearn.preprocessing import OneHotEncoder

def create_linreg_model(inputs: tuple, responses: pd.DataFrame) -> Ridge:
    '''
    :param inputs: exposure values, num_snv map, cancer type map
    :type inputs: tuple of DataFrames
    :param response: pathway scores
    :type response: pd.DataFrame
    '''

    sbs_clr = inputs[0]
    n_snv = inputs[1]
    c_type = inputs[2]

    if not np.allclose(sbs_clr.mean(axis=1), 0.0, atol=1e-6):
        raise ValueError("CLR rows must be centered at 0.")

    predictors = pd.concat(
        [
            sbs_clr,
            np.log10(pd.Series(n_snv, name="n_snv")),
            pd.Series(c_type, name="cancer_type")
        ],
        axis = 1,
        join = "inner"
    )

    common_samples = predictors.index.intersection(responses.index)
    predictors = predictors.loc[common_samples]
    responses = responses.loc[common_samples]

    encoder = OneHotEncoder(
        drop = "first",
        handle_unknown = "ignore",
        sparse_output = False
    )
    c_type_encoded = encoder.fit_transform(predictors[["cancer_type"]])
    c_type_cols = encoder.get_feature_names_out(["cancer_type"])
    c_type_df = pd.DataFrame(c_type_encoded, index=predictors.index, columns=c_type_cols)

    if c_type_df.shape[1] != 3:
        raise ValueError("One-hot encoded c_type must have exactly 3 columns.")
    if not c_type_df.isin([0, 1]).all().all():
        raise ValueError("One-hot encoded c_type must contain only 0/1 values.")
    if not (c_type_df.sum(axis=1) == 0).any():
        raise ValueError("Expected samples with all-zero one-hot c_type rows.")

    predictors = pd.concat(
        [
            predictors.drop(columns=["cancer_type"]),
            c_type_df
        ],
        axis=1
    )

    print(f"X matrix: {predictors.shape}")

    model = Ridge()
    model.fit(predictors, responses)
    return model
