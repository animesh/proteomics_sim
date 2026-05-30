
import numpy as np

def cv_table(df, col, groupby="peptide_id"):
    x = df.groupby(groupby)[col].agg(["mean", "std"])

    x["cv"] = np.where(
        x["mean"] > 0,
        100 * x["std"] / x["mean"],
        np.nan
    )

    return x
