import pandas as pd
def cv_table(df,col):
    x=df.groupby("peptide_id")[col].agg(["mean","std"])
    x["cv"]=100*x["std"]/x["mean"]
    return x
