import pacmap
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA


input_file = "DENV.csv"
df = pd.read_csv(input_file)
print(df.columns)
print(df.shape)

#columns = df.columns
#print(columns[1:])
data = df[df.columns[1:]]
print(data.shape)

## PacMAP model
embedding = pacmap.PaCMAP(n_components=2, n_neighbors=None, MN_ratio=0.5, FP_ratio=2.0)

#PCA
pca = PCA(n_components=2)
x_pca = pca.fit_transform(data)

X_transformed = embedding.fit_transform(data.values, init="pca")
print(X_transformed)
