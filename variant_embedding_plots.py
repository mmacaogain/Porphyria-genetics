
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.impute import SimpleImputer
from sklearn.manifold import TSNE
import umap.umap_ as umap

# Load data
file_path = 'PCA_gpt.csv.xlsx'  # Ensure this file is in the same folder
df = pd.read_excel(file_path, sheet_name='Sheet1')

# Separate features and labels
features = df.drop(columns=["Gene"])
labels = df["Gene"]

# Impute missing values
imputer = SimpleImputer(strategy="mean")
features_imputed = imputer.fit_transform(features)

# PCA
pca = PCA(n_components=2)
pca_result = pca.fit_transform(features_imputed)
pca_df = pd.DataFrame(pca_result, columns=["PC1", "PC2"])
pca_df["Gene"] = labels

plt.figure(figsize=(8, 6))
sns.scatterplot(data=pca_df, x="PC1", y="PC2", hue="Gene", s=100, marker='o', edgecolor='k')
plt.title("PCA of Variant Prediction Scores Colored by Gene")
plt.grid(True)
plt.tight_layout()
plt.savefig("PCA_plot.png")
plt.close()

# t-SNE
tsne = TSNE(n_components=2, random_state=42, perplexity=20, n_iter=1000)
tsne_result = tsne.fit_transform(features_imputed)
tsne_df = pd.DataFrame(tsne_result, columns=["TSNE1", "TSNE2"])
tsne_df["Gene"] = labels

plt.figure(figsize=(8, 6))
sns.scatterplot(data=tsne_df, x="TSNE1", y="TSNE2", hue="Gene", s=100, marker='o', edgecolor='k')
plt.title("t-SNE of Variant Prediction Scores Colored by Gene")
plt.grid(True)
plt.tight_layout()
plt.savefig("tSNE_plot.png")
plt.close()

# UMAP
reducer = umap.UMAP(random_state=42, n_neighbors=15, min_dist=0.1)
umap_result = reducer.fit_transform(features_imputed)
umap_df = pd.DataFrame(umap_result, columns=["UMAP1", "UMAP2"])
umap_df["Gene"] = labels

plt.figure(figsize=(8, 6))
sns.scatterplot(data=umap_df, x="UMAP1", y="UMAP2", hue="Gene", s=100, marker='o', edgecolor='k')
plt.title("UMAP of Variant Prediction Scores Colored by Gene")
plt.grid(True)
plt.tight_layout()
plt.savefig("UMAP_plot.png")
plt.close()

print("Plots saved as: PCA_plot.png, tSNE_plot.png, UMAP_plot.png")
