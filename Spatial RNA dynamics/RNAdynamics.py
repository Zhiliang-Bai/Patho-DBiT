import scvelo as scv
import pandas as pd
import numpy as np
import os
import re
import scanpy as sc
import anndata as ad
import scipy.sparse as sp


# input matrix
# exon.tsv: gene-pixel count matrix which mapped to exon region of genes
# intron.tsv: gene-pixel count matrix which mapped to intron region of genes 
exonin = pd.read_csv('exon.tsv', sep='\t', index_col=0)
intronin = pd.read_csv('intron.tsv', sep='\t', index_col=0)
#np.unique(intronin.columns == exonin.columns)





X_umap = pd.read_csv("UMAP_embeddings.csv", header = 0, names=['pixel','UMAP_1','UMAP_2'])
# An example of UMAP_embeddings.csv
#"","UMAP_1","UMAP_2"
#"100x1",-0.119082775087174,1.58331056589699
#"100x2",-5.7310330682705,0.902558437534336
#"100x3",-2.7046613508175,-5.02175579314613
#"100x4",-2.94043849048119,-4.62671194320106
#"100x5",0.059147853284065,2.16501929516411
#"100x6",4.54619862499732,1.10318996662712
#"100x7",6.4721877759983,3.13371910805321
#"100x8",-6.07275984820824,-0.355289259008404
#"100x9",4.79389597836036,3.4553140615902
clusters = pd.read_csv("Cluster_and_spatial.id.csv", header = 0, names=['pixel','cluster'])
# An example of Cluster_and_spatial.id.csv
#"","SCT_snn_res.1.2"
#"100x1","8"
#"100x2","4"
#"100x3","15"
#"100x4","15"
#"100x5","4"
#"100x6","3"
#"100x7","18"
#"100x8","8"
#"100x9","3"
exonin = exonin[np.isin(exonin.index, clusters['pixel'])]
intronin = intronin[np.isin(intronin.index, clusters['pixel'])]
exonin_index_df     = exonin.index.to_frame(name='pixel')
X_umap_sorted   = exonin_index_df.merge(X_umap, on = "pixel").iloc[:,1:]
clusters_sorted = exonin_index_df.merge(clusters, on = "pixel").iloc[:,1:]


genecolnames = exonin.columns
exonin = sp.csr_matrix(exonin.values)
intronin = sp.csr_matrix(intronin.values)


loomdata = ad.AnnData(X=exonin)
loomdata.layers["spliced"] = exonin
loomdata.layers["unspliced"] = intronin
#loom_data.obsm['Space'] = barcodes_ordered.iloc[:,1:3].values.astype('float64')
loomdata.obsm['X_umap'] = X_umap_sorted.values
loomdata.obs['clusters'] = clusters_sorted.astype('str').values
loomdata.var.index = genecolnames
loomdata.obs.index = exonin_index_df.index

positions = pd.DataFrame(loomdata.obs.index)[0].str.split('x', expand=True)
positions = positions.astype('int')
aa = positions.copy()
aa.iloc[:,1:2] = 101 - aa.iloc[:,1:2]
#aa = (aa-25)/100
loomdata.obsm['X_Space'] = aa.values.astype('float64')

adata2 = loomdata.copy()

scv.pp.filter_and_normalize(adata2,min_shared_counts=20,enforce=True)#
scv.pp.moments(adata2,n_pcs=30, n_neighbors=30)#
scv.tl.recover_dynamics(adata2)
scv.tl.velocity(adata2, mode = "steady_state")
scv.tl.velocity_graph(adata2)
scv.tl.latent_time(adata2)
scv.pl.velocity_embedding_stream(adata2, basis='X_umap', save='Sample_name_umap_all_arrow2.pdf' ,color = "clusters")

adata3 = adata2[np.isin(adata2.obs['clusters'], ['6','7','14','18']),]
scv.pp.neighbors(adata3)
scv.tl.velocity_graph(adata3)
scv.pl.velocity_embedding_stream(adata3, basis='X_umap',color = "clusters", save='Sample_name_umap_tumor_arrow2.svg', figsize=(7,5), legend_fontsize = 9, density = 0.75, title='') 
scv.pl.velocity_embedding_grid(adata3, fontsize = 0, size = 300, color = "clusters", arrow_size= 2, arrow_length= 1.5, save='Sample_name_slide_tumor_grid_GridArrow2.pdf', arrow_color='black', scale = 2, basis='X_Space',legend_loc='none', smooth = 0.8,n_neighbors=8)


from matplotlib.colors import LinearSegmentedColormap

colors = ["green", "red"]
n_bins = 100  # Increase this number for a smoother transition
cmap_name = "Pseudotime"
cm = LinearSegmentedColormap.from_list(cmap_name, colors, N=n_bins)

scv.pl.scatter(adata2, color='velocity_pseudotime', basis='X_umap', save='Sample_name_umap_all_PseudoTime2.pdf', color_map=cm, perc=[2, 98], colorbar=True, rescale_color=[0,1])
scv.pl.scatter(adata2, color='velocity_pseudotime', basis='X_Space', save='Sample_name_slide_all_PseudoTime2.pdf', color_map=cm, perc=[2, 98], size = 80, colorbar=True, rescale_color=[0,1])
scv.pl.scatter(adata3, color='velocity_pseudotime', basis='X_umap', save='Sample_name_umap_tumor_PseudoTime2.pdf', color_map=cm, perc=[2, 98], colorbar=True, rescale_color=[0,1])
scv.pl.scatter(adata3, color='velocity_pseudotime', basis='X_Space', save='Sample_name_slide_tumor_PseudoTime2.pdf', color_map=cm, perc=[2, 98], size = 80, colorbar=True, rescale_color=[0,1])



# could use these commands to rank genes
#scv.tl.rank_velocity_genes(adata3,groupby='clusters')
#import pandas as pd
#pd.DataFrame(adata3.uns['rank_velocity_genes']['names']).head()
var_names = ['IGHA2','SYTL2','MALAT1','MEF2C']# genes you want to see
scv.pl.velocity(adata2, basis='X_umap',color = "clusters", var_names=var_names, colorbar=True, ncols=2, add_density = True)
