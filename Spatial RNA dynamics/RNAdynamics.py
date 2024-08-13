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

adata = loomdata.copy()


scv.pp.filter_and_normalize(adata,min_shared_counts=20, enforce=True)#
scv.pp.moments(adata,n_pcs=30, n_neighbors=30)#
scv.tl.recover_dynamics(adata)
scv.tl.velocity(adata, mode = "steady_state")
scv.tl.velocity_graph(adata)
scv.tl.latent_time(adata)

ident_colours = {'C0':'#8497B0', 'C1':'#878787', 'C2':'#F0CE58', 'C3':'#EB545C','C4':'#D7EF9B', 'C5':'#DCEAF7', 'C6':'#FFEE00','C7':'#00475F',
'C8':'#F297A7','C9':'#ABDDDE','C10':'#97CCF6','C11':'#B487B7','C12':'#1D65A6','C13':'#F2A104','C14':'#FC6FCF','C15':'#0FFFFF',
'C16':'#289E92','C17':'#C3FF00','C18':'#BC589B','C19':'#8ED973'}
scv.pl.velocity_embedding_stream(adata, basis='X_umap', save='Sample_Name_umap_all_arrow2.pdf' ,color = "clusters",palette = ident_colours)


scv.pl.velocity_embedding_stream(adata, figsize = (8,6),basis='X_umap', save='Sample_Name_UMAP.svg' ,legend_fontsize = 12,color = "clusters",palette = ident_colours)
scv.pl.velocity_embedding_stream(adata, fontsize = 0, size = 200, color = "clusters", alpha= 0.8, palette = ident_colours, save='Sample_Name_spatial_UMAP_V4.svg', figsize = (10,10), arrow_color='black', density=4, basis='X_Space',legend_loc='none')

scv.tl.score_genes_cell_cycle(adata)
scv.pl.scatter(adata, color_gradients=['S_score', 'G2M_score'], figsize = (8,6), save='Sample_Name_Large_cell_cycle.svg',smooth=True, perc=[5, 95])



scv.tl.velocity_confidence(adata)
keys = 'velocity_length', 'velocity_confidence'
scv.pl.scatter(adata, c=keys, cmap='coolwarm',figsize = (8,6), save='Sample_Name_Large_confidence_all_clusters.svg',perc=[5, 95])



adata_tumor = adata[np.isin(adata.obs['clusters'], ['C6','C7','C14','C18']),]
scv.pp.neighbors(adata_tumor)
scv.tl.velocity_graph(adata_tumor)
scv.pl.velocity_embedding_stream(adata_tumor, basis='X_umap', save='Sample_Name_UMAP_Tumor.svg' ,legend_fontsize = 12,color = "clusters",palette = ident_colours)
scv.pl.velocity_embedding_stream(adata_tumor, fontsize = 0, size = 400, color = "clusters", alpha= 0.5, palette = ident_colours, save='Sample_Name_spatial_UMAP__Tumor.svg', arrow_color='black', density=4, basis='X_Space',legend_loc='none')





from matplotlib.colors import LinearSegmentedColormap

colors = ['#E0F3DB', '#CCEBC5', '#A8DDB5', '#7BCCC4', '#4EB3D3', '#2B8CBE', '#0868AC', '#084081']
n_bins = 200  # Increase this number for a smoother transition
cmap_name = "Pseudotime"
cm = LinearSegmentedColormap.from_list(cmap_name, colors, N=n_bins)

scv.pl.scatter(adata, figsize = (8,6),color='velocity_pseudotime', basis='X_umap', save='Sample_Name_large_UMAP_PseudoTime.svg', color_map=cm, perc=[2, 98], colorbar=True, rescale_color=[0,1])
scv.pl.scatter(adata, figsize = (10,10), color='velocity_pseudotime', basis='X_Space', save='Sample_Name_large_slide_Spatial_UMAP_PseudoTime.svg', color_map=cm, perc=[2, 98], size = 100, colorbar=True, rescale_color=[0,1])


scv.tl.latent_time(adata)
scv.pl.scatter(adata, color='latent_time', color_map=cm, save='Sample_Name_latent_time_all_clusters.pdf',size=80)
top_genes = adata.var['fit_likelihood'].sort_values(ascending=False).index[:300]
scv.pl.heatmap(adata, var_names=top_genes, sortby='latent_time', col_color='clusters', save='Sample_Name_latent_time_topgenes.pdf',n_convolve=100)



np.sum(adata.var['velocity_genes'])
scv.tl.rank_velocity_genes(adata,groupby='clusters')
pd.DataFrame(adata.uns['rank_velocity_genes']['names']).head()
var_names = ["DOCK8", "PTPRC",'CIITA','IKZF1',"SLC38A1", "SYK",'TCF4','ARHGAP44']
scv.pl.velocity(adata_tumor, basis='X_umap',color = "clusters", add_density=True,save='Sample_Name_Large_selected_genes.pdf',var_names=var_names, colorbar=True, ncols=2)
scv.pl.velocity(adata_tumor, basis='X_umap',color = "clusters", save='Sample_Name_Large_selected_genes_no_density.pdf',var_names=var_names, colorbar=True, ncols=2)


scv.pl.proportions(adata_tumor)

