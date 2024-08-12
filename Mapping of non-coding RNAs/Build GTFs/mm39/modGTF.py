import pandas as pd
import re
import argparse
import os
import subprocess
import io
import numpy as np
import tempfile
import argparse




#python3 modGTF.py -i 3.gtf -T gene_type -G gene_name -f 123
ap = argparse.ArgumentParser()
ap.add_argument("-i", "--input", required=True, help="Path to input GTF file")
ap.add_argument("-f", "--outputfile", required=True, help="Path to output GTF file")
ap.add_argument("-T", "--gene_typeName", required=False, default="type", help="Gene type name")
ap.add_argument("-G", "--geneName", required=False, default="Name", help="Gene type name")
args = vars(ap.parse_args())


gtf_filename = args["input"]
outputfile = args["outputfile"]
gene_typeName = args["gene_typeName"]
geneName = args["geneName"]





gtf_data = pd.read_csv(gtf_filename, sep='\t', header=None, comment='#', names=['seqname', 'source', 'feature', 'start', 'end','score', 'strand', 'frame', 'attribute'])
gtf_data['attribute'] = [';' + attr for attr in gtf_data['attribute']]
gene_types = gtf_data['attribute'].apply(lambda x: re.sub(r'[\'"]', '',re.sub(r'[=\s;]+.*$','',re.sub(fr'^.*[=\s;]+{gene_typeName}[=\s;]+', '', x))))
gtf_data['gene_types'] = gene_types
gene_names = gtf_data['attribute'].apply(lambda x: re.sub(r'[\'"]', '',re.sub(r'[=\s;]+.*$','',re.sub(fr'^.*[=\s;]+{geneName}[=\s;]+', '', x))))
gtf_data['gene_names'] = gene_names





gtf_data['aa'] = gtf_data['seqname'].astype(str) + ':' + gtf_data['start'].astype(str) + ':' + gtf_data['end'].astype(str) + ':' + gtf_data['strand'].astype(str) + ':' + gtf_data['gene_names'].astype(str)
gtf_data = gtf_data.drop_duplicates(subset='aa', keep='last')
gtf_data = gtf_data.drop(columns='aa')

gtf_data = gtf_data[gtf_data['feature'].isin(['transcript', 'exon'])]
special_types = ['miRNA', 'rRNA_pseudogene', 'rRNA', 'scaRNA', 'scRNA', 'snoRNA', 'snRNA', 'vault_RNA', 'tRNA', 'misc_RNA', 'piRNA', 'Y_RNA']
special_cols = gene_types.isin(special_types)
gtf_data.loc[special_cols, 'attribute'] = gtf_data.loc[special_cols,'gene_names'] + '__' + gtf_data.loc[special_cols,'gene_types']
gtf_data.loc[~special_cols, 'attribute'] = gtf_data.loc[~special_cols,'gene_names'] + '__' + gtf_data.loc[~special_cols, 'feature']




genewithtrans = gtf_data.loc[~special_cols,'gene_names'][gtf_data.loc[~special_cols, 'feature'] == 'transcript']
genewithexon = gtf_data.loc[~special_cols,'gene_names'][gtf_data.loc[~special_cols, 'feature'] == 'exon']
genewithtrans = set(genewithtrans.unique())
genewithexon = set(genewithexon.unique())
needchange = genewithtrans - genewithexon
needchange = pd.Series(list(needchange))
gtf_data.loc[gtf_data['gene_names'].isin(needchange), 'attribute'] = gtf_data.loc[gtf_data['gene_names'].isin(needchange), 'gene_names'] + '__exon'

gtf_data = gtf_data.drop(columns='gene_types')
gtf_data = gtf_data.drop(columns='gene_names')
gtf_data.to_csv(outputfile, sep="\t", index=False, header=False)

