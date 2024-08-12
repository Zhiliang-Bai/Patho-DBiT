import pandas as pd
import re
import argparse
import os
import subprocess
import io
import numpy as np
import tempfile
import argparse




#python3 div_clpsGTF.py -i 3.gtf -T gene_type -G gene_name -f 123
ap = argparse.ArgumentParser()
ap.add_argument("-i", "--inputfiles", required=True, help="Path to input GTF file")
ap.add_argument("-o", "--outputfile", required=True, help="Path to output GTF file")
ap.add_argument("-C", "--whetherchr", required=False, action='store_true', help="whether include chr in chr name")

args = vars(ap.parse_args())


gtf_filename = args["inputfiles"]
outputfile = args["outputfile"]
whetherchr = args["whetherchr"]


gtf_filename = gtf_filename.split(':')

gtf_data = pd.concat(
    [
        pd.read_csv(
            file_path,
            sep='\t',
            header=None,
            comment='#',
            names=['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
        )
        for file_path in gtf_filename
    ],
    ignore_index=True 
)

gtf_data['seqname'] = [re.sub('^chr','',str(attr)) for attr in gtf_data['seqname']]
if whetherchr:
    gtf_data['seqname'] = ['chr' + attr for attr in gtf_data['seqname']]



original_column_order = gtf_data.columns.tolist()
merged_df = gtf_data.groupby(['seqname', 'start', 'end', 'strand']).agg({
    'source': 'first', 
    'feature': 'first',  
    'score': 'first', 
    'frame': 'first',
    'attribute': lambda x: '--'.join(x)  # 
}).reset_index()
merged_df = merged_df[original_column_order]

merged_df.to_csv(outputfile, sep='\t', index=False, header=False)

