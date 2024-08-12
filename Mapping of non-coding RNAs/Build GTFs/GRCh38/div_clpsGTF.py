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
ap.add_argument("-i", "--input", required=True, help="Path to input GTF file")
ap.add_argument("-f", "--outputfolder", required=True, help="Path to output GTF file")
ap.add_argument("-T", "--gene_typeName", required=False, default="type", help="Gene type name")
ap.add_argument("-G", "--geneName", required=False, default="Name", help="Gene type name")
args = vars(ap.parse_args())


gtf_filename = args["input"]
outputfolder = args["outputfolder"]
gene_typeName = args["gene_typeName"]
geneName = args["geneName"]

outputfolder = re.sub('/$','',outputfolder) + '/'
if not os.path.exists(outputfolder):
    os.makedirs(outputfolder)




def collapseGTF(gtfin,outgtfname):
    def uniquenames(names):
        name_counts = {} 
        outnames = []
        for name in names:
            if name in name_counts:
                name2 = f"{name}_{name_counts[name]}"
                name_counts[name] += 1
            else:
                name_counts[name] = 1
                name2 = name
            outnames.append(name2)
        return outnames

    class UnionFind:
        def __init__(self, dfin):
          genes = set(np.concatenate((dfin.iloc[:,0], dfin.iloc[:,1])))
          genes = list(genes)
          self.parent = {gene: gene for gene in genes}
          self.rank = {genes[genei]: genei for genei in range(len(genes))}
          for genei in range(dfin.shape[0]):
            gene1 = dfin.iloc[genei,0]
            gene2 = dfin.iloc[genei,1]
            root1, root2 = self.find(gene1), self.find(gene2)
            if self.rank[root1] < self.rank[root2]:
              self.parent[root1] = root2
            else:
              self.parent[root2] = root1
          groups = {}
          for gene in genes:
            root = self.find(gene)
            if root not in groups:
                groups[root] = []
            groups[root].append(gene)
          self.groups = list(groups.values())
        def find(self, gene):
            if self.parent[gene] != gene:
                self.parent[gene] = self.find(self.parent[gene])
            return self.parent[gene]

    gtfin.iloc[:,8] = uniquenames(gtfin.iloc[:,8])

    with tempfile.NamedTemporaryFile(mode='w+', delete=True) as tmpfile:
        gtfin.to_csv(tmpfile.name, sep="\t", index=False, header=False)
        command = "bedtools intersect -a " + tmpfile.name + " -b " + tmpfile.name + " -bed -wao"
        output = subprocess.check_output(command, shell=True, universal_newlines=True)


    output_temp = io.StringIO(output)
    df = pd.read_csv(output_temp, sep="\t", header=None)
    groupout = []
    if not df.empty:
    
        lengthf = df.iloc[:,4]-df.iloc[:,3]+1
        lengthb = df.iloc[:,13]-df.iloc[:,12]+1
        efficentlen = np.where(lengthf < lengthb, lengthf, lengthb)
        whetherconnect = np.where(df.iloc[:,18] > efficentlen*0.8, True, False)
        df = df.loc[whetherconnect,:]
        df = df[df.iloc[:,8] != df.iloc[:,17]]
        df = df.iloc[:,[8,17]]
        if not df.empty:
            deletrows = []
            groups = UnionFind(df).groups
            for groupi in groups:
                tempin = gtfin[gtfin.iloc[:,8].isin(groupi)]
                fis = np.where(gtfin.iloc[:,8].isin(groupi))[0]
                gtfinmin = np.min(tempin.iloc[:,3])
                gtfinmax = np.max(tempin.iloc[:,4])
                gtfin.iloc[fis[0],3] = gtfinmin
                gtfin.iloc[fis[0],4] = gtfinmax
                groupout.append("\t".join(groupi))
                deletrows.extend(list(fis[1:]))
            outgtf = gtfin.drop(gtfin.index[deletrows], axis='index')
            outgtf.to_csv(outgtfname, sep="\t", index=False, header=False)
    if df.empty:
        gtfin.to_csv(outgtfname, sep="\t", index=False, header=False)
    return groupout




gtf_data = pd.read_csv(gtf_filename, sep='\t', header=None, comment='#', names=['seqname', 'source', 'feature', 'start', 'end','score', 'strand', 'frame', 'attribute'])
gtf_data['attribute'] = [';' + attr for attr in gtf_data['attribute']]
gene_types = gtf_data['attribute'].apply(lambda x: re.sub(r'[\'"]', '',re.sub(r'[=\s;]+.*$','',re.sub(fr'^.*[=\s;]+{gene_typeName}[=\s;]+', '', x))))
gene_names = gtf_data['attribute'].apply(lambda x: re.sub(r'[\'"]', '',re.sub(r'[=\s;]+.*$','',re.sub(fr'^.*[=\s;]+{geneName}[=\s;]+', '', x))))

gtf_data['attribute'] = gene_names

groupouts = []
for typei in set(gene_types):
    tempfilename = f"{outputfolder}/{typei}.gtf"
    tempframe = gtf_data[gene_types == typei]
    indices = tempframe.index
    tempframe.loc[indices, 'attribute'] = 'gene_type "' + typei + '"; gene_name "' + tempframe.loc[indices, 'attribute'].astype(str)
    groupout = collapseGTF(tempframe,tempfilename)
    groupouts.extend(groupout)


with open(f"{outputfolder}/group.txt", "w", encoding="utf-8") as file:
    for string in groupouts:
        file.write(str(string) + "\n")
