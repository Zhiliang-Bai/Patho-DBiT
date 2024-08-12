from collections import defaultdict
import argparse
import re

ap = argparse.ArgumentParser()
ap.add_argument("-i", "--input", required=True, help="Path to input GTF file")
ap.add_argument("-r", "--ref", required=True, default="type", help="Gene type name")
ap.add_argument("-o", "--output", required=True, default="Name", help="Gene type name")
args = vars(ap.parse_args())

input = args["input"]
ref = args["ref"]
output = args["output"]


def parse_fasta(file_path):
    fasta_dict = defaultdict(str)
    with open(file_path, 'r') as file:
        sequence_name = ""
        for line in file:
            line = line.strip()
            if line.startswith(">"):  # 序列名行
                sequence_name = line[1:]  # 移除‘>’获取序列名
            else:
                fasta_dict[sequence_name] += line.upper()  # 追加序列到相应的序列名
    return fasta_dict
    
    
fastaref = parse_fasta(ref)
fastain  = parse_fasta(input)



with open(output, 'w') as output_file:
    for name in fastain:
        seq_name, ref_name, positions, strand, _ = re.split(fr'[:()]+', name)
        start, end = positions.split("-")
        start = str(int(start) + 1) 
        source = 'piRBase'
        feature_type = 'exon'
        attributes = f'gene_id "{seq_name}"; gene_type "piRNA";'
        gtf_line = f'{ref_name}\t{source}\t{feature_type}\t{start}\t{end}\t.\t{strand}\t.\t{attributes}\n'
        output_file.write(gtf_line)
        
        if not (seq_name in fastaref and fastain[name] == fastaref[seq_name]):
            if len(fastain[name]) != len(fastaref[seq_name]):
                print(name)
                print(len(fastain[name]) - len(fastaref[seq_name]))
                print(fastain[name])
                print(fastaref[seq_name])
                print('--')
 
