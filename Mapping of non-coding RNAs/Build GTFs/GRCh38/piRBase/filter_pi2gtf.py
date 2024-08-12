import argparse

#python3 filter_pi2gtf.py -o hsa.piRNA.gtf -r hsa.gold.fa -i hsa.align.bed
#python3 filter_pi2gtf.py -o mmu.piRNA.gtf -r mmu.gold.fa -i mmu.align.bed
  
ap = argparse.ArgumentParser()
ap.add_argument("-i", "--input", required=True, help="Path to input GTF file")
ap.add_argument("-r", "--ref", required=True, default="type", help="Gene type name")
ap.add_argument("-o", "--output", required=True, default="Name", help="Gene type name")
ap.add_argument("-b", "--bed", required=False, action='store_true', help="Output in BED format if specified")
args = vars(ap.parse_args())

input = args["input"]
ref = args["ref"]
output = args["output"]
whetherbed = args["bed"]

def extract_sequence_names(fasta_path):
    sequence_names = set()
    with open(fasta_path, 'r') as fasta_file:
        for line in fasta_file:
            if line.startswith('>'):
                sequence_name = line[1:].split()[0]
                sequence_names.add(sequence_name)
    return sequence_names
    
def bed_to_gtf_line(bed_line):
    fields = bed_line.strip().split()
    ref_name = fields[0]
    seq_name = fields[3]
    start = str(int(fields[1]) + 1)  # BED到GTF的位置转换
    end = fields[2]
    source = 'piRBase'
    feature_type = 'exon'
    strand = fields[5]
    attributes = f'gene_id "{seq_name}"; gene_type "piRNA";'
    gtf_line = f'{ref_name}\t{source}\t{feature_type}\t{start}\t{end}\t.\t{strand}\t.\t{attributes}\n'
    return gtf_line

def filter_bed_file(bed_path, sequence_names, output_path):
    with open(bed_path, 'r') as bed_file, open(output_path, 'w') as output_file:
        for line in bed_file:
            sequence_name = line.split()[3]
            if sequence_name in sequence_names:
                if not whetherbed:
                    line = bed_to_gtf_line(line)
                output_file.write(line)


if __name__ == '__main__':
    seq_names = extract_sequence_names(ref)
    filter_bed_file(input, seq_names, output)

