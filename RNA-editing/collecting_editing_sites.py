import os
import sys
import argparse
from collections import defaultdict as dd
import cgranges as cr

min_total_read_count = 10
# min_editing_ratio = 0.1
min_edited_read_count = 0


def get_spot_to_clu(spot_to_clu_fn, clu_to_name_tsv):
    spot_to_clu = dd(lambda: 'NotClustered')
    clu_to_name = dict()
    with open(spot_to_clu_fn) as spot_fp, open(clu_to_name_tsv) as clu_fp:
        for line in clu_fp:
            ele = line.strip().rsplit()
            if len(ele) < 2:
                continue
            clu, name = ele[0], ele[1]
            clu_to_name[clu] = name
        for line in spot_fp:
            ele = line.strip().rsplit()
            if len(ele) < 2:
                continue
            spot, clu = ele[0], ele[1]
            if clu not in clu_to_name:
                name = 'NoClusterName'
                # sys.stderr.write('Cluster not in clu_to_name.tsv: {}\n'.format(clu))
                # continue
            else:
                name = clu_to_name[clu]
            spot_to_clu[spot] = name
    return spot_to_clu

# collect coordiantes of each gene, crate cgranges object for each gene
def get_gene_cr(anno_gtf):
    gene_cr = cr.cgranges()
    gene_idx = dict()
    n_gene = 0
    with open(anno_gtf) as fp:
        for line in fp:
            if line.startswith('#'):
                continue
            ele = line.strip().rsplit('\t')
            if ele[2] != 'gene':
                continue
            chrom, start, end = ele[0], int(ele[3]), int(ele[4])
            gene_id, gene_name = '', ''
            if 'gene_id' in ele[8]:
                gene_id = ele[8][ele[8].index('gene_id') + 9:].split('\"')[0]
            if 'gene_name' in ele[8]:
                gene_name = ele[8][ele[8].index('gene_name') + 11:].split('\"')[0]
            if gene_name != '':
                gene_id = gene_name
            if gene_id != '':
                gene_cr.add(chrom, start-1, end, n_gene)
                gene_idx[n_gene] = gene_id
                n_gene += 1
    gene_cr.index()
    return gene_cr, gene_idx

def get_editing_sites(site_file):
    #Region	Position	Ref	Ed	Strand	db	type	dbsnp	repeat
    editing_sites = dd(lambda: dd(lambda: '')) # {(chr,pos) -> {ref/alt/strand/type}}
    with open(site_file) as f:
        header = True
        header_idx = {}
        for line in f:
            ele = line.strip().rsplit('\t')
            if  header:
                header = False
                header_idx = {ele[i]:i for i in range(len(ele))}
                continue
            site = ele[header_idx['#Region']] + ':' + ele[header_idx['Position']]
            editing_sites[site]['ref'] = ele[header_idx['Ref']]
            editing_sites[site]['alt'] = ele[header_idx['Ed']]
            editing_sites[site]['strand'] = ele[header_idx['Strand']]
            editing_sites[site]['type'] = ele[header_idx['type']]
    return editing_sites

# spot_to_editing: {spot -> {site -> count/ratio}}
def collect_editing_count(editing_count_fn, min_read_count, min_edit_count, 
                          site_wise_spot_to_editing, gene_wise_spot_to_editing,
                          site_to_gene,
                          gene_cr, gene_idx):
    with open(editing_count_fn) as f: # header: chr/pos/ref/total_count/bases/base_quality/read_names
        for line in f:
            ele = line.strip().rsplit('\t')
            chrom, pos, ref, total_count1, bases, base_quality, read_names = ele
            read_names = read_names.rsplit(',')
            if total_count1 == '0':
                continue
            ref, bases = ref.upper(), bases.upper()
            if ref == 'A': # search G in bases
                # if 'G' not in bases:
                    # continue
                target = 'G'
            elif ref == 'T':
                # if 'C' not in bases:
                    # continue
                target = 'C'
            else:
                print('Error: ref not A or T')
                sys.exit(1)
            if len(bases) != len(read_names):
                print('Error: bases and read_names not match')
                sys.exit(1)
            ref_count = bases.count('.')
            alt_count = bases.count(target)
            total_count = ref_count + alt_count
            if total_count < min_read_count or alt_count < min_edit_count:
                continue
            genes = [gene_idx[r[2]] for r in gene_cr.overlap(chrom, int(pos)-1, int(pos))]
            site = chrom + ':' + pos
            # collect site/gene wise editing count
            for base, read_name in zip(bases, read_names):
                spot = read_name.rsplit('|')[-1].replace('_', 'x')
                site_wise_spot_to_editing[spot][site]['total_count'] += 1
                for gene in genes:
                    gene_wise_spot_to_editing[spot][gene]['total_count'] += 1
                if base == target:
                    site_wise_spot_to_editing[spot][site]['ed_count'] += 1
                    for gene in genes:
                        gene_wise_spot_to_editing[spot][gene]['ed_count'] += 1
                        site_to_gene[site] = gene
    return site_wise_spot_to_editing, gene_wise_spot_to_editing

def write_spot_stats(out_dir, site_wise_spot_to_editing, gene_wise_spot_to_editing, spot_to_clu):
    ed_site_count_fn = os.path.join(out_dir, 'spot_site_editing_count.txt')
    ed_site_ratio_fn = os.path.join(out_dir, 'spot_site_editing_ratio.txt')
    ed_gene_count_fn = os.path.join(out_dir, 'spot_gene_editing_count.txt')
    ed_gene_ratio_fn = os.path.join(out_dir, 'spot_gene_editing_ratio.txt')
    stat_fn = os.path.join(out_dir, 'spot_editing_stat.txt')
    all_sites = set()
    all_genes = set()
    spot_to_stat = dd(lambda: dd(lambda: 0.0)) # {spot -> n_edit_site/gene/read/ratio}

    for spot in site_wise_spot_to_editing:
        all_sites.update(site_wise_spot_to_editing[spot].keys())
        all_genes.update(gene_wise_spot_to_editing[spot].keys())
    
    with open(ed_site_count_fn, 'w') as ed_site_cnt_fp, open(ed_site_ratio_fn, 'w') as ed_site_ratio_fp, \
         open(ed_gene_count_fn, 'w') as ed_gene_cnt_fp, open(ed_gene_ratio_fn, 'w') as ed_gene_ratio_fp:
        ed_site_cnt_fp.write('spot\t{}\n'.format('\t'.join(all_sites)))
        ed_site_ratio_fp.write('spot\t{}\n'.format('\t'.join(all_sites)))
        ed_gene_cnt_fp.write('spot\t{}\n'.format('\t'.join(all_genes)))
        ed_gene_ratio_fp.write('spot\t{}\n'.format('\t'.join(all_genes)))
        for spot in list(site_wise_spot_to_editing.keys()):
            ed_site_cnt_fp.write('{}'.format(spot))
            ed_site_ratio_fp.write('{}'.format(spot))
            ed_gene_cnt_fp.write('{}'.format(spot))
            ed_gene_ratio_fp.write('{}'.format(spot))
            for site in all_sites:
                edit_count = site_wise_spot_to_editing[spot][site]['ed_count']
                total_count = site_wise_spot_to_editing[spot][site]['total_count']
                edit_ratio = 0 if total_count == 0 else edit_count / total_count
                ed_site_cnt_fp.write('\t{}'.format(edit_count))
                ed_site_ratio_fp.write('\t{:.4f}'.format(edit_ratio))
                spot_to_stat[spot]['n_edit_site'] += 1 if edit_count > 0 else 0
                spot_to_stat[spot]['n_edit_read'] += edit_count
                spot_to_stat[spot]['n_total_read'] += total_count
            ed_site_cnt_fp.write('\n')
            ed_site_ratio_fp.write('\n')
            
            for gene in all_genes:
                edit_count = gene_wise_spot_to_editing[spot][gene]['ed_count']
                total_count = gene_wise_spot_to_editing[spot][gene]['total_count']
                edit_ratio = 0 if total_count == 0 else edit_count / total_count
                ed_gene_cnt_fp.write('\t{}'.format(edit_count))
                ed_gene_ratio_fp.write('\t{:.4f}'.format(edit_ratio))
                spot_to_stat[spot]['n_edit_gene'] += 1 if edit_count > 0 else 0
            ed_gene_cnt_fp.write('\n')
            ed_gene_ratio_fp.write('\n')
    with open(stat_fn, 'w') as fp:
        fp.write('spot\tn_edit_sites\tn_edit_genes\tn_edit_reads\tedit_ratio\n')
        for spot in spot_to_stat:
            n_edit_read, n_total_read = spot_to_stat[spot]['n_edit_read'], spot_to_stat[spot]['n_total_read']
            edit_ratio = 0 if n_total_read == 0 else n_edit_read / n_total_read
            fp.write('{}\t{}\t{}\t{}\t{:.4f}\n'.format(spot, spot_to_stat[spot]['n_edit_site'], spot_to_stat[spot]['n_edit_gene'], n_edit_read, edit_ratio))

def write_clu_stats(out_dir, site_wise_spot_to_editing, gene_wise_spot_to_editing, spot_to_clu):
    ed_site_count_fn = os.path.join(out_dir, 'clu_site_editing_count.txt')
    ed_site_ratio_fn = os.path.join(out_dir, 'clu_site_editing_ratio.txt')
    ed_gene_count_fn = os.path.join(out_dir, 'clu_gene_editing_count.txt')
    ed_gene_ratio_fn = os.path.join(out_dir, 'clu_gene_editing_ratio.txt')
    stat_fn = os.path.join(out_dir, 'clu_editing_stat.txt')
    all_sites = set()
    all_genes = set()
    all_clus = set(spot_to_clu.values())
    clu_to_stat = dd(lambda: dd(lambda: 0.0)) # {clu -> n_edit_site/gene/read/ratio}
    for spot in site_wise_spot_to_editing:
        all_sites.update(site_wise_spot_to_editing[spot].keys())
        all_genes.update(gene_wise_spot_to_editing[spot].keys())
    # calculate cluster-wise editing count/ratio
    site_wise_clu_to_editing, gene_wise_clu_to_editing = dd(lambda: dd(lambda: dd(lambda: 0.0))), dd(lambda: dd(lambda: dd(lambda: 0.0)))
    for spot in site_wise_spot_to_editing:
        if spot not in spot_to_clu:
            continue
        clu = spot_to_clu[spot]
        for site in site_wise_spot_to_editing[spot]:
            ed_count, total_count = site_wise_spot_to_editing[spot][site]['ed_count'], site_wise_spot_to_editing[spot][site]['total_count']
            site_wise_clu_to_editing[clu][site]['ed_count'] += ed_count
            site_wise_clu_to_editing[clu][site]['total_count'] += total_count
        for gene in gene_wise_spot_to_editing[spot]:
            ed_count, total_count = gene_wise_spot_to_editing[spot][gene]['ed_count'], gene_wise_spot_to_editing[spot][gene]['total_count']
            gene_wise_clu_to_editing[clu][gene]['ed_count'] += ed_count
            gene_wise_clu_to_editing[clu][gene]['total_count'] += total_count
    
    with open(ed_site_count_fn, 'w') as site_cnt_fp, open(ed_site_ratio_fn, 'w') as site_ratio_fp, \
         open(ed_gene_count_fn, 'w') as gene_cnt_fp, open(ed_gene_ratio_fn, 'w') as gene_ratio_fp:
        site_cnt_fp.write('clu\t{}\n'.format('\t'.join(all_sites)))
        site_ratio_fp.write('clu\t{}\n'.format('\t'.join(all_sites)))

        for clu in all_clus:
            site_cnt_fp.write('{}'.format(clu))
            site_ratio_fp.write('{}'.format(clu))
            for site in all_sites:
                edit_count = site_wise_clu_to_editing[clu][site]['ed_count']
                total_count = site_wise_clu_to_editing[clu][site]['total_count']
                edit_ratio = 0 if total_count == 0 else edit_count / total_count
                site_cnt_fp.write('\t{}'.format(edit_count))
                site_ratio_fp.write('\t{:.4f}'.format(edit_ratio))
                clu_to_stat[clu]['n_edit_site'] += 1 if edit_count > 0 else 0
                clu_to_stat[clu]['n_edit_read'] += edit_count
                clu_to_stat[clu]['n_total_read'] += total_count
            site_cnt_fp.write('\n')
            site_ratio_fp.write('\n')
        
        gene_cnt_fp.write('clu\t{}\n'.format('\t'.join(all_genes)))
        gene_ratio_fp.write('clu\t{}\n'.format('\t'.join(all_genes)))
        for clu in all_clus:
            gene_cnt_fp.write('{}'.format(clu))
            gene_ratio_fp.write('{}'.format(clu))
            for gene in all_genes:
                edit_count = gene_wise_clu_to_editing[clu][gene]['ed_count']
                total_count = gene_wise_clu_to_editing[clu][gene]['total_count']
                edit_ratio = 0 if total_count == 0 else edit_count / total_count
                gene_cnt_fp.write('\t{}'.format(edit_count))
                gene_ratio_fp.write('\t{:.4f}'.format(edit_ratio))
                clu_to_stat[clu]['n_edit_gene'] += 1 if edit_count > 0 else 0
            gene_cnt_fp.write('\n')
            gene_ratio_fp.write('\n')
    with open(stat_fn, 'w') as fp:
        fp.write('clu\tn_edit_sites\tn_edit_genes\tn_edit_reads\tedit_ratio\n')
        for clu in clu_to_stat:
            n_edit_read, n_total_read = clu_to_stat[clu]['n_edit_read'], clu_to_stat[clu]['n_total_read']
            edit_ratio = 0 if n_total_read == 0 else n_edit_read / n_total_read
            fp.write('{}\t{}\t{}\t{}\t{:.4f}\n'.format(clu, clu_to_stat[clu]['n_edit_site'], clu_to_stat[clu]['n_edit_gene'], n_edit_read, edit_ratio))

def write_site_stats(out_dir, site_to_gene, site_wise_spot_to_editing):
    site_stat_fn = os.path.join(out_dir, 'site_editing_count_ratio.txt')
    site_to_count = dd(lambda: dd(lambda: 0)) # {site -> {total_count/ed_count}}
    for spot in site_wise_spot_to_editing:
        for site in site_wise_spot_to_editing[spot]:
            site_to_count[site]['total_count'] += site_wise_spot_to_editing[spot][site]['total_count']
            site_to_count[site]['ed_count'] += site_wise_spot_to_editing[spot][site]['ed_count']
    with open(site_stat_fn, 'w') as out_fp:
        out_fp.write('site\ttotal_count\ted_count\tedit_ratio\n')
        for site in site_to_count:
            ed_count, total_count = site_to_count[site]['ed_count'], site_to_count[site]['total_count']
            edit_ratio = 0 if total_count == 0 else ed_count / total_count
            gene = site_to_gene[site] if site in site_to_gene else 'NA'
            out_fp.write('{}:{}\t{}\t{}\t{:.4f}\n'.format(gene, site, int(total_count), int(ed_count), edit_ratio))
    
def parser_argv():
    # parse command line arguments
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('editing_site', metavar='edit_tsv', type=str, help='editing site file')
    parser.add_argument('pixel_to_clu_ID_tsv', metavar='pixel_to_clu_ID.tsv', type=str, help='Spatial pixel to cluster ID')
    parser.add_argument('clu_ID_to_region_tsv', metavar='clu_ID_to_region.tsv', type=str, help='Cluster to regioni name')
    parser.add_argument('anno_gtf', metavar='anno.gtf', type=str, help='Annotation file')
    parser.add_argument('output_dir', metavar='out_dir', type=str, help='Output directory')
    
    option_args = parser.add_argument_group('optional arguments')
    option_args.add_argument('-m', '--min-read-count', type=int, default=min_total_read_count, help='Minimum total read count')
    option_args.add_argument('-c', '--min-edit-count', type=int, default=min_edited_read_count, help='Minimum edited read count')
    return parser.parse_args()

if __name__ == '__main__':
    args = parser_argv()
    edit_site_tsv, spot_to_clu_tsv, clu_to_name_tsv, out_dir = args.editing_site, args.pixel_to_clu_ID_tsv, args.clu_ID_to_region_tsv, args.output_dir
    anno_gtf = args.anno_gtf
    min_read_count, min_edit_count = args.min_read_count, args.min_edit_count
    # use_all_spot = args.use_all_spot
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    
    # read editing sites
    # editing_sites = get_editing_sites(site_file)
    # read spot to clu
    spot_to_clu = get_spot_to_clu(spot_to_clu_tsv, clu_to_name_tsv)
    gene_cr, gene_idx = get_gene_cr(anno_gtf)
    
    # read editing count
    site_wise_spot_to_editing = dd(lambda: dd(lambda: dd(lambda: 0.0))) # {spot -> {site -> total_count/ed_count/ratio}}
    # site_wise_clu_to_editing = dd(lambda: dd(lambda: dd(lambda: 0.0))) # {clu -> {site -> total_count/ed_count/ratio}}
    gene_wise_spot_to_editing = dd(lambda: dd(lambda: dd(lambda: 0.0))) # {spot -> {gene -> total_count/ed_count/ratio}}
    site_to_gene = dict()
    # gene_wise_clu_to_editing = dd(lambda: dd(lambda: dd(lambda: 0.0))) # {clu -> {gene -> total_count/ed_count/ratio}}
    collect_editing_count(edit_site_tsv,
                          min_read_count, min_edit_count,
                          site_wise_spot_to_editing, 
                          gene_wise_spot_to_editing, 
                          site_to_gene,
                          gene_cr, gene_idx)
    # write combined out
    write_spot_stats(out_dir, site_wise_spot_to_editing, gene_wise_spot_to_editing, spot_to_clu)
    write_clu_stats(out_dir, site_wise_spot_to_editing, gene_wise_spot_to_editing, spot_to_clu)
    write_site_stats(out_dir, site_to_gene, site_wise_spot_to_editing)    