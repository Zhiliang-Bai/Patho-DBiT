import sys
import os
import pysam as ps
import mappy as mp
from collections import defaultdict as dd

def get_bc_to_clu(bc_to_spot_fn, spot_to_clu_fn):
    bc_to_clu = dict()
    with open(bc_to_spot_fn) as bc_fp, open(spot_to_clu_fn) as spot_fp:
        spot_to_clu = dict()
        for line in spot_fp:
            spot, clu = line.strip().split()
            spot_to_clu[spot] = clu
        for line in bc_fp:
            bc, idx1, idx2 = line.strip().split()
            spot = idx1 + 'x' + idx2
            if spot not in spot_to_clu:
                continue
            clu = spot_to_clu[spot]
            bc_to_clu[bc] = clu
    return bc_to_clu

def get_spot_to_clu(spot_to_clu_fn, clu_to_name_tsv):
    spot_to_clu = dict()
    clu_to_name = dict()
    with open(spot_to_clu_fn) as spot_fp, open(clu_to_name_tsv) as clu_fp:
        for line in clu_fp:
            clu, name = line.strip().split()
            clu_to_name[clu] = name
        for line in spot_fp:
            spot, clu = line.strip().split()
            if clu not in clu_to_name:
                sys.stderr.write('Cluster not in clu_to_name.tsv: {}\n'.format(clu))
                continue
            spot_to_clu[spot] = clu_to_name[clu]
    return spot_to_clu

def get_bam_fps(bam_fn, spot_to_clu, split_bam_dir, collect_non_clu=False):
    bam_fps = dd(lambda: None)
    all_clus = set(spot_to_clu.values())
    for clu in all_clus:
        if not collect_non_clu:
            split_bam_fn = split_bam_dir + '/clu' + clu + '.bam'
            split_bam_fp = ps.AlignmentFile(split_bam_fn, 'wb', template=ps.AlignmentFile(bam_fn))
        else:
            split_bam_fn = split_bam_dir + '/non_clu' + clu + '.bam'
            split_bam_fp = ps.AlignmentFile(split_bam_fn, 'wb', template=ps.AlignmentFile(bam_fn))
        bam_fps[clu] = split_bam_fp
    return bam_fps


def get_bc_from_tag(read):
    if read.has_tag('B0'):
        bc = read.get_tag('B0')
    else:
        return None
    return bc

def get_spot_from_name(qname):
    spots = qname.rsplit('|')[-1]
    spot1, spot2 = spots.split('_')[0], spots.split('_')[1]
    spot = spot1 + 'x' + spot2
    return spot

def split_bam(bam_fn, bam_fps, spot_to_clu, bc_to_clu, bc_in_name=True, collect_non_clu=False):
    if not collect_non_clu:
        with ps.AlignmentFile(bam_fn) as bam_fp:
            for read in bam_fp:
                if read.is_unmapped:
                    # sys.stderr.write('Unmapped\t{}\n'.format(read.query_name))
                    continue
                if bc_in_name:
                    spot = get_spot_from_name(read.query_name)
                    if spot not in spot_to_clu:
                        sys.stderr.write('Unclustered\t{}\n'.format(read.query_name))
                        continue
                    clu = spot_to_clu[spot]
                else:
                    bc = get_bc_from_tag(read)
                    if bc not in bc_to_clu:
                        sys.stderr.write('Unclustered\t{}\n'.format(read.query_name))
                        continue
                    clu = bc_to_clu[bc]
                bam_fps[clu].write(read)
        return
    for to_output_clu in bam_fps:
        fp = bam_fps[to_output_clu]
        with ps.AlignmentFile(bam_fn) as bam_fp:
            for read in bam_fp:
                if read.is_unmapped:
                    # sys.stderr.write('Unmapped\t{}\n'.format(read.query_name))
                    continue
                if bc_in_name:
                    spot = get_spot_from_name(read.query_name)
                    if spot not in spot_to_clu:
                        sys.stderr.write('Unclustered\t{}\n'.format(read.query_name))
                        continue
                    clu = spot_to_clu[spot]
                else:
                    bc = get_bc_from_tag(read)
                    if bc not in bc_to_clu:
                        sys.stderr.write('Unclustered\t{}\n'.format(read.query_name))
                        continue
                    clu = bc_to_clu[bc]
                if to_output_clu != clu:
                    fp.write(read)
        fp.close()

if __name__ == '__main__':
    if len(sys.argv) != 8:
        print('Usage: python {} bam_file barcode_to_pixel.tsv pixel_to_clu_ID.tsv clu_ID_to_region.tsv split_bam_dir'.format(sys.argv[0]))
        sys.exit(1)
    bam_fn, bc_to_spot_fn, spot_to_clu_fn, clu_to_name_tsv, split_bam_dir, bc_in_name_tag, collect_non_bam = sys.argv[1:]
    bc_in_name = True
    collect_non_bam = False
    bc_to_clu = get_bc_to_clu(bc_to_spot_fn, spot_to_clu_fn)
    spot_to_clu = get_spot_to_clu(spot_to_clu_fn, clu_to_name_tsv)
    if not os.path.exists(split_bam_dir):
        os.mkdir(split_bam_dir)
    
    bam_fps = get_bam_fps(bam_fn, spot_to_clu, split_bam_dir, collect_non_bam)
    split_bam(bam_fn, bam_fps, spot_to_clu, bc_to_clu, bc_in_name, collect_non_bam)