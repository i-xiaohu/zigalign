import gzip
import sys

from Bio import SeqIO


def main(seq_fn: str, ann_fn: str, chr_name: str, out_fn: str):
    cen_start, cen_end = -1, -1
    with open(ann_fn, 'r') as f:
        for line in f:
            cols = line.split()
            chrom, start, end = cols[0], int(cols[1]), int(cols[2])
            name, score = cols[3], int(cols[4])
            rgb = [int(i) for i in cols[-1].split(',')]
            # Pure red indicates active centromere regions; choose the leftmost and the rightmost positions
            if chrom == chr_name and rgb[0] == 153 and rgb[1] == 0 and rgb[2] == 0:
                if cen_start == -1:
                    cen_start = start
                cen_end = end
                print('Active region [%d,%d), name=%s, length=%d'
                      % (start, end, name, end - start))
    if cen_end - cen_start == 0:
        print('Found no centromere on %s' % chr_name)
        return
    print('Choose the centromere region [%d,%d), length %d' % (cen_start, cen_end, cen_end - cen_start))

    with gzip.open(seq_fn, 'rt') as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if record.id == chr_name:
                print('Chromosome %s, length: %d, cen_ratio: %.2f' %
                      (chr_name, len(record.seq), 100 * (cen_end - cen_start) / len(record.seq)))
                with open(out_fn, 'w') as f:
                    f.write('>cen_%s:%d_%d\n' % (chr_name, cen_start, cen_end))
                    i = cen_start
                    seg_len = 100
                    while i + seg_len <= cen_end:
                        segment = str(record.seq[i: i + seg_len]).upper()
                        f.write('%s\n' % segment)
                        i += seg_len
                    if i < cen_end:
                        segment = str(record.seq[i: cen_end]).upper()
                        f.write('%s\n' % segment)
                break


if __name__ == '__main__':
    if len(sys.argv) != 5:
        print('Usage: extract_centromere.py <FASTA> <BED> <Chromosome> <Output>')
        exit(1)
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])