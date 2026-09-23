
if __name__ == '__main__':
    # Unialigner's CIGAR doesn't match the original sequence?

    work_dir = '/mnt/d/ClionProjects/zigalign/pairwise/unialigner'
    cen_id_list = [str(i) for i in range(1, 23)]

    out = open(work_dir + '/variation.csv', 'w')
    out.write('cen,aligner,match,aln_rate[%d],mismatch,mut_rate[%],insertion,deletion\n')
    for cid in cen_id_list:
        with open('%s/cigar/cen%s_cigar.txt' % (work_dir, cid), 'r') as f:
            d = dict()
            line = f.readline().strip()
            num = 0
            op = ' '
            d['M'], d['X'], d['I'], d['D'] = 0, 0, 0, 0
            for c in line:
                if c.isdigit():
                    num = num * 10 + int(c)
                else:
                    op = c
                    if op not in d:
                        d[op] = 0
                    d[op] += num
                    num = 0
            M, X, I, D = d['M'], d['X'], d['I'], d['D']
            L = max(M + X + I, M + X + D)
            out.write('%s,unialigner,%d,%.2f,%d,%.2f,%d,%d\n' % (cid, M + X, 100 * (M + X) / L, X, 100 * X / (M + X), I, D))
    out.close()