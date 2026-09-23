

def main():
    print('cen', 'chm13_beg', 'chm13_end', 'chm13_len', 'hg002M_beg', 'hg002M_end', 'hg002M_len', 'diff[%]', sep='\t')
    for i in range(1, 23):
        print('cen%d' % i, end='\t')
        chm13_len = 0
        with open("chm13_cen%d.fa" % i, 'r') as f:
            s = f.readline().strip()
            assert s[0] == '>'
            r = s.split(':')[1]
            start, end = r.split('_')
            start, end = int(start), int(end)
            chm13_len = end - start
            print(start, end, chm13_len, sep='\t', end='\t')
        hg002_len = 0
        with open("hg002_Mcen%d.fa" % i, 'r') as f:
            s = f.readline().strip()
            assert s[0] == '>'
            r = s.split(':')[1]
            start, end = r.split('_')
            start, end = int(start), int(end)
            hg002_len += end - start
            print(start, end, hg002_len, sep='\t', end='\t')
        print('%.2f' % (100 * abs(hg002_len - chm13_len) / min(hg002_len, chm13_len)))


if __name__ == '__main__':
    main()