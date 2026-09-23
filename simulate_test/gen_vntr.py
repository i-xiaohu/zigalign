import random
import pandas as pd


def gen_motif(motif_len: int):
    alphabet = 'ACGT'
    pattern = list()
    for i in range(motif_len):
        pattern.append(alphabet[random.randint(0, 3)])
    pattern = ''.join(pattern)
    return pattern


def add_variation(motif: str, var_rate: float):
    alphabet = 'ACGT'
    var_n = 0
    if var_rate < 1e-9:  # Avoid division by zero
        number_pool = -1
    else:
        number_pool = round(1 / var_rate)
    result = ''
    for p in motif:
        # If the picked number from the pool is 1, then a mutation happens
        if number_pool == -1 or random.randint(1, number_pool) != 1:
            result += p
        else:
            var_n += 1
            var_type = random.randint(0, 2)
            if var_type == 0:  # Mismatch
                v = alphabet[random.randint(0, 3)]
                while v == p:
                    v = alphabet[random.randint(0, 3)]
                result += v
                # print('  Mismatch at %d: %s -> %s' % (j, pattern[j], v))
            elif var_type == 1:  # Deletion
                pass
                # print('  Deletion at %d' % j)
            else:  # Insertion
                v = alphabet[random.randint(0, 3)]
                # print('  Insertion at %d: %s' % (j, v))
                result += v
                result += p
    result = ''.join(result)
    return result


def gen_vntr(motif: str, repeat_n: int, var_rate: float):
    vntr = ''
    for i in range(repeat_n):
        vntr += add_variation(motif, var_rate)
    return vntr


def gen_test_seqs(prefix: str):
    var_rate = 0.02
    seq1, seq2 = '', ''
    inter = [100, 200, 300, 200, 300]
    units = [5, 10, 20, 40, 100]
    size1 = [10, 7, 21, 9, 20]
    size2 = [7, 10, 15, 19, 30]

    for i in range(0, 5):
        flank1 = gen_motif(inter[i])
        flank2 = add_variation(flank1, var_rate)
        seq1 += flank1
        seq2 += flank2
        mid = gen_motif(units[i])
        vntr1 = gen_vntr(mid, size1[i], var_rate)
        vntr2 = gen_vntr(mid, size2[i], var_rate)
        seq1 += vntr1
        seq2 += vntr2
    name1 = '>%s|%s' % ('_'.join([str(i) for i in units]), '_'.join([str(i) for i in size1]))
    name2 = '>%s|%s' % ('_'.join([str(i) for i in units]), '_'.join([str(i) for i in size2]))
    with open('%s_1.fa' % prefix, 'w') as f:
        f.write('%s\n' % name1)
        f.write('%s\n' % seq1)
    with open('%s_2.fa' % prefix, 'w') as f:
        f.write('%s\n' % name2)
        f.write('%s\n' % seq2)


if __name__ == '__main__':
    gen_test_seqs('complex')
