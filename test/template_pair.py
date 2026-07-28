import random

def gen_template(motif_len: int):
    alphabet = 'ACGT'
    # Generate a random motif with the given length
    pattern = list()
    for i in range(motif_len):
        pattern.append(alphabet[random.randint(0, 3)])
    pattern = ''.join(pattern)
    return pattern


def add_mutation(text: str, mutation_rate: float):
    alphabet = 'ACGT'
    text_len = len(text)
    mutated_n = 0
    if mutation_rate < 1e-9: # Avoid division by zero
        number_pool = -1
    else:
        number_pool = round(1 / mutation_rate)
    result = []
    for j in range(text_len):
        # If the picked number from the pool is 1, then a mutation happens
        if number_pool == -1 or random.randint(1, number_pool) != 1:
            result.append(text[j])
        else:
            mutated_n += 1
            var_type = random.randint(0, 2)
            if var_type == 0:  # Mismatch
                v = alphabet[random.randint(0, 3)]
                while v == text[j]:
                    v = alphabet[random.randint(0, 3)]
                result.append(v)
                # print('  Mismatch at %d: %s -> %s' % (j, pattern[j], v))
            elif var_type == 1:  # Deletion
                pass
                # print('  Deletion at %d' % j)
            else:  # Insertion
                v = alphabet[random.randint(0, 3)]
                # print('  Insertion at %d: %s' % (j, v))
                result.append(v)
                result.append(text[j])
    result = ''.join(result)
    return result


def work(prefix: str):
    template = gen_template(200)
    all_units = []
    copy_n = 30
    delete_n = 10
    for i in range(copy_n):
        t = add_mutation(template, 0.05)
        all_units.append(t)
    # Generate two different sets of units to delete
    deleted1 = sorted(random.sample(range(0, copy_n+1), delete_n))
    deleted2 = sorted(random.sample(range(0, copy_n+1), delete_n))
    units1, units2 = [], []
    kept1, kept2 = [], []
    for i in range(copy_n):
        if not (i in deleted1):
            units1.append(all_units[i])
            kept1.append(i)
        if not (i in deleted2):
            units2.append(all_units[i])
            kept2.append(i)
    print('Seq1 keeps', kept1)
    print('Seq2 keeps', kept2)
    flank_l = gen_template(300)
    flank_r = gen_template(300)
    seq1 = flank_l + ''.join(units1) + flank_r
    seq2 = flank_l + ''.join(units2) + flank_r

    with open(prefix + '_1.fa', 'w') as f:
        f.write('>Simulated_1\n')
        f.write(seq1)

    # TODO: add small mutations
    with open(prefix + '_2.fa', 'w') as f:
        f.write('>Simulated_2\n')
        f.write(seq2)

    # Ground Truth
    marker1 = [True] * copy_n
    marker2 = [True] * copy_n
    for i in range(copy_n):
        if i in deleted1:
            marker1[i] = False
        if i in deleted2:
            marker2[i] = False
    print('300M', end=' ')
    for i in range(copy_n):
        if marker1[i]:
            if marker2[i]:
                print("200M", end=' ')
            else:
                print("200D", end=' ')
        else:
            if marker2[i]:
                print("200I", end=' ')
    print('300M')


if __name__ == '__main__':
    work('toy')



