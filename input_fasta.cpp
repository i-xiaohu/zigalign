#include <cassert>
#include "input_fasta.h"

double cputime()
{
    struct rusage r;
    getrusage(RUSAGE_SELF, &r);
    return r.ru_utime.tv_sec + r.ru_stime.tv_sec + 1e-6 * (r.ru_utime.tv_usec + r.ru_stime.tv_usec);
}

double realtime()
{
    struct timeval tp;
    struct timezone tzp;
    gettimeofday(&tp, &tzp);
    return tp.tv_sec + tp.tv_usec * 1e-6;
}

FastaSeq input_fasta_seq(const char *fn)
{
    FastaSeq ret;
    gzFile f = gzopen(fn, "r");
    assert(f != nullptr);
    kseq_t *ks = kseq_init(f);
    // Only process the first sequence
    string &name = ret.name, &seq = ret.seq;
    if (kseq_read(ks) >= 0) {
        fprintf(stderr, "Input sequence %s (length=%ld) from `%s`\n", ks->name.s, ks->seq.l, fn);
        int l = ks->seq.l;
        char *s = ks->seq.s;
        seq.resize(l);
        for (int i = 0; i < l; i++) {
            s[i] = (char)toupper(s[i]);
            if (s[i] !=  'A' and s[i] != 'C' and s[i] != 'G' and s[i] != 'T') {
                fprintf(stderr, "Unsupported character `%c`\n", s[i]);
                abort();
            }
            seq[i] = s[i];
        }
        l = ks->name.l;
        s = ks->name.s;
        name.resize(l);
        for (int i = 0; i < l; i++) {
            name[i] = s[i];
        }
    } else {
        fprintf(stderr, "No sequence found in `%s`\n", fn);
        abort();
    }
    kseq_destroy(ks);
    gzclose(f);
    return ret;
}
