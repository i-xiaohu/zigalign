#include <getopt.h>
#include "zigalign.h"
#include "input_fasta.h"

const string VERSION = "v0.2-023";

static int str2int(const char* s) {
    return (int)strtol(s, nullptr, 10);
}

int usage(const ZigOptions &o) {
    fprintf(stderr, "Usage: zigalign [options] seq1.fa seq2.fa > aln.paf\n");
    fprintf(stderr, "Version: %s\n", VERSION.data());
    fprintf(stderr, "  Common options:\n");
    fprintf(stderr, "    -t [INT]  number of threads\n");
    fprintf(stderr, "    -v [STR]  intermediate results prefix\n");
    fprintf(stderr, "  Scoring parameters for pairwise alignment:\n");
    fprintf(stderr, "    -A [INT]  match score [%d]\n", o.mat_score);
    fprintf(stderr, "    -B [INT]  mismatch penalty [%d]\n", o.mis_pen);
    fprintf(stderr, "    -O [INT]  open gap(indel) penalty [%d]\n", o.gap_o);
    fprintf(stderr, "    -E [INT]  extend gap penalty [%d]\n", o.gap_e);
    fprintf(stderr, "    -D [INT]  repeat unit deletion penalty [%d]\n", o.del_unit_pen);
    fprintf(stderr, "  Scoring options for self-alignment:\n");
    fprintf(stderr, "    -u [INT]  minimum repeat unit size [%d]\n", o.min_unit_size);
    fprintf(stderr, "    -d [INT]  open tandem repeat penalty [%d]\n", o.open_tr_pen);
    fprintf(stderr, "    -p [INT]  close tandem repeat penalty [%d]\n", o.close_tr_pen);
    fprintf(stderr, "    -a [INT]  match score [%d]\n", o.sa_mat_score);
    fprintf(stderr, "    -b [INT]  mismatch penalty [%d]\n", o.sa_mis_pen);
    fprintf(stderr, "    -o [INT]  open gap(indel) penalty [%d]\n", o.sa_gap_o);
    fprintf(stderr, "    -e [INT]  extend gap penalty [%d]\n", o.sa_gap_e);
    fprintf(stderr, "Note: self-alignment scoring matrix must reward more and/or \n"
                    "  penalize less than regular matrix to discover tandem repeats.\n"
                    "  Pairwise scoring matrix must penalize discrepancies harder as DSI models do.\n");
    return 1;
}


int main(int argc, char *argv[]) {
    double ctime = cputime(), rtime = realtime();
    ZigOptions opt;
    if (argc == 1) return usage(opt);
    int c;
    while ((c = getopt(argc, argv, "A:B:O:E:D:u:d:p:a:b:o:e:k:v:t:@:")) >= 0) {
        switch (c) {
        case 'A':
            opt.mat_score = abs(str2int(optarg));
            break;
        case 'B':
            opt.mis_pen = -abs(str2int(optarg));
            break;
        case 'O':
            opt.gap_o = -abs(str2int(optarg));
            break;
        case 'E':
            opt.gap_e = -abs(str2int(optarg));
            break;
        case 'D':
            opt.del_unit_pen = -abs(str2int(optarg));
            break;
        case 'u':
            opt.min_unit_size = abs(str2int(optarg));
            break;
        case 'd':
            opt.open_tr_pen = -abs(str2int(optarg));
            break;
        case 'p':
            opt.close_tr_pen = -abs(str2int(optarg));
            break;
        case 'a':
            opt.sa_mat_score = abs(str2int(optarg));
            break;
        case 'b':
            opt.sa_mis_pen = -abs(str2int(optarg));
            break;
        case 'o':
            opt.sa_gap_o = -abs(str2int(optarg));
            break;
        case 'e':
            opt.sa_gap_e = -abs(str2int(optarg));
            break;
        case 'k':
            opt.band_width = abs(str2int(optarg));
            break;
        case 'v':
            opt.log_prefix = optarg;
            break;
        case 't':
            opt.n_threads = abs(str2int(optarg));
            break;
        case '@':
            opt.dev_label = optarg;
            break;
        default:
            fprintf(stderr, "Unrecognized option `%c`\n", c);
            return 1;
        }
    }
    omp_set_num_threads(opt.n_threads);

    if (argc - optind == 2) {
        align_etr(opt, argv[optind], argv[optind+1]);
    } else {
        fprintf(stderr, "Two FASTA files are required\n");
        return 1;
    }
    fprintf(stderr, "Program finishes in %.3f CPU seconds, %.3f real seconds\n", cputime()-ctime, realtime()-rtime);
    return 0;
}