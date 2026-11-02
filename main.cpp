#include <getopt.h>
#include "zigalign.h"
#include "input_fasta.h"

const string VERSION = "v0.3";

static int str2int(const char* s) {
    return (int)strtol(s, nullptr, 10);
}

int usage(const ZigOptions &o) {
    fprintf(stderr, "Usage: zigalign [options] seq1.fa seq2.fa > aln.paf\n");
    fprintf(stderr, "Version: %s\n", VERSION.data());
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "    -t  number of threads, default=%d\n", o.n_threads);
    fprintf(stderr, "    -o  output alignment in PAF format, default=stdout\n");
    fprintf(stderr, "    -c  replace M with =X in CIGAR, default=false\n");
    fprintf(stderr, "    -y  output detailed alignment in TSV format, recommended to set\n");
    fprintf(stderr, "    -S  [mat,mis,gap_o,gap_e] score matrix of SI model, default={%d,%d,%d,%d}\n",
        o.si.mat, o.si.mis, o.si.gap_o, o.si.gap_e);
    fprintf(stderr, "    -D  penalty of D transfer, default=%d\n", o.open_tr_pen);
    fprintf(stderr, "    -m  minimum repeat unit size, default=%d\n", o.min_unit_size);
    fprintf(stderr, "    -F  [mat,mis,gap_o,gap_e] score matrix of self-alignment, default={%d,%d,%d,%d}\n",
        o.self_aln.mat, o.self_aln.mis, o.self_aln.gap_o, o.self_aln.gap_e);
    fprintf(stderr, "    -U  [mat,mis,gap_o,gap_e] score matrix of pair-alignment for repeat units, default={%d,%d,%d,%d}\n",
        o.unit_aln.mat, o.unit_aln.mis, o.unit_aln.gap_o, o.unit_aln.gap_e);
    fprintf(stderr, "    -z  maximum deleted repeat units, default=%d\n", o.max_units_del);
    fprintf(stderr, "    -B  penalty of B transfer, default=%d\n", o.close_tr_pen);
    fprintf(stderr, "    -k  minimum repeat count; filtered if fewer, default=%d\n", o.min_reps_cnt);
    fprintf(stderr, "\n");
    return 1;
}

ScoreParameters parse_scores(const char *s)
{
    assert(s[0] == '[');
    string x;
    vector<int> nums;
    for (int i = 1; s[i] != ']'; i++) {
        if (s[i] == ',') {
            nums.push_back(str2int(x.data()));
            x = "";
        } else {
            x += s[i];
        }
    }
    nums.push_back(str2int(x.data()));
    assert(nums.size() == 4);
    ScoreParameters ret(nums[0], nums[1], nums[2], nums[3]);
    return ret;
}

int main(int argc, char *argv[]) {
    double ctime = cputime(), rtime = realtime();
    ZigOptions opt;
    if (argc == 1) return usage(opt);
    int c;
    while ((c = getopt(argc, argv, "t:o:cy:S:D:m:F:U:z:B:k:l:d:")) >= 0) {
        switch (c) {
        case 't':
            opt.n_threads = abs(str2int(optarg));
            break;
        case 'o':
            opt.paf_fn = optarg;
            break;
        case 'c':
            opt.cigar_neq = true;
            break;
        case 'y':
            opt.dsi_fn = optarg;
            break;
        case 'S':
            opt.si = parse_scores(optarg);
            break;
        case 'D':
            opt.open_tr_pen = str2int(optarg);
            break;
        case 'm':
            opt.min_unit_size = str2int(optarg);
            break;
        case 'F':
            opt.self_aln = parse_scores(optarg);
            break;
        case 'U':
            opt.unit_aln = parse_scores(optarg);
            break;
        case 'z':
            opt.max_units_del = str2int(optarg);
            break;
        case 'B':
            opt.close_tr_pen = str2int(optarg);
            break;
        case 'k':
            opt.min_reps_cnt = str2int(optarg);
            break;
        case 'l':
            opt.log_prefix = optarg;
            break;
        case 'd':
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