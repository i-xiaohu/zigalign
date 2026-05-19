//
// Created by ixiaohu on 2026/1/20.
//

#ifndef TRDP_UTILS_H
#define TRDP_UTILS_H

#include <string>
#include <sys/resource.h>
#include <sys/time.h>
#include <zlib.h>
#include "kseq.h"
using namespace std;

KSEQ_INIT(gzFile, gzread)

double cputime();
double realtime();

static inline int str2int(const char* s) {
	return (int)strtol(s, nullptr, 10);
}

pair<string,string> input_fasta_seq(const char *fn);

struct TestEntity {
	string motif;
	int period;
	int mutation;
	int flank_l;
	int flank_r;
	string seq;
};

TestEntity input_csv_test_seq(int n, const char *fn);

#endif //TRDP_UTILS_H
