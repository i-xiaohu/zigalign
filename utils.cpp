//
// Created by ixiaohu on 2026/1/20.
//

#include <fstream>
#include <cassert>
#include <sstream>
#include "utils.h"

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

pair<string,string> input_fasta_seq(const char *fn)
{
	gzFile f = gzopen(fn, "r");
	assert(f != nullptr);
	kseq_t *ks = kseq_init(f);
	// Only process the first sequence
	string name, seq;
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
	return make_pair(name, seq);
}

TestEntity input_csv_test_seq(int n, const char *fn)
{
	ifstream in(fn);
	assert(in.is_open());
	string line;
	getline(in, line); // Header
	for (int i = 0; i < n-1; i++) {
		getline(in, line);
	}
	getline(in, line);
	for (char &a : line) {
		if (a == ',') {
			a = ' ';
		}
	}
	stringstream ss(line);
	int id; ss >> id;
	string motif; ss >> motif;
	int period; ss >> period;
	int mutation; ss >> mutation;
	int flank_l; ss >> flank_l;
	int flank_r; ss >> flank_r;
	string seq; ss >> seq;
	in.close();
	return TestEntity{motif, period, mutation, flank_l, flank_r, seq};
}
