//
// Created by ixiaohu on 2025/12/30.
//

#include <vector>
#include <cstdint>
#include <cstdio>
#include <string>
#include <fstream>
#include <cassert>
#include <iostream>
#include <getopt.h>
#include <algorithm>
#include <numeric>
#include "utils.h"
using namespace std;

// Zigalign parameters
struct ZigOptions {
	// Scoring matrix for non-repeat regions
	int mat_score;
	int mis_pen;
	int gap_o;
	int gap_e;
	// Minimum repeat unit size (to reduce random repeats)
	int min_unit_size;
	// Scoring matrix for tandem repeats in self-alignment (reward more and/or penalize less)
	int tr_mat_score;
	int tr_mis_pen;
	int tr_gap_o;
	int tr_gap_e;
	// Penalty for opening/closing tandem repeats
	int open_tr_pen;
	int close_tr_pen;
	// Penalty for deleting duplications
	int del_dup_o;
	int del_dup_e;
	// Visualization
	const char *vis_prefix;

	ZigOptions() {
		// Scoring matrix for SI model
		mat_score = 1;
		mis_pen = -4;
		gap_o = -6;
		gap_e = -1;
		// Scoring matrix for DSI model should be adjusted by duplication and variation rate
		min_unit_size = 5;
		tr_mat_score = 2;
		tr_mis_pen = -3;
		tr_gap_o = -3;
		tr_gap_e = -1;
		open_tr_pen = -2;
		close_tr_pen = -6;
		del_dup_o = -6;
		del_dup_e = -2;
		vis_prefix = nullptr;
	}
};

const int INF = 100000000;
const int NORMAL = 0;
const int START_REP = 1;
const int NEW_COPY = 2;
const int WITHIN_REP = 3;
const int END_REP = 4;

int DEBUG = 0;

struct Dp1Cell {
	int E, F, H; // The original SW matrix
	int D_gate; // Gate of duplication
	int D_end; // Where duplication ends
	int D_beg; // Where duplication or new copy begins
	int de, df, dh; // Sub matrix of duplication
	int pi, pj, event; // Backtrace

	Dp1Cell() {
		E = F = H = -INF;
		D_gate = -INF;
		D_end = -1;
		D_beg = -1;
		de = df = dh = -INF;
		pi = pj = event = -1;
	}
};

struct RepUnit {
	int tb, te; // 1-based
	int qb, qe; // [qb, qe) is a tandem repeat of [tb, te)
	int score; // Alignment score
	int match, mis, gap;

	RepUnit() {
		tb = te = 0;
		qb = qe = 0;
		score = -INF;
		match = mis = gap = 0;
	}
};

// Stage 1: identify breakpoints of tandem repeats using self alignment
static vector<RepUnit> self_alignment(const ZigOptions &o, int n, const char *seq, const string &vis_fn)
{
	double ctime = cputime();
	const int MAT_SCORE = o.mat_score;
	const int MIS_PEN = o.mis_pen;
	const int GAP_O = o.gap_o;
	const int GAP_E = o.gap_e;
	const int MIN_UNIT = o.min_unit_size;
	const int OPEN_TR = o.open_tr_pen;
	const int CLOSE_TR = o.close_tr_pen;
	const int TR_MAT_SCORE = o.tr_mat_score;
	const int TR_MIS_PEN = o.tr_mis_pen;
	const int TR_GAP_O = o.tr_gap_o;
	const int TR_GAP_E = o.tr_gap_e;

	// TODO: reduce memory consumption
	vector<vector<Dp1Cell>> dp;
	dp.resize(n + 1);
	for (int i = 0; i <= n; i++) {
		dp[i].resize(i + 1);
	}
	// Global alignment in left-down triangle
	dp[0][0].H = 0;
	for (int i = 1; i <= n; i++) {
		dp[i][0].E = dp[i][0].H = GAP_O + GAP_E * i;
		for (int j = 1; j < i; j++) {
			dp[i][j].E = max(dp[i-1][j].H + GAP_O, dp[i-1][j].E) + GAP_E;
			dp[i][j].F = max(dp[i][j-1].H + GAP_O, dp[i][j-1].F) + GAP_E;
			int M = dp[i-1][j-1].H + (seq[i-1] == seq[j-1] ? MAT_SCORE : MIS_PEN);
			if (dp[i][j].E > dp[i][j].H) {
				dp[i][j].H = dp[i][j].E;
				dp[i][j].pi = i-1;
				dp[i][j].pj = j;
				dp[i][j].event = NORMAL;
			}
			if (dp[i][j].F > dp[i][j].H) {
				dp[i][j].H = dp[i][j].F;
				dp[i][j].pi = i;
				dp[i][j].pj = j-1;
				dp[i][j].event = NORMAL;
			}
			if (M > dp[i][j].H) {
				dp[i][j].H = M;
				dp[i][j].pi = i-1;
				dp[i][j].pj = j-1;
				dp[i][j].event = NORMAL;
			}
		}
		int f = max(dp[i][i-1].H + GAP_O, dp[i][i-1].F) + GAP_E;
		int m = dp[i-1][i-1].H + MAT_SCORE;
		if (f > dp[i][i].H) {
			dp[i][i].H = f;
			dp[i][i].pi = i;
			dp[i][i].pj = i-1;
			dp[i][i].event = NORMAL;
		}
		if (m > dp[i][i].H) {
			dp[i][i].H = m;
			dp[i][i].pi = i-1;
			dp[i][i].pj = i-1;
			dp[i][i].event = NORMAL;
		}

		if (i > MIN_UNIT) {
			// D gate comes from the last row
			int max_value = dp[i-1][i-1].H; // Diagonal must be the maximum in regular matrix
			int D_end = i - 1;
			int D_beg = -1;
			int event = START_REP;
			for (int j = 1; j < i-1; j++) {
				// Only start new copy after the end of duplication to prevent illegal path
				if (dp[i-1][j].D_end == j and dp[i-1][j].dh > max_value) { // If multiple maximums exist, choose the first one
					max_value = dp[i-1][j].dh;
					D_beg = dp[i-1][j].D_beg;
					D_end = j;
					event = NEW_COPY;
				}
			}

			// 0 is excluded because a match/mismatch is mandatory
			if (event == START_REP) {
				// Open duplication
				for (int j = 1; j <= D_end - MIN_UNIT + 1; j++) {
					int tmp = (seq[i-1] == seq[j-1] ? TR_MAT_SCORE : TR_MIS_PEN) + OPEN_TR;
					dp[i][j].D_gate = max_value + tmp;
					dp[i][j].D_beg = j;
					dp[i][j].D_end = D_end;
					dp[i][j].pi = i-1;
					dp[i][j].pj = D_end;
					dp[i][j].event = START_REP;
				}
			} else {
				assert(D_beg > 0);
				for (int j = D_beg; j <= D_end - MIN_UNIT + 1; j++) {
					// No penalty for new copy
					int tmp = (seq[i-1] == seq[j-1] ? TR_MAT_SCORE : TR_MIS_PEN);
					dp[i][j].D_gate = max_value + tmp;
					dp[i][j].D_beg = D_beg; // Inherit from previous gate
					dp[i][j].D_end = D_end;
					dp[i][j].pi = i-1;
					dp[i][j].pj = D_end;
					dp[i][j].event = NEW_COPY;
				}
			}
		}

		// Repetition alignment starting from D gates
		for (int j = 1; j < i; j++) {
			int v_score = -INF, h_score = -INF, d_score = -INF;
			if (j >= dp[i-1][j].D_beg and j <= dp[i-1][j].D_end) {
				v_score = max(max(dp[i-1][j].D_gate, dp[i-1][j].dh) + TR_GAP_O, dp[i-1][j].de) + TR_GAP_E;
			}
			if (j >= dp[i][j-1].D_beg and j <= dp[i][j-1].D_end) {
				h_score = max(max(dp[i][j-1].D_gate, dp[i][j-1].dh) + TR_GAP_O, dp[i][j-1].df) + TR_GAP_E;
			}
			if (j >= dp[i-1][j-1].D_beg and j <= dp[i-1][j-1].D_end) {
				int tmp = (seq[i-1] == seq[j-1] ? TR_MAT_SCORE : TR_MIS_PEN);
				d_score = max(dp[i-1][j-1].D_gate, dp[i-1][j-1].dh) + tmp;
			}

			// Here, >= prefers another copy instead of a new copy; it can generate longer tandem repeats
			if (v_score >= dp[i][j].D_gate and v_score > dp[i][j].dh) {
				dp[i][j].dh = v_score;
				dp[i][j].D_beg = dp[i - 1][j].D_beg;
				dp[i][j].D_end = dp[i - 1][j].D_end;
				dp[i][j].pi = i-1;
				dp[i][j].pj = j;
				dp[i][j].event = WITHIN_REP;
			}
			if (h_score >= dp[i][j].D_gate and h_score > dp[i][j].dh) {
				dp[i][j].dh = h_score;
				dp[i][j].D_beg = dp[i][j - 1].D_beg;
				dp[i][j].D_end = dp[i][j - 1].D_end;
				dp[i][j].pi = i;
				dp[i][j].pj = j-1;
				dp[i][j].event = WITHIN_REP;
			}
			if (d_score >= dp[i][j].D_gate and d_score > dp[i][j].dh) {
				dp[i][j].dh = d_score;
				dp[i][j].D_beg = dp[i - 1][j - 1].D_beg;
				dp[i][j].D_end = dp[i - 1][j - 1].D_end;
				dp[i][j].pi = i-1;
				dp[i][j].pj = j-1;
				dp[i][j].event = WITHIN_REP;
			}
		}
		// The alignment above can't reach the diagonal

		// B transfer: close a repetition
		int max_value = -INF, max_j = -1;
		for (int j = 1; j < i; j++) {
			// Only return to diagonal if sub-matrix reaches the lower-right corner
			if (dp[i][j].dh > max_value and dp[i][j].D_end == j) {
				max_value = dp[i][j].dh;
				max_j = j;
			}
		}
		if (max_value + CLOSE_TR > dp[i][i].H) {
			dp[i][i].H = max_value + CLOSE_TR;
			dp[i][i].pi = i;
			dp[i][i].pj = max_j;
			dp[i][i].event = END_REP;
		}
	}

	if (not vis_fn.empty()) {
		ofstream out(vis_fn);
		assert(out.is_open());
		int ti = n, tj = n, te = -1;
		while (ti > 0 and tj > 0) {
			const Dp1Cell &c = dp[ti][tj];
			if (c.event != te) {
				out << ti << "\t" << tj << "\t" << c.event << endl;
			}
			ti = c.pi;
			tj = c.pj;
			te = c.event;
		}
		out << 0 << "\t" << 0 << "\t" << te << endl;
		out.close();
	}

	// Trace back the optimal path
	int ti = n, tj = n;
	vector<RepUnit> repetitions;
	while (ti > 0 and tj > 0) {
		Dp1Cell &t = dp[ti][tj];
		if (t.event == END_REP) {
			assert(ti == tj); // Only main diagonal closes repetitions
			if (DEBUG) fprintf(stderr, "Close: %d -> %d (Diagonal)\n", t.pj, ti);

			while (t.event != START_REP) {
				ti = t.pi;
				tj = t.pj;
				t = dp[ti][tj]; // Lower-right corner of the sub-matrix
				RepUnit u;
				u.score = t.dh;
				u.qe = ti + 1;
				u.te = tj + 1;
				while (t.event != NEW_COPY and t.event != START_REP) {
					if (ti == t.pi + 1 and tj == t.pj + 1) {
						if (seq[ti - 1] == seq[tj - 1]) u.match++;
						else u.mis++;
					} else {
						u.gap++;
					}
					ti = t.pi;
					tj = t.pj;
					t = dp[ti][tj];
				}
				// Pointer is now at the upper-left corner
				u.score -= dp[ti][tj].D_gate;
				u.qb = ti;
				u.tb = tj;
				if (seq[ti - 1] == seq[tj - 1]) {
					u.match++;
					u.score += o.tr_mat_score;
				} else {
					u.mis++;
					u.score += o.tr_mis_pen;
				}
				repetitions.push_back(u);

				if (t.event == START_REP) {
					if (DEBUG) fprintf(stderr, "Open: %d (Diagonal) -> %d, unit_size=%d\n", tj, t.pj, t.pj - tj + 1);
				} else if (t.event == NEW_COPY) {
					if (DEBUG) fprintf(stderr, "Copy: %d (i=%d) -> %d, unit_size=%d\n", tj, ti, t.pj, t.pj - tj + 1);
				}
			}
		}
		ti = t.pi;
		tj = t.pj;
	}

	reverse(repetitions.begin(), repetitions.end());
	fprintf(stderr, "Self-alignment found %ld duplications in %.3f CPU seconds\n", repetitions.size(), cputime() - ctime);
	if (DEBUG) {
		for (int i = 0; i < repetitions.size(); i++) {
			const RepUnit &u = repetitions[i];
			fprintf(stdout, "[%d] Tandem repeat between [%d,%d) and [%d,%d), score=%d (matches=%d, mismatches=%d, gaps=%d)\n",
			        i+1, u.qb, u.qe, u.tb, u.te, u.score, u.match, u.mis, u.gap);
		}
	}
	return repetitions;
}

struct Dp2Cell {
	int E, F, B1, B2, H;
	int pi, pj;
	Dp2Cell() {
		E = F = H = B1 = B2 = -INF;
		pi = pj = -1;
	}
};

void align_with_dups(const ZigOptions &opt, const char *fn1, const char *fn2)
{
	string seq1 = input_fasta_seq(fn1);
	string seq2 = input_fasta_seq(fn2);
	const int n = seq1.length();
	const int m = seq2.length();
	const char *a = seq1.data();
	const char *b = seq2.data();
	string vis_self_fn1, vis_self_fn2;
	if (opt.vis_prefix) {
		vis_self_fn1 = string(opt.vis_prefix) + "_s1.txt";
		vis_self_fn2 = string(opt.vis_prefix) + "_s2.txt";
	}
	vector<RepUnit> rep_a = self_alignment(opt, n, a, vis_self_fn1);
	vector<RepUnit> rep_b = self_alignment(opt, m, b, vis_self_fn2);
	// Set end positions (exclusive) as break points
	vector<bool> bp_a(n + 1, false);
	vector<bool> bp_b(m + 1, false);
	for (const RepUnit &u : rep_a) {
		bp_a[u.qe-1] = true; // 1-based; open right interval
		bp_a[u.te-1] = true; // It allows the deletion of the second unit
	}
	for (const RepUnit &u : rep_b) {
		bp_b[u.qe-1] = true;
		bp_b[u.te-1] = true;
	}
	if (DEBUG) {
		for (int i = 0; i < n; i++) {
			if (bp_a[i+1]) fprintf(stderr, "|");
			fprintf(stderr, "%c", a[i]);
		}
		fprintf(stderr, "\n");
		for (int i = 0; i < m; i++) {
			if (bp_b[i+1]) fprintf(stderr, "|");
			fprintf(stderr, "%c", b[i]);
		}
		fprintf(stderr, "\n");
	}
	int sum1 = accumulate(bp_a.begin(), bp_a.end(), 0);
	int sum2 = accumulate(bp_b.begin(), bp_b.end(), 0);
	fprintf(stderr, "Identified %d and %d break points in two sequences\n", sum1, sum2);

	const int MAT_SCORE = opt.mat_score;
	const int MIS_PEN = opt.mis_pen;
	const int GAP_O = opt.gap_o;
	const int GAP_E = opt.gap_e;
	const int DEL_DUP_O = opt.del_dup_o;
	const int DEL_DUP_E = opt.del_dup_e;
	vector<vector<Dp2Cell>> dp;
	dp.resize(n + 1);
	for (int i = 0; i <= n; i++) {
		dp[i].resize(m + 1);
	}
	dp[0][0].H = 0;
	for (int j = 1; j <= m; j++) {
		dp[0][j].F = dp[0][j].H = GAP_O + j * GAP_E;
	}
	int last_row = -1;
	for (int i = 1; i <= n; i++) {
		dp[i][0].E = dp[i][0].H = GAP_O + i * GAP_E;
		int last_col = -1;
		for (int j = 1; j <= m; j++) {
			dp[i][j].E = max(dp[i-1][j].H + GAP_O, dp[i-1][j].E) + GAP_E;
			dp[i][j].F = max(dp[i][j-1].H + GAP_O, dp[i][j-1].F) + GAP_E;
			int M = dp[i-1][j-1].H + (a[i-1] == b[j-1] ? MAT_SCORE : MIS_PEN);
			if (dp[i][j].E > dp[i][j].H) {
				dp[i][j].H = dp[i][j].E;
				dp[i][j].pi = i-1;
				dp[i][j].pj = j;
			}
			if (dp[i][j].F > dp[i][j].H) {
				dp[i][j].H = dp[i][j].F;
				dp[i][j].pi = i;
				dp[i][j].pj = j-1;
			}
			if (M > dp[i][j].H) {
				dp[i][j].H = M;
				dp[i][j].pi = i-1;
				dp[i][j].pj = j-1;
			}
			// New path between two break points
			if (bp_b[j] and last_col != -1) {
				dp[i][j].B2 = max(dp[i][last_col].H + DEL_DUP_O, dp[i][last_col].B2) + DEL_DUP_E;
				if (dp[i][j].B2 > dp[i][j].H) {
					dp[i][j].H = dp[i][j].B2;
					dp[i][j].pi = i;
					dp[i][j].pj = last_col;
				}
			}
			if (bp_b[j]) last_col = j;

			if (bp_a[i] and last_row != -1) {
				dp[i][j].B1 = max(dp[last_row][j].H + DEL_DUP_O, dp[last_row][j].B1) + DEL_DUP_E;
				if (dp[i][j].B1 > dp[i][j].H) {
					dp[i][j].H = dp[i][j].B1;
					dp[i][j].pi = last_row;
					dp[i][j].pj = j;
				}
			}
		}
		if (bp_a[i]) last_row = i;
	}

	int ti = n, tj = m;
	int del_n = 0, ins_n = 0, mat_n = 0, mis_n = 0, dup_n = 0;
	string ext_a, ext_b;
	while (ti > 0 and tj > 0) {
		const Dp2Cell &t = dp[ti][tj];
		if (t.pi == ti - 1 and t.pj == tj) {
			del_n++;
			ext_a += a[ti-1];
			ext_b += '-';
		} else if (t.pi == ti and t.pj == tj - 1) {
			ins_n++;
			ext_a += '-';
			ext_b += b[tj-1];
		} else if (t.pi == ti - 1 and t.pj == tj - 1) {
			if (a[ti-1] == b[tj-1]) mat_n++;
			else mis_n++;
			ext_a += a[ti-1];
			ext_b += b[tj-1];
		} else {
//			fprintf(stderr, "%d %d -> %d %d\n", ti, tj, t.pi, t.pj);
			if (ti == t.pi) {
				for (int i = tj-1; i >= t.pj+1; i--) {
					ext_a += '+';
					ext_b += b[i-1];
				}
			} else {
				for (int i = ti-1; i >= t.pi+1; i--) {
					ext_a += a[i-1];
					ext_b += '+';
				}
			}
			dup_n++;
		}
		ti = t.pi;
		tj = t.pj;
	}
	if (ti > 0) del_n += ti;
	if (tj > 0) ins_n += tj;
	fprintf(stderr, "%d deletions, %d insertions, %d mismatches and %d duplications\n", del_n, ins_n, mis_n, dup_n);
	reverse(ext_a.begin(), ext_a.end());
	reverse(ext_b.begin(), ext_b.end());
	fprintf(stderr, "%s\n", ext_a.data());
	fprintf(stderr, "%s\n", ext_b.data());

	if (opt.vis_prefix) {
		string pair_vis_fn = string(opt.vis_prefix) + "_p.txt";
		ofstream out(pair_vis_fn);
		assert(out.is_open());

		out.close();
	}
}

int usage(const ZigOptions &o) {
	fprintf(stderr, "Usage: zigalign [options] seq1.fa seq2.fa\n");
	fprintf(stderr, "  Regular Scoring options:\n");
	fprintf(stderr, "    -A [INT]  match score [%d]\n", o.mat_score);
	fprintf(stderr, "    -B [INT]  mismatch penalty [%d]\n", o.mis_pen);
	fprintf(stderr, "    -O [INT]  open gap(indel) penalty [%d]\n", o.gap_o);
	fprintf(stderr, "    -E [INT]  extend gap penalty [%d]\n", o.gap_e);
	fprintf(stderr, "    -J [INT]  open delete duplication penalty [%d]\n", o.del_dup_o);
	fprintf(stderr, "    -j [INT]  extend delete duplication penalty [%d]\n", o.del_dup_e);
	fprintf(stderr, "  Scoring options for self-alignment:\n");
	fprintf(stderr, "    -u [INT]  minimum repeat unit size [%d]\n", o.min_unit_size);
	fprintf(stderr, "    -d [INT]  open tandem repeat penalty [%d]\n", o.open_tr_pen);
	fprintf(stderr, "    -p [INT]  close tandem repeat penalty [%d]\n", o.close_tr_pen);
	fprintf(stderr, "    -a [INT]  match score [%d]\n", o.tr_mat_score);
	fprintf(stderr, "    -b [INT]  mismatch penalty [%d]\n", o.tr_mis_pen);
	fprintf(stderr, "    -o [INT]  open gap(indel) penalty [%d]\n", o.tr_gap_o);
	fprintf(stderr, "    -e [INT]  extend gap penalty [%d]\n", o.tr_gap_e);
	fprintf(stderr, "Note: self-alignment scoring matrix must reward more and/or \n"
	                "  penalize less than regular matrix to discover tandem repeats.\n");
	return 1;
}

int main(int argc, char *argv[]) {
	double ctime = cputime(), rtime = realtime();
	ZigOptions opt;
	if (argc == 1) return usage(opt);
	int c;
	while ((c = getopt(argc, argv, "A:B:O:E:J:j:u:d:p:a:b:o:e:v:")) >= 0) {
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
			case 'J':
				opt.del_dup_o = -abs(str2int(optarg));
				break;
			case 'j':
				opt.del_dup_e = -abs(str2int(optarg));
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
				opt.tr_mat_score = abs(str2int(optarg));
				break;
			case 'b':
				opt.tr_mis_pen = -abs(str2int(optarg));
				break;
			case 'o':
				opt.tr_gap_o = -abs(str2int(optarg));
				break;
			case 'e':
				opt.tr_gap_e = -abs(str2int(optarg));
				break;
			case 'v':
				opt.vis_prefix = optarg;
				break;
			default:
				fprintf(stderr, "Unrecognized option `%c`\n", c);
				return 1;
		}
	}

	if (argc - optind == 2) {
		align_with_dups(opt, argv[optind], argv[optind+1]);
	} else {
		fprintf(stderr, "Two FASTA files are required\n");
		return 1;
	}
	fprintf(stderr, "Program finishes in %.3f CPU seconds, %.3f real seconds\n", cputime()-ctime, realtime()-rtime);
	return 0;
}