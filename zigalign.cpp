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
#include <ratio>

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
	// K-band size
	int band_width;
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
		band_width = 300;
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

struct Interval {
	int l, r; // [l, r)
};

struct TandemRepeat {
	Interval temp; // Templates
	vector<Interval> units;
};

// Stage 1: identify breakpoints of tandem repeats using self alignment
static vector<int> self_alignment(const ZigOptions &o, int n, const char *seq, const string &vis_fn)
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
	dp[0][0].D_beg = 1;
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
		dp[i][i].D_beg = dp[i-1][i-1].D_beg;

		if (i > MIN_UNIT) {
			// D gate comes from the last row
			int max_value = dp[i-1][i-1].H; // Diagonal must be the maximum in regular matrix
			int D_end = i - 1;
			int D_beg = dp[i-1][i-1].D_beg; // To prevent illegal path
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
				for (int j = D_beg; j <= D_end - MIN_UNIT + 1; j++) {
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
					dp[i][j].D_beg = j; // Does it make more sense?
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
				dp[i][j].de = v_score;
			}
			if (j >= dp[i][j-1].D_beg and j <= dp[i][j-1].D_end) {
				h_score = max(max(dp[i][j-1].D_gate, dp[i][j-1].dh) + TR_GAP_O, dp[i][j-1].df) + TR_GAP_E;
				dp[i][j].df = h_score;
			}
			if (j >= dp[i-1][j-1].D_beg and j <= dp[i-1][j-1].D_end) {
				int tmp = (seq[i-1] == seq[j-1] ? TR_MAT_SCORE : TR_MIS_PEN);
				d_score = max(dp[i-1][j-1].D_gate, dp[i-1][j-1].dh) + tmp;
			}

			// Here, >= prefers continue an existing copy, so it can generate longer tandem repeats
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
		int max_value = -INF, max_j = -1, limit = -1;
		for (int j = 1; j < i; j++) {
			// Only return to diagonal if sub-matrix reaches the lower-right corner
			if (dp[i][j].dh > max_value and dp[i][j].D_end == j) {
				max_value = dp[i][j].dh;
				max_j = j;
				limit = dp[i][j].D_beg;
			}
		}
		if (max_value + CLOSE_TR > dp[i][i].H) {
			dp[i][i].H = max_value + CLOSE_TR;
			dp[i][i].pi = i;
			dp[i][i].pj = max_j;
			dp[i][i].D_beg = limit;
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

	// Trace back the optimal path to find tandem repeats
	int ti = n, tj = n;
	vector<TandemRepeat> tr_array;
	while (ti > 0 and tj > 0) {
		Dp1Cell &t = dp[ti][tj];
		if (t.event == END_REP) {
			assert(ti == tj); // Only main diagonal closes repetitions
			TandemRepeat group;
			while (t.event != START_REP) {
				ti = t.pi;
				tj = t.pj;
				t = dp[ti][tj]; // Lower-right corner of the sub-matrix
				Interval intv;
				intv.r = ti + 1;
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
				intv.l = ti;
				group.units.push_back(intv);
				if (seq[ti - 1] == seq[tj - 1]) {
					u.match++;
					u.score += o.tr_mat_score;
				} else {
					u.mis++;
					u.score += o.tr_mis_pen;
				}
			}
			group.temp.l = tj;
			group.temp.r = ti;
			reverse(group.units.begin(), group.units.end());
			tr_array.push_back(group);
		}
		ti = t.pi;
		tj = t.pj;
	}
	reverse(tr_array.begin(), tr_array.end());
	if (DEBUG) {
		fprintf(stdout, "Tandem repeat groups\n");
		for (const TandemRepeat &group: tr_array) {
			fprintf(stdout, "Template: [%d, %d), Units:", group.temp.l, group.temp.r);
			for (const Interval &intv: group.units) {
				fprintf(stdout, " [%d, %d)", intv.l, intv.r);
			}
			fprintf(stdout, "\n");
		}
	}

	// Merge tandem groups
	double MIN_OVERLAP_COEF = 0.9;
	int new_size = 0;
	for (int i = 0; i < tr_array.size(); ) {
		TandemRepeat &tr = tr_array[i];
		int range_r = tr.units.back().r;
		int next = tr_array.size();
		for (int j = i + 1; j < tr_array.size(); j++) {
			int l = tr_array[j].temp.l;
			int r = tr_array[j].temp.r;
			int len = r - l;
			int overlap = max(range_r - l, 0);
			if (overlap > MIN_OVERLAP_COEF * len) {
				// Merge into previous group
				for (const Interval &intv: tr_array[j].units) {
					tr.units.push_back(intv);
				}
				range_r = tr.units.back().r;
			} else {
				next = j;
				break;
			}
		}
		tr_array[new_size++] = tr;
		i = next;
	}

	// FIXME: This is incorrect because it mixes different groups of tandem repeats
	int rep_end = -1;
	vector<int> ret;
	for (const TandemRepeat &group: tr_array) {
		if (group.temp.l-1 > rep_end) ret.push_back(group.temp.l-1);
		if (group.temp.r-1 > rep_end) ret.push_back(group.temp.r-1);
		for (int i = 0; i < group.units.size(); i++) {
			if (group.units[i].r-1 > rep_end) {
				ret.push_back(group.units[i].r-1);
			}
		}
		if (not ret.empty()) rep_end = ret.back();
	}

	if (ret.empty()) return ret;

	// Breakpoints must be naturally sorted
	for (int i = 1; i < ret.size(); i++) {
		assert(ret[i] > ret[i-1]);
	}

	// Breakpoints are 1-based

	if (1) {
		fprintf(stdout, "Breakpoints:");
		for (int i : ret) {
			fprintf(stdout, " %d", i);
		}
		fprintf(stdout, "\n");
	}

	fprintf(stderr, "Self-alignment identified %ld breakpoints in %.3f CPU seconds\n", ret.size(), cputime() - ctime);
	return ret;
}

void paf_format(const string &q_name, const string &que, const string &t_name, string tar, const vector<int> &cv)
{
	const int COP_M = 0;
	const int COP_I = 1;
	const int COP_D = 2;
	const int DUP_I = 3;
	const int DUP_D = 4;
	const string OP_CHAR = "MIDID";

	int q_len = que.length(), t_len = tar.length();
	char strand = '+';
	int matches_n = 0, mismatches_n = 0, extended_length = 0;
	int mapq = 60;
	int edit_distance = 0;
	int ti = 0, qi = 0;
	string cigar;
	for (int x : cv) {
		int op_type = x & 15, op_cnt = x >> 4;
		if (op_type == COP_M) {
			for (int i = 0; i < op_cnt; i++) {
				if (que[qi + i] == tar[ti + i]) {
					matches_n++;
				} else {
					mismatches_n++;
					edit_distance++;
				}
			}
			ti += op_cnt;
			qi += op_cnt;
		} else if (op_type == COP_I) {
			qi += op_cnt;
			edit_distance += op_cnt;
		} else if (op_type == DUP_I) {
			qi += op_cnt;
			edit_distance += 1; // Tandem duplication has an edit distance of 1
		} else if (op_type == COP_D) {
			ti += op_cnt;
			edit_distance += op_cnt;
		} else {
			ti += op_cnt;
			edit_distance += 1; // Tandem deletion has an edit distance of 1
		}
		extended_length += op_cnt;
		// TODO: identify matches between units
		if (op_type == DUP_D or op_type == DUP_I) cigar += 'U';
		cigar += to_string(op_cnt);
		cigar += OP_CHAR[op_type];
	}
	assert(ti == t_len and qi == q_len);

	// CIGAR validation
	if (DEBUG) {
		string modified_tar;
		qi = 0;
		for (int x : cv) {
			int op_type = x & 15, op_cnt = x >> 4;
			if (op_type == COP_M) {
				for (int i = 0; i < op_cnt; i++) {
					modified_tar += que[qi++];
				}
			} else if (op_type == COP_I or op_type == DUP_I) {
				for (int i = 0; i < op_cnt; i++) {
					modified_tar += que[qi++];
				}
			}
		}
		assert(modified_tar == que);
	}
	fprintf(stdout, "%s\t%d\t%d\t%d\t%c\t", q_name.c_str(), q_len, 0, q_len, strand);
	fprintf(stdout, "%s\t%d\t%d\t%d\t", t_name.c_str(), t_len, 0, t_len);
	fprintf(stdout, "%d\t%d\t%d\t", matches_n+mismatches_n, extended_length, mapq);
	fprintf(stdout, "NM:i:%d\t", edit_distance);
	fprintf(stdout, "cg:Z:%s\n", cigar.c_str());
}

struct Dp2Cell {
	int E, F, B1, B2, H;
	int pi, pj;
	Dp2Cell() {
		E = F = H = B1 = B2 = -INF;
		pi = pj = -1;
	}
};

void align_with_dups(const ZigOptions &opt,
	int t_len, const char *t, const vector<int> &t_bp,
	int q_len, const char *q, const vector<int> &q_bp)
{
	// Breakpoints from self-alignment are not optimal in pairwise alignment
	const int BP_HALF_KMER = 5; // 5 bp before and after the breakpoint, i.e., 11-mer
	vector<Interval> t_bp_intv;
	for (int i = 0; i < t_bp.size(); i++) {
		// If breakpoint kmers overlap with other kmers, then use the breakpoint itself.
		bool overlapped = false;
		if (i > 0 and t_bp[i] - BP_HALF_KMER <= t_bp[i-1] + BP_HALF_KMER) overlapped = true;
		if (i+1 < t_bp.size() and t_bp[i] + BP_HALF_KMER >= t_bp[i+1] - BP_HALF_KMER) overlapped = true;
		Interval intv;
		if (not overlapped) {
			intv.l = max(0, t_bp[i] - BP_HALF_KMER);
			intv.r = min(t_len, t_bp[i] + BP_HALF_KMER);
		} else {
			intv.l = t_bp[i];
			intv.r = t_bp[i];
		}
		t_bp_intv.push_back(intv);
	}
	vector<Interval> q_bp_intv;
	for (int i = 0; i < q_bp.size(); i++) {
		bool overlapped = false;
		if (i > 0 and q_bp[i] - BP_HALF_KMER <= q_bp[i-1] + BP_HALF_KMER) overlapped = true;
		if (i+1 < q_bp.size() and q_bp[i] + BP_HALF_KMER >= q_bp[i+1] - BP_HALF_KMER) overlapped = true;
		Interval intv;
		if (not overlapped) {
			intv.l = max(0, q_bp[i] - BP_HALF_KMER);
			intv.r = min(q_len, q_bp[i] + BP_HALF_KMER);
		} else {
			intv.l = q_bp[i];
			intv.r = q_bp[i];
		}
		q_bp_intv.push_back(intv);
	}

	// Pairwise alignment
	const int MAT_SCORE = opt.mat_score;
	const int MIS_PEN = opt.mis_pen;
	const int GAP_O = opt.gap_o;
	const int GAP_E = opt.gap_e;
	const int DEL_DUP_O = opt.del_dup_o;
	const int DEL_DUP_E = opt.del_dup_e;
	vector<vector<Dp2Cell>> dp;
	dp.resize(t_len + 1);
	for (int i = 0; i <= t_len; i++) {
		dp[i].resize(q_len + 1);
	}
	dp[0][0].H = 0;
	for (int j = 1; j <= q_len; j++) {
		dp[0][j].F = dp[0][j].H = GAP_O + j * GAP_E;
	}
	int t_pointer = 0;
	for (int i = 1; i <= t_len; i++) {
		dp[i][0].E = dp[i][0].H = GAP_O + i * GAP_E;
		while (t_pointer < t_bp_intv.size() and t_bp_intv[t_pointer].r < i) t_pointer++;
		bool in_t_bp = false;
		if (t_pointer < t_bp_intv.size() and t_bp_intv[t_pointer].l <= i) in_t_bp = true;
		if (in_t_bp) {
			assert(t_pointer < t_bp_intv.size() and t_bp_intv[t_pointer].l <= i and t_bp_intv[t_pointer].r >= i);
			// cout << i << " is in the bp interval " << t_bp_intv[t_pointer].l << " " << t_bp_intv[t_pointer].r << endl;
			// if (i == t_bp_intv[t_pointer].l) {
			// 	cout << i << " enter the interval " << t_bp_intv[t_pointer].l << " " << t_bp_intv[t_pointer].r << endl;
			// }
			// if (i == t_bp_intv[t_pointer].r) {
			// 	cout << i << " leave the interval " << t_bp_intv[t_pointer].l << " " << t_bp_intv[t_pointer].r << endl;
			// }
		}

		int q_pointer = 0;
		for (int j = 1; j <= q_len; j++) {
			while (q_pointer < q_bp_intv.size() and q_bp_intv[q_pointer].r < j) q_pointer++;
			bool in_q_bp = false;
			if (q_pointer < q_bp_intv.size() and q_bp_intv[q_pointer].l <= j) in_q_bp = true;
			if (in_q_bp) {
				assert(q_pointer < q_bp_intv.size() and q_bp_intv[q_pointer].l <= j and q_bp_intv[q_pointer].r >= j);
			}

			if (in_t_bp and in_q_bp) {
				if (q_pointer > 0) {
					// Jump from the last breakpoint kmer
					int l1 = q_bp_intv[q_pointer-1].l, r1 = q_bp_intv[q_pointer-1].r;
					int max_score = -INF, max_id = -1;
					// Redundant caculation
					for (int k = l1; k <= r1; k++) {
						int tmp = max(dp[i][k].H + DEL_DUP_O, dp[i][k].B2) + DEL_DUP_E;
						if (tmp > max_score) {
							max_score = tmp;
							max_id = k;
						}
					}
					dp[i][j].B2 = max_score;
					if (dp[i][j].B2 > dp[i][j].H) {
						dp[i][j].H = dp[i][j].B2;
						dp[i][j].pi = i;
						dp[i][j].pj = max_id;
					}
				}
				if (t_pointer > 0) {
					int l1 = t_bp_intv[t_pointer-1].l, r1 = t_bp_intv[t_pointer-1].r;
					int max_score = -INF, max_id = -1;
					for (int k = l1; k <= r1; k++) {
						int tmp = max(dp[k][j].H + DEL_DUP_O, dp[k][j].B1) + DEL_DUP_E;
						if (tmp > max_score) {
							max_score = tmp;
							max_id = k;
						}
					}
					dp[i][j].B1 = max_score;
					if (dp[i][j].B1 > dp[i][j].H) {
						dp[i][j].H = dp[i][j].B1;
						dp[i][j].pi = max_id;
						dp[i][j].pj = j;
					}
				}
			}

			dp[i][j].E = max(dp[i-1][j].H + GAP_O, dp[i-1][j].E) + GAP_E;
			dp[i][j].F = max(dp[i][j-1].H + GAP_O, dp[i][j-1].F) + GAP_E;
			int M = dp[i-1][j-1].H + (t[i-1] == q[j-1] ? MAT_SCORE : MIS_PEN);
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
		}
	}

	cout << "Truth score: " << dp[t_len][q_len].H << endl;

	// CIGAR generation
	int ti = t_len, tj = q_len;
	int del_n = 0, ins_n = 0, mat_n = 0, mis_n = 0, dup_n = 0;
	string ext_t, ext_q;
	vector<int> cv;
	// All operations occur on target sequence _t_
	const int COP_M = 0;
	const int COP_I = 1;
	const int COP_D = 2;
	const int DUP_I = 3;
	const int DUP_D = 4;
	const string OP_CHAR = "MIDID";
	int op_type = -1, last_op = -1, op_cnt = 0;
	while (ti > 0 and tj > 0) {
		const Dp2Cell &p = dp[ti][tj];
		if (p.pi == ti - 1 and p.pj == tj) { // Deletion from _t_
			op_type = COP_D;
			del_n++;
			ext_t += t[ti-1];
			ext_q += '-';
		} else if (p.pi == ti and p.pj == tj - 1) { // Insertion into _t_
			op_type = COP_I;
			ins_n++;
			ext_t += '-';
			ext_q += q[tj-1];
		} else if (p.pi == ti - 1 and p.pj == tj - 1) { // Match/Mismatch
			op_type = COP_M;
			if (t[ti-1] == q[tj-1]) mat_n++;
			else mis_n++;
			ext_t += t[ti-1];
			ext_q += q[tj-1];
		} else {
			if (last_op != -1) cv.push_back(op_cnt << 4 | last_op);
//			fprintf(stderr, "%d %d -> %d %d\n", ti, tj, t.pi, t.pj);
			if (ti == p.pi) { // Duplication insertion into _t_
				op_type = DUP_I;
				op_cnt = tj - p.pj;
				for (int j = tj; j > p.pj; j--) {
					ext_t += '+';
					ext_q += q[j-1];
				}
			} else { // Duplication deletion from _t_
				op_type = DUP_D;
				op_cnt = ti - p.pi;
				for (int i = ti; i > p.pi; i--) {
					ext_t += t[i-1];
					ext_q += '+';
				}
			}

			// TODO: is it necessary to align template unit to copied units?
			// If so, which unit is the best template?
			cv.push_back(op_cnt << 4 | op_type); // Do not merge duplication indels
			dup_n++;
			op_type = last_op = -1;
			op_cnt = -1;
		}
		// Merge match/mismatch/indels
		if (last_op != -1 and op_type != last_op) {
			cv.push_back(op_cnt << 4 | last_op);
			op_cnt = 0;
		}
		last_op = op_type;
		op_cnt++;
		ti = p.pi;
		tj = p.pj;
	}
	if (last_op != -1) cv.push_back(op_cnt << 4 | last_op);
	if (ti > 0) {
		cv.push_back(ti << 4 | COP_D);
	}
	if (tj > 0) {
		cv.push_back(tj << 4 | COP_I);
	}
	reverse(cv.begin(), cv.end());
	reverse(ext_t.begin(), ext_t.end());
	reverse(ext_q.begin(), ext_q.end());
	fprintf(stderr, "%d deletions, %d insertions, %d matches, %d mismatches and %d duplication indels\n", del_n, ins_n, mat_n, mis_n, dup_n);

	// Sanity check
	if (1) {
		if (DEBUG) {
			fprintf(stderr, "%s\n", ext_t.data());
			fprintf(stderr, "%s\n", ext_q.data());
		}
		string non_t, non_q;
		for (char c: ext_t) {
			if (c != '-' and c != '+') {
				non_t += c;
			}
		}
		for (char c: ext_q) {
			if (c != '-' and c != '+') {
				non_q += c;
			}
		}
		assert(non_t.length() == t_len);
		assert(non_t == string(t));
		assert(non_q.length() == q_len);
		assert(non_q == string(q));
	}

	if (1) {
		string pair_vis_fn = "test/dev_p.txt";
		ofstream out(pair_vis_fn);
		assert(out.is_open());
		// Meta information and breakpoints
		out << "Seq1 length: " << t_len << ", Seq2 length: " << q_len << endl;
		out << "Breakpoints1:" << endl;
		for (int i: t_bp) out << i << " ";
		out << endl;
		out << "Breakpoints2:" << endl;
		for (int i: q_bp) out << i << " ";
		out << endl;

		ti = t_len; tj = q_len;
		int last_type = -1;
		const int TYPE_DEL = 0;
		const int TYPE_INS = 1;
		const int TYPE_MAT = 2;
		const int TYPE_DUP = 3;
		while (ti > 0 and tj > 0) {
			const Dp2Cell &p = dp[ti][tj];
			int type;
			if (p.pi == ti - 1 and p.pj == tj) {
				type = TYPE_DEL;
			} else if (p.pi == ti and p.pj == tj - 1) {
				type = TYPE_INS;
			} else if (p.pi == ti - 1 and p.pj == tj - 1) {
				type = TYPE_MAT;
			} else {
				type = TYPE_DUP;
			}
			// Do not merge tandem duplications
			if (type != last_type or type == TYPE_DUP) {
				out << ti << "\t" << tj << "\t" << type << endl;
			}
			last_type = type;
			ti = p.pi;
			tj = p.pj;
		}
		out << ti << "\t" << tj << "\t" << last_type << endl;
		out.close();
	}
}

struct Dp3Cell {
	int ci, cj; // Coordinates in DP matrix
	int E, F, B1, B2, H;
	int pi, pj; // Traceback in stored matrix
	Dp3Cell() {
		E = F = H = B1 = B2 = -INF;
		pi = pj = -1;
	}
};

vector<Interval> merge_intervals(const vector<Interval> &v) {
	// Input intervals must be l-sorted
	assert(v.size() > 0);
	for (int i = 1; i < v.size(); i++) {
		assert(v[i].l >= v[i-1].l);
	}
	vector<Interval> ret;
	Interval a = v[0];
	for (int i = 1; i < v.size(); i++) {
		const Interval &t = v[i];
		if (t.l > a.r + 1) {
			ret.push_back(a);
			a = t;
		} else {
			a.r = max(t.r, a.r); // Merge overlapping intervals
		}
	}
	ret.push_back(a);
	return ret;
}

int MAX_BAND_WIDTH = 100;

void align_with_dups2(const ZigOptions &opt,
	int t_len, const char *t, const vector<int> &t_bp,
	int q_len, const char *q, const vector<int> &q_bp)
{
	// Breakpoints from self-alignment are not optimal in pairwise alignment
	const int BP_HALF_KMER = 5; // 5 bp before and after the breakpoint, i.e., 11-mer
	vector<Interval> t_bp_intv;
	for (int i = 0; i < t_bp.size(); i++) {
		// If breakpoint kmers overlap with other kmers, then use the breakpoint itself.
		bool overlapped = false;
		if (i > 0 and t_bp[i] - BP_HALF_KMER <= t_bp[i-1] + BP_HALF_KMER) overlapped = true;
		if (i+1 < t_bp.size() and t_bp[i] + BP_HALF_KMER >= t_bp[i+1] - BP_HALF_KMER) overlapped = true;
		Interval intv;
		if (not overlapped) {
			intv.l = max(0, t_bp[i] - BP_HALF_KMER);
			intv.r = min(t_len, t_bp[i] + BP_HALF_KMER);
		} else {
			intv.l = t_bp[i];
			intv.r = t_bp[i];
		}
		t_bp_intv.push_back(intv);
	}
	vector<Interval> q_bp_intv;
	for (int i = 0; i < q_bp.size(); i++) {
		bool overlapped = false;
		if (i > 0 and q_bp[i] - BP_HALF_KMER <= q_bp[i-1] + BP_HALF_KMER) overlapped = true;
		if (i+1 < q_bp.size() and q_bp[i] + BP_HALF_KMER >= q_bp[i+1] - BP_HALF_KMER) overlapped = true;
		Interval intv;
		if (not overlapped) {
			intv.l = max(0, q_bp[i] - BP_HALF_KMER);
			intv.r = min(q_len, q_bp[i] + BP_HALF_KMER);
		} else {
			intv.l = q_bp[i];
			intv.r = q_bp[i];
		}
		q_bp_intv.push_back(intv);
	}

	// Pairwise alignment
	const int MAT_SCORE = opt.mat_score;
	const int MIS_PEN = opt.mis_pen;
	const int GAP_O = opt.gap_o;
	const int GAP_E = opt.gap_e;
	const int DEL_DUP_O = opt.del_dup_o;
	const int DEL_DUP_E = opt.del_dup_e;

	int bw = MAX_BAND_WIDTH; // Bandwidth
	vector<vector<Dp3Cell>> dp;
	dp.resize(t_len + 1);

	// The first row
	// In each row, I obtain the cells that will be calculated
	vector<Interval> calc_band;
	Interval cb;
	cb.l = -bw; // For movement of the band
	cb.r = bw;
	calc_band.push_back(cb);
//	if (t_bp_intv.size() > 0 and t_bp_intv[0].l <= 0 and t_bp_intv[0].r >= 0) {
//		for (const Interval &v: q_bp_intv) {
//			cb.l = v.l - bw;
//			cb.r = v.r + bw;
//			calc_band.push_back(cb);
//		}
//	}
	calc_band = merge_intervals(calc_band);
	for (int i = 0; i < calc_band.size(); i++) {
		const Interval &v = calc_band[i];
		// printf("[%d, %d]\n", v.l, v.r);
		int l = max(0, v.l);
		int r = min(q_len, v.r);
		for (int j = l; j <= r; j++) {
			Dp3Cell c;
			c.ci = 0;
			c.cj = j;
			dp[0].push_back(c);
		}
	}
	assert(dp[0][0].ci == 0 and dp[0][0].cj == 0);
	dp[0][0].H = 0; // Left-top corner of the partial matrix
	dp[0][0].E = dp[0][0].F = GAP_O; // For convenient calculation
	for (int j = 1; j < dp[0].size(); j++) {
		if (dp[0][j].cj == dp[0][j-1].cj + 1 and dp[0][j-1].F != -INF) {
			dp[0][j].F = dp[0][j].H = dp[0][j-1].F + GAP_E;
		} // else the previous cell is not stored
		// FIXME: I should process breakpoints in the first row
	}

	int t_pointer = 0;
	for (int i = 1; i <= t_len; i++) {
		while (t_pointer < t_bp_intv.size() and t_bp_intv[t_pointer].r < i) t_pointer++;
		bool in_t_bp = (t_pointer < t_bp_intv.size() and t_bp_intv[t_pointer].l <= i);

		// Shift the band
		for (Interval &v: calc_band) {
			v.l++;
			v.r++;
		}
		// Add bands of breakpoints
		if (in_t_bp) {
			for (const Interval &bp: q_bp_intv) {
				Interval v{bp.l - bw, bp.r + bw};
				calc_band.push_back(v);
			}
			sort(calc_band.begin(), calc_band.end(), [](const Interval &a, const Interval &b) -> bool { return a.l <= b.l; });
			calc_band = merge_intervals(calc_band);
		}
		for (Interval &v: calc_band) {
			int l = max(0, v.l);
			int r = min(q_len, v.r);
			for (int j = l; j <= r; j++) {
				Dp3Cell c;
				c.ci = i;
				c.cj = j;
				dp[i].push_back(c);
			}
		}

		int mon_j = 0;
		int q_pointer = 0;
		for (int j = 0; j < dp[i].size(); j++) {
			int curr_j = dp[i][j].cj;
			// Diagonal transfer
			while (mon_j < dp[i-1].size() and dp[i-1][mon_j].cj < curr_j-1) mon_j++;
			if (mon_j < dp[i-1].size() and dp[i-1][mon_j].cj == curr_j-1) {
				int M = dp[i-1][mon_j].H + (t[i-1] == q[curr_j-1] ? MAT_SCORE : MIS_PEN);
				if (M > dp[i][j].H) {
					dp[i][j].H = M;
					dp[i][j].pi = i-1;
					dp[i][j].pj = mon_j;
				}
			}

			// Vertical transfer
			while (mon_j < dp[i-1].size() and dp[i-1][mon_j].cj < curr_j) mon_j++;
			if (mon_j < dp[i-1].size() and dp[i-1][mon_j].cj == curr_j) {
				dp[i][j].E = max(dp[i-1][mon_j].H + GAP_O, dp[i-1][mon_j].E) + GAP_E;
				if (dp[i][j].E > dp[i][j].H) {
					dp[i][j].H = dp[i][j].E;
					dp[i][j].pi = i-1;
					dp[i][j].pj = mon_j;
				}
			}

			// Horizontal transfer
			if (j > 0 and curr_j == dp[i][j-1].cj + 1) {
				dp[i][j].F = max(dp[i][j-1].H + GAP_O, dp[i][j-1].F) + GAP_E;
				if (dp[i][j].F > dp[i][j].H) {
					dp[i][j].H = dp[i][j].F;
					dp[i][j].pi = i;
					dp[i][j].pj = j-1;
				}
			}

			while (q_pointer < q_bp_intv.size() and q_bp_intv[q_pointer].r < curr_j) q_pointer++;
			bool in_q_bp = (q_pointer < q_bp_intv.size() and q_bp_intv[q_pointer].l <= curr_j);
			if (in_t_bp and in_q_bp) {
				// TODO: remove the redundant calculation
				if (q_pointer > 0) {
					// Jump from the last breakpoint kmer
					int l1 = q_bp_intv[q_pointer-1].l, r1 = q_bp_intv[q_pointer-1].r;
					int max_score = -INF, max_id = -1;
					assert(j > 0);
					int low = 0, high = j-1, ans = -1;
					while (low <= high) {
						int mid = (low + high) >> 1;
						if (dp[i][mid].cj >= l1) {
							ans = mid;
							high = mid - 1;
						} else {
							low = mid + 1;
						}
					}
					if (ans != -1) {
						assert(dp[i][ans].cj == l1);
						for (int k = ans; k < j; k++) {
							if (dp[i][k].cj > r1) break;
							assert(dp[i][k].cj >= l1 and dp[i][k].cj <= r1);
							int tmp = max(dp[i][k].H + DEL_DUP_O, dp[i][k].B2) + DEL_DUP_E;
							if (tmp > max_score) {
								max_score = tmp;
								max_id = k;
							}
						}
						dp[i][j].B2 = max_score;
						if (dp[i][j].B2 > dp[i][j].H) {
							dp[i][j].H = dp[i][j].B2;
							dp[i][j].pi = i;
							dp[i][j].pj = max_id;
						}
					}
				}

				if (t_pointer > 0) {
					int l1 = t_bp_intv[t_pointer-1].l, r1 = t_bp_intv[t_pointer-1].r;
					int max_score = -INF, max_id = -1;
					for (int k = l1; k <= r1; k++) {
						int low = 0, high = dp[k].size()-1, ans = -1;
						while (low <= high) {
							int mid = (low + high) >> 1;
							if (dp[k][mid].cj > curr_j) {
								high = mid - 1;
							} else if (dp[k][mid].cj < curr_j) {
								low = mid + 1;
							} else {
								ans = mid;
								break;
							}
						}
						assert(ans != -1);
						if (ans != -1) {
							int tmp = max(dp[k][ans].H + DEL_DUP_O, dp[k][ans].B1) + DEL_DUP_E;
							if (tmp > max_score) {
								max_score = tmp;
								max_id = k;
							}
						}
					}
					dp[i][j].B1 = max_score;
					if (dp[i][j].B1 > dp[i][j].H) {
						dp[i][j].H = dp[i][j].B1;
						dp[i][j].pi = max_id;
						dp[i][j].pj = j;
					}
				}
			}
		}
	}
	cout << "DSI result " << dp[t_len].back().H << endl;
	int raw_size = (t_len + 1) * (q_len + 1);
	int part_size = 0;
	for (int i = 0; i <= t_len; i++) {
		part_size += dp[i].size();
	}
	cout << "Band size: " << bw << endl;
	cout << "Calculation rate [%%]: " << 100.0 * part_size / raw_size << endl;
}

static vector<int> input_break_points(const char *fn)
{
	ifstream in(fn);
	assert(in.is_open());
	vector<int> ret;
	int num;
	while (in >> num) {
		ret.push_back(num);
	}
	in.close();
	return ret;
}

int usage(const ZigOptions &o) {
	fprintf(stderr, "Usage: zigalign [options] seq1.fa seq2.fa > aln.paf\n");
	fprintf(stderr, "  Regular Scoring options:\n");
	fprintf(stderr, "    -A [INT]  match score [%d]\n", o.mat_score);
	fprintf(stderr, "    -B [INT]  mismatch penalty [%d]\n", o.mis_pen);
	fprintf(stderr, "    -O [INT]  open gap(indel) penalty [%d]\n", o.gap_o);
	fprintf(stderr, "    -E [INT]  extend gap penalty [%d]\n", o.gap_e);
	fprintf(stderr, "    -J [INT]  open delete duplication penalty [%d]\n", o.del_dup_o);
	fprintf(stderr, "    -j [INT]  extend delete duplication penalty [%d]\n", o.del_dup_e);
	fprintf(stderr, "    -v [STR]  output file prefix\n");
	fprintf(stderr, "              matrix backtrace paths will be written to prefix_s1.txt, prefix_s2.txt and prefix_p.txt\n");
	fprintf(stderr, "              use vis_self.py to render prefix_s1.txt and prefix_s2.txt\n");
	fprintf(stderr, "              use vis_pair.py to render prefix_p.txt\n");
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

void process_long(const char *fn) {
	pair<string, string> pair = input_fasta_seq(fn);
	cout << pair.first << endl;
	cout << pair.second.length() << endl;
	const char *seq = pair.second.data();
	int n = pair.second.length();

	// NOTE: I got different results from n = 20000; n might affect the accuracy of breakpoints
	// NOTE: serial execution because of cyclic repeats
	const int part_len = 40000;
	ZigOptions opt;
	opt.close_tr_pen = -15; // It will increase the continuity of zigzag lines

	vector<int> all_bps;
	all_bps.push_back(0);
	int offset = 0;
	int cnt = 0;
	while (offset < n) {
		int len = min(n - offset, part_len);
		fprintf(stderr, "offset=%d, length=%d\n", offset, len);
		string vis_fn = "pairwise/hors/hor" + to_string(++cnt) + ".txt";
		vector<int> bps = self_alignment(opt, len, seq + offset, vis_fn);
		for (int i = 1; i < bps.size(); i++) {
			all_bps.push_back(bps[i] + offset);
		}

		// Continue from the penultimate breakpoints
		if (bps.size() > 2) {
			offset += bps[bps.size() - 2] + 1;
		} else {
			offset += len;
		}
	}
	fprintf(stdout, "Stitched %ld breakpoints\n", all_bps.size());
	for (int b: all_bps) {
		fprintf(stdout, "%d ", b);
	}
	fprintf(stdout, "\n");
}

void dev_pairwise(const char *fn1, const char *fn2) {
	pair<string, string> pair1 = input_fasta_seq(fn1);
	pair<string, string> pair2 = input_fasta_seq(fn2);
	string name1 = pair1.first, seq1 = pair1.second;
	string name2 = pair2.first, seq2 = pair2.second;
	int t_len = seq1.length(), q_len = seq2.length();
	cout << t_len << "\t" << q_len << endl;
	const char *t = seq1.data(), *q = seq2.data();
	ZigOptions opt;
	vector<int> t_bp = self_alignment(opt, t_len, t, "test/dev1.txt");
	vector<int> q_bp = self_alignment(opt, q_len, q, "test/dev2.txt");

	align_with_dups2(opt, t_len, t, t_bp, q_len, q, q_bp);
	align_with_dups(opt, t_len, t, t_bp, q_len, q, q_bp);
}

int main(int argc, char *argv[]) {
	MAX_BAND_WIDTH = atoi(argv[3]);
	dev_pairwise(argv[1], argv[2]);

	// process_long(argv[1]);
	return 1;

	double ctime = cputime(), rtime = realtime();
	ZigOptions opt;
	if (argc == 1) return usage(opt);
	int c;
	while ((c = getopt(argc, argv, "A:B:O:E:J:j:u:d:p:a:b:o:e:k:v:")) >= 0) {
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
			case 'k':
				opt.band_width = abs(str2int(optarg));
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
		// align_with_dups(opt, argv[optind], argv[optind+1]);
	} else {
		fprintf(stderr, "Two FASTA files are required\n");
		return 1;
	}
	fprintf(stderr, "Program finishes in %.3f CPU seconds, %.3f real seconds\n", cputime()-ctime, realtime()-rtime);
	return 0;
}