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
#include <omp.h>

#include "utils.h"
using namespace std;

// Zigalign parameters
struct ZigOptions {
	// SI Scoring matrix
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
	int n_threads;
	// Visualization
	const char *vis_prefix;

	ZigOptions() {
		// Scoring parameters should be adjusted by duplication and variation rate
		mat_score = 1;
		mis_pen = -12; // Increased penalty for small variants, i.e., duplication indels are preferred.
		gap_o = -18;
		gap_e = -3;
		min_unit_size = 5;
		tr_mat_score = 2;
		tr_mis_pen = -3;
		tr_gap_o = -3;
		tr_gap_e = -1;
		open_tr_pen = -2;
		close_tr_pen = -6;
		del_dup_o = -6;
		del_dup_e = -2;
		band_width = 50;
		n_threads = 8;
		vis_prefix = nullptr;
	}
};

const int INF = 100000000;

const int PART_LEN = 40000;
const int STEP_LEN = 20000;

// Move them to code block
const int NORMAL = 0;
const int START_REP = 1;
const int NEW_COPY = 2;
const int WITHIN_REP = 3;
const int END_REP = 4;

int DEBUG = 0;

struct Dp1Cell {
	int D, E, F, H; // D: score for duplication events
	int beg, end; // The range of duplication or new copy
	Dp1Cell() {
		D = E = F = H = -INF;
		beg = end = -1;
	}
};

// Compact structure for backtrace
struct Bt1Cell {
	uint8_t event;
	uint8_t pi;
	uint16_t pj;
	Bt1Cell() {
		event = -1;
		pi = -1;
		pj = -1;
	}
};

struct RepUnit {
	int qb, qe; // [qb, qe) is a copy of [tb, te)
	int match, mis, gap;

	RepUnit() {
		qb = qe = 0;
		match = mis = gap = 0;
	}
};

struct TandemGroup {
	int tb, te; // [tb, te) is the template
	vector<RepUnit> units;

	TandemGroup() {
		tb = te = 0;
	}
};

struct Interval {
	int l, r; // [l, r)
};

struct TandemRepeat {
	Interval temp; // Templates
	vector<Interval> units;
};

const uint16_t SA_MAX_LEN = 50000;

 void self_alignment2(const ZigOptions &o, uint16_t n, const char *seq, vector<vector<Bt1Cell>> &bt)
{
	if (n > SA_MAX_LEN) {
		fprintf(stderr, "Sequence length %d exceeds the limit of self alignment %d\n", n, SA_MAX_LEN);
		abort();
	}

	const int MAT_SCORE = o.mat_score;
	const uint16_t MIN_UNIT = o.min_unit_size;
	const int OPEN_TR = o.open_tr_pen;
	const int CLOSE_TR = o.close_tr_pen;
	const int TR_MAT_SCORE = o.tr_mat_score;
	const int TR_MIS_PEN = o.tr_mis_pen;
	const int TR_GAP_O = o.tr_gap_o;
	const int TR_GAP_E = o.tr_gap_e;

	if (MAT_SCORE >= TR_MAT_SCORE) {
		fprintf(stderr, "Warning: setting higher score for diagonal match can't detect tandem repeats\n");
	}

	if (n < MIN_UNIT) {
		bt.resize(n + 1);
		for (int i = 0; i <= n; i++) {
			bt[i].resize(i + 1);
			if (i > 0) {
				bt[i][i].event = NORMAL;
				bt[i][i].pi = 1;
				bt[i][i].pj = i-1;
			}
		}
	}

	// Memory allocation
	vector<Dp1Cell> prev_dp(n + 1);
	vector<Dp1Cell> curr_dp(n + 1);
 	bt.resize(n + 1);
	for (int i = 0; i <= n; i++) {
		bt[i].resize(i + 1);
	}

	// Initialization
	prev_dp[0].H = 0;
	prev_dp[0].beg = 1;
	for (int i = 1; i <= MIN_UNIT; i++) {
		bt[i][i].event = NORMAL;
		bt[i][i].pi = 1;
		bt[i][i].pj = i-1;
	}
	prev_dp[MIN_UNIT].H = MIN_UNIT * MAT_SCORE;
	prev_dp[MIN_UNIT].beg = 1;

	// Main loop
	for (int i = MIN_UNIT + 1; i <= n; i++) {
		for (int j = 0; j <= i; j++) curr_dp[j] = Dp1Cell();
		curr_dp[i].H = prev_dp[i-1].H + MAT_SCORE;
		bt[i][i].event = NORMAL;
		bt[i][i].pi = 1;
		bt[i][i].pj = i-1;
		curr_dp[i].beg = prev_dp[i-1].beg;

		int max_value = (i-1) * MAT_SCORE; // From non-repetitive region (NOTE: decreased score)
		int t_beg = prev_dp[i-1].beg; // To prevent illegal path
		int t_end = i - 1;
		uint8_t event = START_REP;
		for (int j = 1; j < i-1; j++) {
			// Start new copy at the end of duplications
			if (prev_dp[j].end == j and prev_dp[j].H > max_value) {
				max_value = prev_dp[j].H;
				t_beg = prev_dp[j].beg;
				t_end = j;
				event = NEW_COPY;
			}
		}

		// FIXME: maximum value should be chosen from overlapping template intervals
		// D transfer
		for (int j = t_beg; j <= t_end - MIN_UNIT + 1; j++) {
			int tmp = (seq[i-1] == seq[j-1]) ? TR_MAT_SCORE : TR_MIS_PEN;
			int pen = (event == START_REP) ? OPEN_TR : 0; // No penalty for new copy
			curr_dp[j].D = max_value + tmp + pen;
			curr_dp[j].beg = j;
			curr_dp[j].end = t_end;
			bt[i][j].event = event;
			bt[i][j].pi = 1;
			bt[i][j].pj = t_end;
		}
		for (int j = t_end - MIN_UNIT + 2; j <= i - MIN_UNIT; j++) {
			int tmp = (seq[i-1] == seq[j-1]) ? TR_MAT_SCORE : TR_MIS_PEN;
			curr_dp[j].D = (i - 1) * MAT_SCORE + tmp + OPEN_TR;
			curr_dp[j].beg = j;
			curr_dp[j].end = i - 1;
			bt[i][j].event = START_REP;
			bt[i][j].pi = 1;
			bt[i][j].pj = i - 1;
		}

		// Sub matrix of repetition
		for (int j = 1; j < i; j++) {
			int v_score = -INF, h_score = -INF, d_score = -INF;
			if (j >= prev_dp[j].beg and j <= prev_dp[j].end) {
				v_score = max(max(prev_dp[j].D, prev_dp[j].H) + TR_GAP_O, prev_dp[j].E) + TR_GAP_E;
				curr_dp[j].E = v_score;
			}
			if (j >= curr_dp[j-1].beg and j <= curr_dp[j-1].end) {
				h_score = max(max(curr_dp[j-1].D, curr_dp[j-1].H) + TR_GAP_O, curr_dp[j-1].F) + TR_GAP_E;
				curr_dp[j].F = h_score;
			}
			if (j >= prev_dp[j-1].beg and j <= prev_dp[j-1].end) {
				int tmp = (seq[i-1] == seq[j-1] ? TR_MAT_SCORE : TR_MIS_PEN);
				d_score = max(prev_dp[j-1].D, prev_dp[j-1].H) + tmp;
			}

			// Here, >= prefers continue an existing copy, so it can generate longer tandem repeats
			if (v_score >= curr_dp[j].D and v_score > curr_dp[j].H) {
				curr_dp[j].H = v_score;
				curr_dp[j].beg = prev_dp[j].beg;
				curr_dp[j].end = prev_dp[j].end;
				bt[i][j].event = WITHIN_REP;
				bt[i][j].pi = 1;
				bt[i][j].pj = j;
			}
			if (h_score >= curr_dp[j].D and h_score > curr_dp[j].H) {
				curr_dp[j].H = h_score;
				curr_dp[j].beg = curr_dp[j-1].beg;
				curr_dp[j].end = curr_dp[j-1].end;
				bt[i][j].event = WITHIN_REP;
				bt[i][j].pi = 0;
				bt[i][j].pj = j-1;
			}
			if (d_score >= curr_dp[j].D and d_score > curr_dp[j].H) {
				curr_dp[j].H = d_score;
				curr_dp[j].beg = prev_dp[j-1].beg;
				curr_dp[j].end = prev_dp[j-1].end;
				bt[i][j].event = WITHIN_REP;
				bt[i][j].pi = 1;
				bt[i][j].pj = j-1;
			}
		}

		// B transfer
		max_value = -INF;
		t_beg = -1;
		t_end = -1;
		for (int j = 1; j < i; j++) {
			// Only return to diagonal if sub-matrix reaches the lower-right corner
			if (curr_dp[j].H > max_value and curr_dp[j].end == j) {
				max_value = curr_dp[j].H;
				t_beg = curr_dp[j].beg;
				t_end = curr_dp[j].end;
			}
		}
		if (max_value + CLOSE_TR > curr_dp[i].H) {
			curr_dp[i].H = max_value + CLOSE_TR;
			curr_dp[i].beg = t_beg; // This variable is reused to prevent illegal path
			bt[i][i].event = END_REP;
			bt[i][i].pi = 0;
			bt[i][i].pj = t_end;
		}

		swap(prev_dp, curr_dp);
	}
}

struct RepInterval {
	int beg, end; // 1-based [beg, end)
	int mis, gap; // #mismatches and #gaps

	RepInterval() {
		beg = end = -1;
		mis = gap = 0;
	}

	bool operator < (const RepInterval &b) const {
		return (this->beg != b.beg) ? this->beg < b.beg : this->end < b.end;
	}
};

vector<RepInterval> get_reps(uint16_t n, const char *seq, const vector<vector<Bt1Cell>> &bt)
{
	int ti = n, tj = n;
	vector<RepInterval> reps;
	while (ti > 0 and tj > 0) {
		Bt1Cell t = bt[ti][tj];
		if (t.event == END_REP) {
			assert(ti == tj); // Only main diagonal closes repetitions
			while (t.event != START_REP) {
				ti = t.pi == 1 ?ti-1 :ti;
				tj = t.pj;
				t = bt[ti][tj]; // Lower-right corner of the sub-matrix
				RepInterval u;
				u.end = ti + 1;
				while (t.event != NEW_COPY and t.event != START_REP) {
					if (t.pi == 1 and tj == t.pj + 1) {
						u.mis += (seq[ti - 1] != seq[tj - 1]);
					} else {
						u.gap++;
					}
					ti = t.pi == 1 ?ti-1 :ti;
					tj = t.pj;
					t = bt[ti][tj];
				}
				// Pointer is now at the upper-left corner
				u.mis += (seq[ti - 1] != seq[tj - 1]);
				u.beg = ti;
				reps.push_back(u);
			}
			// Template for all copy units above
			RepInterval tem;
			tem.beg = tj;
			tem.end = ti;
			tem.mis = tem.gap = -1;
			reps.push_back(tem);
		}
		ti = t.pi == 1 ?ti-1 :ti;
		tj = t.pj;
	}

	if (reps.empty()) return reps;

	// de-overlap repeat units (tandem repeats should be neatly stacked)
	sort(reps.begin(), reps.end());
	int new_size = 0;
	for (int i = 0; i < reps.size(); ) {
		const RepInterval &x = reps[i];
		bool kept = true;
		int j = i + 1;
		for (; j < reps.size(); j++) {
			const RepInterval &y = reps[j];
			if (y.beg >= x.end) break;
			// Keep the longer one
			if (y.end - y.beg >= x.end - x.beg) {
				kept = false;
				break;
			}
		}
		if (kept) reps[new_size++] = x;
		i = j;
	}
	reps.resize(new_size);
	// TODO: merge or split repeats

	// for (const RepInterval &r: reps) {
	// 	printf("[%d, %d) mis=%d, gap=%d\n", r.beg, r.end, r.mis, r.gap);
	// }
	return reps;
}

vector<int> self_alignment(const ZigOptions &o, uint16_t n, const char *seq, const string &vis_fn = "")
{
	vector<vector<Bt1Cell>> bt;
	self_alignment2(o, n, seq, bt);
	if (not vis_fn.empty()) {
		ofstream out(vis_fn);
		assert(out.is_open());
		int ti = n, tj = n, te = -1;
		while (ti > 0 and tj > 0) {
			const Bt1Cell &c = bt[ti][tj];
			if (c.event != te) {
				out << ti << "\t" << tj << "\t" << (int)c.event << endl;
			}
			ti = (c.pi == 1) ?ti-1 :ti;
			tj = c.pj;
			te = c.event;
		}
		out << 0 << "\t" << 0 << "\t" << te << endl;
		out.close();
	}

	// Trace back the optimal path to find tandem repeats
	int ti = n, tj = n;
	vector<TandemGroup> groups;
	while (ti > 0 and tj > 0) {
		Bt1Cell &t = bt[ti][tj];
		if (t.event == END_REP) {
			assert(ti == tj); // Only main diagonal closes repetitions
			TandemGroup g;
			while (t.event != START_REP) {
				ti = t.pi == 1 ?ti-1 :ti;
				tj = t.pj;
				t = bt[ti][tj]; // Lower-right corner of the sub-matrix
				RepUnit u;
				u.qe = ti + 1;
				while (t.event != NEW_COPY and t.event != START_REP) {
					if (t.pi == 1 and tj == t.pj + 1) {
						if (seq[ti - 1] == seq[tj - 1]) u.match++;
						else u.mis++;
					} else {
						u.gap++;
					}
					ti = t.pi == 1 ?ti-1 :ti;
					tj = t.pj;
					t = bt[ti][tj];
				}
				// Pointer is now at the upper-left corner
				u.qb = ti;
				g.units.push_back(u);
				if (seq[ti - 1] == seq[tj - 1]) {
					u.match++;
				} else {
					u.mis++;
				}
			}
			g.tb = tj;
			g.te = ti;
			reverse(g.units.begin(), g.units.end());
			groups.push_back(g);
		}
		ti = t.pi == 1 ?ti-1 :ti;
		tj = t.pj;
	}
	// It is difficult and inaccurate to directly merge the tandem groups.
	// An easy way is to break them down to non-overlapping intervals, then merge them based on length and base identity.
	// But it is more time-consuming if base identity is calculated by DP.

	reverse(groups.begin(), groups.end());
	for (const TandemGroup &g: groups) {
		assert(g.units.size() > 0);
		assert(g.units[0].qb >= g.te);
		printf("Template: [%d, %d)\n", g.tb, g.te);
		// Tandem repeats must be stacked
		for (int i = 1; i < g.units.size(); i++) {
			assert(g.units[i].qb >= g.units[i-1].qe);
		}
		printf("Units: ");
		int total_match = 0, total_mis = 0, total_gap = 0;
		for (const RepUnit &u: g.units) {
			printf(" [%d, %d)", u.qb, u.qe);
			total_match += u.match;
			total_mis += u.mis;
			total_gap += u.gap;
		}
		printf("\n");
		double ave_match = 100.0 * total_match / g.units.size() / (g.te - g.tb);
		double ave_mis = 1.0 * total_mis / g.units.size();
		double ave_gap = 1.0 * total_gap / g.units.size();
		printf("Identity: %.2f %%, ave_mis: %.2f, ave_gap: %.2f\n", ave_match, ave_mis, ave_gap);
	}
	for (int i = 1; i < groups.size(); i++) {
		assert(groups[i].units.back().qe > groups[i-1].units.back().qe);
	}

	vector<int> ret;
	if (ret.empty()) return ret;
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

vector<int> global_pairwise(const ZigOptions &opt,
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
			int M = dp[i-1][j-1].H + (t[i-1] == q[j-1] ? MAT_SCORE : MIS_PEN);
			if (M > dp[i][j].H) {
				dp[i][j].H = M;
				dp[i][j].pi = i-1;
				dp[i][j].pj = j-1;
			}
			dp[i][j].E = max(dp[i-1][j].H + GAP_O, dp[i-1][j].E) + GAP_E;
			if (dp[i][j].E > dp[i][j].H) {
				dp[i][j].H = dp[i][j].E;
				dp[i][j].pi = i-1;
				dp[i][j].pj = j;
			}
			dp[i][j].F = max(dp[i][j-1].H + GAP_O, dp[i][j-1].F) + GAP_E;
			if (dp[i][j].F > dp[i][j].H) {
				dp[i][j].H = dp[i][j].F;
				dp[i][j].pi = i;
				dp[i][j].pj = j-1;
			}

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
					// Redundant calculation
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
		}
	}

	cout << "DIS alignment score: " << dp[t_len][q_len].H << endl;

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

	if (opt.vis_prefix) {
		string pair_vis_fn = string(opt.vis_prefix) + "_p.txt";
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
	return cv;
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

vector<int> banded_pairwise(const ZigOptions &opt,
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

	const int bw = opt.band_width; // Bandwidth
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
					assert(ans != -1);
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
					int max_score = -INF, max_id = -1, mate_id;
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
								mate_id = ans;
							}
						}
					}
					dp[i][j].B1 = max_score;
					if (dp[i][j].B1 > dp[i][j].H) {
						dp[i][j].H = dp[i][j].B1;
						dp[i][j].pi = max_id;
						dp[i][j].pj = mate_id;
					}
				}
			}
		}
	}
	long raw_size = (long)(t_len + 1) * (q_len + 1);
	long part_size = 0;
	for (int i = 0; i <= t_len; i++) {
		part_size += dp[i].size();
	}
	cout << "DSI alignment score: " << dp[t_len].back().H << endl;
	cout << "Band size: " << bw << endl;
	cout << "Calculation rate [%]: " << 100.0 * part_size / raw_size << endl;

	// CIGAR generation
	int ti = t_len, tj = dp[t_len].size()-1;
	int ci = dp[ti][tj].ci, cj = dp[ti][tj].cj;
	assert(ci == t_len and cj == q_len);
	int del_n = 0, ins_n = 0, mat_n = 0, mis_n = 0, dup_n = 0;
	vector<int> cv;
	// All operations occur on target sequence _t_
	const int COP_M = 0;
	const int COP_I = 1;
	const int COP_D = 2;
	const int DUP_I = 3;
	const int DUP_D = 4;
	const string OP_CHAR = "MIDID";
	int op_type = -1, last_op = -1, op_cnt = 0;
	string ext_t, ext_q;
	while (ti > 0 and tj > 0) {
		const Dp3Cell &p = dp[ti][tj];
		assert(p.pi != -1 and p.pj != -1);
		const Dp3Cell &g = dp[p.pi][p.pj];
		int pi = g.ci, pj = g.cj;
		if (pi == ci - 1 and pj == cj) { // Deletion from _t_
			op_type = COP_D;
			del_n++;
			ext_t += t[ci-1];
			ext_q += '-';
		} else if (pi == ci and pj == cj - 1) { // Insertion into _t_
			op_type = COP_I;
			ins_n++;
			ext_t += '-';
			ext_q += q[cj-1];
		} else if (pi == ci - 1 and pj == cj - 1) { // Match/Mismatch
			op_type = COP_M;
			if (t[ci-1] == q[cj-1]) mat_n++;
			else mis_n++;
			ext_t += t[ci-1];
			ext_q += q[cj-1];
		} else {
			if (last_op != -1) cv.push_back(op_cnt << 4 | last_op);
//			fprintf(stderr, "%d %d -> %d %d\n", ti, tj, t.pi, t.pj);
			if (ci == pi) { // Duplication insertion into _t_
				op_type = DUP_I;
				op_cnt = cj - pj;
				for (int j = cj; j > pj; j--) {
					ext_t += '+';
					ext_q += q[j-1];
				}
			} else { // Duplication deletion from _t_
				op_type = DUP_D;
				op_cnt = ci - pi;
				for (int i = ci; i > pi; i--) {
					ext_t += t[i-1];
					ext_q += '+';
				}
			}

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
		ci = g.ci;
		cj = g.cj;
	}
	if (last_op != -1) cv.push_back(op_cnt << 4 | last_op);
	// cout << ti << "\t" << tj << endl;
	if (ti > 0) {
		cv.push_back(ci << 4 | COP_D);
	}
	if (tj > 0) {
		cv.push_back(cj << 4 | COP_I);
	}
	reverse(cv.begin(), cv.end());
	reverse(ext_t.begin(), ext_t.end());
	reverse(ext_q.begin(), ext_q.end());
	fprintf(stderr, "%d deletions, %d insertions, %d matches, %d mismatches and %d duplication indels\n", del_n, ins_n, mat_n, mis_n, dup_n);

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
	return cv;
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

// FIXME: the stitching is inaccurate
vector<int> process_long(const ZigOptions &opt, int n, const char *seq) {
	// NOTE: I got different results from n = 20000; n might affect the accuracy of breakpoints
	// NOTE: serial execution because of cyclic repeats
	const int part_len = 40000;

	vector<int> all_bps;
	all_bps.push_back(0);
	int offset = 0;
	int cnt = 0;
	while (offset < n) {
		int len = min(n - offset, part_len);
		fprintf(stderr, "cnt=%d, offset=%d, length=%d\n", cnt+1, offset, len);
		string vis_fn = "pairwise/hors/hor_b" + to_string(++cnt) + ".txt";
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
	fprintf(stderr, "Stitched %ld breakpoints\n", all_bps.size());
	return all_bps;
}

const int MAX_SEQ_LEN = 40000;
void align_with_dups(const ZigOptions &opt, const char *fn1, const char *fn2) {
	pair<string, string> pair1 = input_fasta_seq(fn1);
	pair<string, string> pair2 = input_fasta_seq(fn2);
	string name1 = pair1.first, seq1 = pair1.second;
	string name2 = pair2.first, seq2 = pair2.second;
	int t_len = seq1.length(), q_len = seq2.length();
	const char *t = seq1.data(), *q = seq2.data();
	if (t_len > MAX_SEQ_LEN or q_len > MAX_SEQ_LEN) {
		vector<int> t_bp = (t_len > MAX_SEQ_LEN) ?process_long(opt, t_len, t) :self_alignment(opt, t_len, t, "");
		vector<int> q_bp = (q_len > MAX_SEQ_LEN) ?process_long(opt, q_len, q) :self_alignment(opt, q_len, q, "");
		// It is slow because the partial matrix is still large
		vector<int> cv2 = banded_pairwise(opt, t_len, t, t_bp, q_len, q, q_bp);
		paf_format(name2, seq2, name1, seq1, cv2);
	} else {
		// banded_pairwise(opt, t_len, t, q_len, q);
		string t_vis_fn = "", q_vis_fn = "";
		if (opt.vis_prefix) {
			t_vis_fn = string(opt.vis_prefix) + "_s1.txt";
			q_vis_fn = string(opt.vis_prefix) + "_s2.txt";
		}
		vector<int> t_bp = self_alignment(opt, t_len, t, t_vis_fn);
		vector<int> q_bp = self_alignment(opt, q_len, q, q_vis_fn);
		vector<int> cv2 = banded_pairwise(opt, t_len, t, t_bp, q_len, q, q_bp);
		paf_format(name2, seq2, name1, seq1, cv2);
	}
}

// Recursive search is not efficient.
inline int us_find(vector<int> &sid, int u) {
	if (sid[u] != u) {
		sid[u] = us_find(sid, sid[u]);
	}
	return sid[u];
}

// Rank is not used.
inline void us_union(vector<int> &sid, int u, int v) {
	int su = us_find(sid, u);
	int sv = us_find(sid, v);
	if (su != sv) {
		sid[sv] = su;
	}
}

vector<RepInterval> shift_reps(int n, const vector<RepInterval> &reps, vector<vector<Bt1Cell>> &bt, int pos)
{
	assert(pos < n);
	vector<bool> vis(n + 1, false);
	for (const RepInterval &r: reps) {
		if (r.mis != -1) {
			for (int i = r.beg; i < r.end; i++) {
				vis[i] = true;
			}
		}
	}
	vector<int> sid(n + 1);
	for (int i = 1; i <= n; i++) {
		sid[i] = i;
	}
	int ti = n, tj = n;
	while (ti > 0 and tj > 0) {
		const Bt1Cell &c = bt[ti][tj];
		// Project copies onto template
		if (vis[ti]) us_union(sid, ti, tj);
		ti = (c.pi == 1) ?ti-1 :ti;
		tj = c.pj;
	}

	vector<RepInterval> ret;
	int shift = pos; // Truncated leading bases
	int k = 0;
	for (; k < reps.size(); k++) {
		const RepInterval &r = reps[k];
		if (r.beg > shift) {
			break;
		}
	}
	if (k == 0 or reps[k-1].end <= shift) {
		// Truncated position is not in any tandem repeat
		for (; k < reps.size(); k++) {
			ret.push_back(reps[k]);
		}
	} else {
		assert(reps[k-1].beg <= shift and reps[k-1].end > shift);
		RepInterval v;
		v.beg = shift;
		int root = us_find(sid, shift);
		for (; k < reps.size(); k++) {
			const RepInterval &r = reps[k];
			bool found = false;
			for (int i = r.beg; i < r.end; i++) {
				if (us_find(sid, i) == root) {
					found = true;
					v.end = i;
					ret.push_back(v);
					v.beg = i;
					break;
				}
			}
			if (not found) break; // Till another group of repeats
		}
		// The rest of repeats are not affected by the truncation
		for (; k < reps.size(); k++) {
			const RepInterval &r = reps[k];
			ret.push_back(r);
		}
	}
	// printf("Shifted reps:\n");
	// for (const RepInterval &r: ret) {
	// 	printf("[%d, %d) mis=%d, gap=%d\n", r.beg, r.end, r.mis, r.gap);
	// }
	return ret;
}

void process_long2(const ZigOptions &opt, int n, const char *seq)
{
	assert(opt.n_threads > 0);
	int n_threads = opt.n_threads;
	int offset = 0;
	vector<vector<Bt1Cell>> bt_bin[n_threads];
	vector<RepInterval> reps_bin[n_threads];
	vector<RepInterval> all_reps;

	int batch_id = 0;
	omp_set_num_threads(n_threads);
	while (offset < n) {
		double r_start = realtime(), c_start = cputime();
		int guard = offset + (n_threads - 1) * STEP_LEN + PART_LEN;
		// fprintf(stderr, "offset = %d, guard = %d\n", offset, guard);
		if (guard >= n) break;
		// Multi-threading DP
		#pragma omp parallel for
		for (int i = 0; i < n_threads; i++) {
			const char *s = seq + offset + i * STEP_LEN;
			self_alignment2(opt, PART_LEN, s, bt_bin[i]);
			reps_bin[i] = get_reps(PART_LEN, s, bt_bin[i]);
		}

		// for (int i = 0; i < n_threads; i++) {
		// 	const vector<RepInterval> &bin = reps_bin[i];
		// 	fprintf(stderr, "thread %d\n", i);
		// 	for (int j = 0; j < bin.size(); j++) {
		// 		fprintf(stderr, "@%d %d %d\n", j, bin[j].beg, bin[j].end);
		// 	}
		// 	fprintf(stderr, "\n");
		// }

		// TODO: use thread pipeline
		// Single-thread Stitching
		int shift_pos = 0;
		int added_reps = 0;
		for (int i = 0; i < n_threads; i++) {
			// fprintf(stderr, "i=%d, shift = %d\n", i, shift_pos);
			if (shift_pos) {
				reps_bin[i] = shift_reps(PART_LEN, reps_bin[i], bt_bin[i], shift_pos);
			}
			for (const RepInterval &r: reps_bin[i]) {
				RepInterval x = r;
				x.beg += offset + i * STEP_LEN;
				x.end += offset + i * STEP_LEN;
				all_reps.push_back(x);
				added_reps++;
				if (r.end >= STEP_LEN) {
					shift_pos = r.end - STEP_LEN;
					// Find the repeat that crosses the next bin to connect them
					if (i + 1 < n_threads and reps_bin[i + 1].size() > 0) {
						bool ok = (reps_bin[i+1][0].beg <= shift_pos and reps_bin[i+1][0].end >= shift_pos);
						if (not ok) continue;
					}
					break;
				}
			}
		}

		// Find the repeat crossing the boundary
		int old_offset = offset;
		bool found = false;
		for (const RepInterval &r: reps_bin[n_threads - 1]) {
			if (r.end >= STEP_LEN) {
				found = true;
				offset += (n_threads - 1) * STEP_LEN + r.end;
				break;
			}
		}
		if (not found) offset = guard;
		fprintf(stderr, "Batch %d, offset: %d, added_reps: %d, real_time: %.2f, CPU_time: %.2f\n",
			++batch_id, old_offset, added_reps, realtime() - r_start, cputime() - c_start);
	}

	// Process the last batch
	int m = n - offset;
	if (m <= PART_LEN) {
		n_threads = 1;
	} else {
		n_threads = m / STEP_LEN;
	}
	double r_start = realtime(), c_start = cputime();
	#pragma omp parallel for num_threads(n_threads)
	for (int i = 0; i < n_threads; i++) {
		int os = offset + i * STEP_LEN;
		const char *s = seq + os;
		int len = min(PART_LEN, n - os);
		self_alignment2(opt, len, s, bt_bin[i]);
		reps_bin[i] = get_reps(len, s, bt_bin[i]);
	}
	int shift_pos = 0;
	int added_reps = 0;
	for (int i = 0; i < n_threads; i++) {
		int os = offset + i * STEP_LEN;
		int len = min(PART_LEN, n - os);
		if (shift_pos) {
			reps_bin[i] = shift_reps(len, reps_bin[i], bt_bin[i], shift_pos);
		}
		for (const RepInterval &r: reps_bin[i]) {
			RepInterval x = r;
			x.beg += offset + i * STEP_LEN;
			x.end += offset + i * STEP_LEN;
			all_reps.push_back(x);
			added_reps++;
			if (i != n_threads - 1 and r.end >= STEP_LEN) {
				shift_pos = r.end - STEP_LEN;
				break;
			}
		}
	}
	fprintf(stderr, "Batch %d, offset: %d, length: %d, threads: %d, added_reps: %d, real_time: %.2f, CPU_time: %.2f\n",
		++batch_id, offset, m, n_threads, added_reps, realtime() - r_start, cputime() - c_start);

	for (const RepInterval &r: all_reps) {
		printf("[%d, %d), length=%d\n", r.beg, r.end, r.end - r.beg);
	}
}

void test(const char *fn1, int offset, int n) {
	pair<string, string> pair1 = input_fasta_seq(fn1);
	string name1 = pair1.first, seq1 = pair1.second;
	int t_len = seq1.length();
	const char *t = seq1.data() + offset;
	t_len -= offset;

	ZigOptions opt;
	process_long2(opt, t_len, t);
}

int main(int argc, char *argv[]) {
	test(argv[1], atoi(argv[2]), atoi(argv[3]));
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
		align_with_dups(opt, argv[optind], argv[optind+1]);
	} else {
		fprintf(stderr, "Two FASTA files are required\n");
		return 1;
	}
	fprintf(stderr, "Program finishes in %.3f CPU seconds, %.3f real seconds\n", cputime()-ctime, realtime()-rtime);
	return 0;
}