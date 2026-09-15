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
#include <omp.h>
#include "utils.h"
using namespace std;

const int INF = 100000000;
const int SA_MAX_LEN = 50000; // Self-alignment max length
const int PART_LEN = 40000;
const int MAX_UNIT_DIS = 100000; // TODO: set it to be a parameter
const double MIN_MATCH_RATIO = 0.6;

const bool DEBUG = false;

// Zigalign parameters
struct ZigOptions {
	// Scoring matrix for repeats in pairwise alignment
	int mat_score;
	int mis_pen;
	int gap_o;
	int gap_e;
	// Minimum repeat unit size
	int min_unit_size;
	// Scoring matrix for repeats in self-alignment
	int sa_mat_score;
	int sa_mis_pen;
	int sa_gap_o;
	int sa_gap_e;
	// Penalty for opening/closing tandem repeats
	int open_tr_pen;
	int close_tr_pen;
	// Penalty for deleting duplications
	int del_unit;
	int del_dup_o;
	int del_dup_e;
	// K-band size
	int band_width;
	// #Threads
	int n_threads;
	// Intermediate results
	const char *log_prefix;

	ZigOptions() {
		// Scoring parameters should be adjusted by duplication and variation rate
		mat_score = 1;
		mis_pen = -160; // Increased penalty for small variants, i.e., duplication indels are preferred.
		gap_o = -240;
		gap_e = -40;
		min_unit_size = 100;
		sa_mat_score = 2;
		sa_mis_pen = -3;
		sa_gap_o = -3;
		sa_gap_e = -1;
		open_tr_pen = -2;
		close_tr_pen = -6;
		del_unit = -20;
		del_dup_o = -6;
		del_dup_e = -2;
		band_width = 50;
		n_threads = 8;
		log_prefix = nullptr;
	}
};

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

struct RepInterval {
	int beg, end; // [beg, end)
	int mat, mis, gap; // #matches, mismatches and gaps

	RepInterval() {
		beg = end = -1;
		mat = mis = gap = 0;
	}

	bool operator < (const RepInterval &b) const {
		return (this->beg != b.beg) ? this->beg < b.beg : this->end < b.end;
	}
};

vector<RepInterval> self_alignment2(const ZigOptions &o, int n, const char *seq, vector<vector<Bt1Cell>> &bt)
{
	if (n > SA_MAX_LEN) {
		fprintf(stderr, "Sequence length %d exceeds the limit of self alignment %d\n", n, SA_MAX_LEN);
		abort();
	}

	const int MIN_UNIT = o.min_unit_size;
	const int MAT_SCORE = o.mat_score;
	const int OPEN_TR = o.open_tr_pen;
	const int CLOSE_TR = o.close_tr_pen;
	const int SA_MAT_SCORE = o.sa_mat_score;
	const int SA_MIS_PEN = o.sa_mis_pen;
	const int SA_GAP_O = o.sa_gap_o;
	const int SA_GAP_E = o.sa_gap_e;

	const int BT_NORMAL = 0;
	const int BT_START_REP = 1;
	const int BT_NEW_COPY = 2;
	const int BT_WITHIN_REP = 3;
	const int BT_END_REP = 4;

	if (MAT_SCORE >= SA_MAT_SCORE) {
		fprintf(stderr, "Warning: setting higher score for diagonal match can't detect tandem repeats\n");
	}

	if (n < MIN_UNIT) {
		bt.resize(n + 1);
		for (int i = 0; i <= n; i++) {
			bt[i].resize(i + 1);
			if (i > 0) {
				bt[i][i].event = BT_NORMAL;
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
		bt[i][i].event = BT_NORMAL;
		bt[i][i].pi = 1;
		bt[i][i].pj = i-1;
	}
	prev_dp[MIN_UNIT].H = MIN_UNIT * MAT_SCORE;
	prev_dp[MIN_UNIT].beg = 1;

	// Main loop
	for (int i = MIN_UNIT + 1; i <= n; i++) {
		for (int j = 0; j <= i; j++) curr_dp[j] = Dp1Cell();
		curr_dp[i].H = prev_dp[i-1].H + MAT_SCORE;
		bt[i][i].event = BT_NORMAL;
		bt[i][i].pi = 1;
		bt[i][i].pj = i-1;
		curr_dp[i].beg = prev_dp[i-1].beg;

		int max_value = (i-1) * MAT_SCORE; // From non-repetitive region (NOTE: decreased score)
		int t_beg = prev_dp[i-1].beg; // To prevent illegal path
		int t_end = i - 1;
		uint8_t event = BT_START_REP;
		for (int j = 1; j < i-1; j++) {
			// Start new copy at the end of duplications
			if (prev_dp[j].end == j and prev_dp[j].H > max_value) {
				max_value = prev_dp[j].H;
				t_beg = prev_dp[j].beg;
				t_end = j;
				event = BT_NEW_COPY;
			}
		}

		// FIXME: maximum value should be chosen from overlapping template intervals
		// D transfer
		for (int j = t_beg; j <= t_end - MIN_UNIT + 1; j++) {
			int tmp = (seq[i-1] == seq[j-1]) ? SA_MAT_SCORE : SA_MIS_PEN;
			int pen = (event == BT_START_REP) ? OPEN_TR : 0; // No penalty for new copy
			curr_dp[j].D = max_value + tmp + pen;
			curr_dp[j].beg = j;
			curr_dp[j].end = t_end;
			bt[i][j].event = event;
			bt[i][j].pi = 1;
			bt[i][j].pj = t_end;
		}
		for (int j = t_end - MIN_UNIT + 2; j <= i - MIN_UNIT; j++) {
			int tmp = (seq[i-1] == seq[j-1]) ? SA_MAT_SCORE : SA_MIS_PEN;
			curr_dp[j].D = (i - 1) * MAT_SCORE + tmp + OPEN_TR;
			curr_dp[j].beg = j;
			curr_dp[j].end = i - 1;
			bt[i][j].event = BT_START_REP;
			bt[i][j].pi = 1;
			bt[i][j].pj = i - 1;
		}

		// Sub matrix of repetition
		for (int j = 1; j < i; j++) {
			int v_score = -INF, h_score = -INF, d_score = -INF;
			if (j >= prev_dp[j].beg and j <= prev_dp[j].end) {
				v_score = max(max(prev_dp[j].D, prev_dp[j].H) + SA_GAP_O, prev_dp[j].E) + SA_GAP_E;
				curr_dp[j].E = v_score;
			}
			if (j >= curr_dp[j-1].beg and j <= curr_dp[j-1].end) {
				h_score = max(max(curr_dp[j-1].D, curr_dp[j-1].H) + SA_GAP_O, curr_dp[j-1].F) + SA_GAP_E;
				curr_dp[j].F = h_score;
			}
			if (j >= prev_dp[j-1].beg and j <= prev_dp[j-1].end) {
				int tmp = (seq[i-1] == seq[j-1] ? SA_MAT_SCORE : SA_MIS_PEN);
				d_score = max(prev_dp[j-1].D, prev_dp[j-1].H) + tmp;
			}

			// Here, >= prefers continue an existing copy, so it can generate longer tandem repeats
			if (v_score >= curr_dp[j].D and v_score > curr_dp[j].H) {
				curr_dp[j].H = v_score;
				curr_dp[j].beg = prev_dp[j].beg;
				curr_dp[j].end = prev_dp[j].end;
				bt[i][j].event = BT_WITHIN_REP;
				bt[i][j].pi = 1;
				bt[i][j].pj = j;
			}
			if (h_score >= curr_dp[j].D and h_score > curr_dp[j].H) {
				curr_dp[j].H = h_score;
				curr_dp[j].beg = curr_dp[j-1].beg;
				curr_dp[j].end = curr_dp[j-1].end;
				bt[i][j].event = BT_WITHIN_REP;
				bt[i][j].pi = 0;
				bt[i][j].pj = j-1;
			}
			if (d_score >= curr_dp[j].D and d_score > curr_dp[j].H) {
				curr_dp[j].H = d_score;
				curr_dp[j].beg = prev_dp[j-1].beg;
				curr_dp[j].end = prev_dp[j-1].end;
				bt[i][j].event = BT_WITHIN_REP;
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
			bt[i][i].event = BT_END_REP;
			bt[i][i].pi = 0;
			bt[i][i].pj = t_end;
		}

		swap(prev_dp, curr_dp);
	}

	// Backtrace
	int ti = n, tj = n;
	vector<RepInterval> reps;
	while (ti > 0 and tj > 0) {
		Bt1Cell t = bt[ti][tj];
		if (t.event == BT_END_REP) {
			assert(ti == tj); // Only main diagonal closes repetitions
			while (t.event != BT_START_REP) {
				ti = t.pi == 1 ?ti-1 :ti;
				tj = t.pj;
				t = bt[ti][tj]; // Lower-right corner of the sub-matrix
				RepInterval u;
				u.end = ti + 1;
				while (t.event != BT_NEW_COPY and t.event != BT_START_REP) {
					if (t.pi == 1 and tj == t.pj + 1) {
						u.mis += (seq[ti - 1] != seq[tj - 1]);
						u.mat += (seq[ti - 1] == seq[tj - 1]);
					} else {
						u.gap++;
					}
					ti = t.pi == 1 ?ti-1 :ti;
					tj = t.pj;
					t = bt[ti][tj];
				}
				// Pointer is now at the upper-left corner
				u.mis += (seq[ti - 1] != seq[tj - 1]);
				u.mat += (seq[ti - 1] == seq[tj - 1]);
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

	// Reset to 0-based index
	for (RepInterval &r: reps) {
		r.beg--;
		r.end--;
	}

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
	return reps;
}

struct WdpCell {
	int E, F, H;

	WdpCell() {
		E = F = H = -INF;
	}
};

vector<RepInterval> wraparound_dp(const int pat_len, const char *pat_seq, const int que_len, const char *que_seq)
{
	// Classical SI score matrix
	const int MAT_S = 1;
	const int MIS_P = -4;
	const int GAP_O = -6;
	const int GAP_E = -1;

	const int VERTICAL = 1;
	const int HORIZONTAL = 2;
	const int DIAGONAL = 3;
	const int WRAP_1 = 4;
	const int WRAP_2 = 5;

	vector<WdpCell> prev_dp( pat_len + 1);
	vector<WdpCell> curr_dp(pat_len + 1);
	vector<vector<uint8_t>> bt(que_len + 1);
	for (int i = 0; i <= que_len; i++) {
		bt[i].resize(pat_len + 1, 0);
	}

	prev_dp[0].H = 0;
	for (int j = 1; j <= pat_len; j++) {
		prev_dp[j].H = prev_dp[j].E = GAP_O + j * GAP_E;
		bt[0][j] = HORIZONTAL;
	}
	for (int i = 1; i <= que_len; i++) {
		// First pass
		curr_dp[0].H = GAP_O + i * GAP_E;
		bt[i][0] = VERTICAL;
		curr_dp[1].F = max(prev_dp[pat_len].H + GAP_O, prev_dp[pat_len].F) + GAP_E; // Wrap
		curr_dp[1].E = max(prev_dp[1].H + GAP_O, prev_dp[1].E) + GAP_E;
		int tmp = que_seq[i-1] == pat_seq[0] ?MAT_S :MIS_P;
		int d1 = prev_dp[0].H + tmp; // Leading deletions
		int d2 = prev_dp[pat_len].H + tmp; // Wrap
		curr_dp[1].H = -INF;
		if (curr_dp[1].F > curr_dp[1].H) {
			curr_dp[1].H = curr_dp[1].F;
			bt[i][1] = WRAP_1;
		}
		if (curr_dp[1].E > curr_dp[1].H) {
			curr_dp[1].H = curr_dp[1].E;
			bt[i][1] = VERTICAL;
		}
		if (d1 > curr_dp[1].H) {
			curr_dp[1].H = d1;
			bt[i][1] = DIAGONAL;
		}
		if (d2 > curr_dp[1].H) {
			curr_dp[1].H = d2;
			bt[i][1] = WRAP_1;
		}
		for (int j = 2; j <= pat_len; j++) {
			curr_dp[j].F = max(curr_dp[j-1].H + GAP_O, curr_dp[j-1].F) + GAP_E;
			curr_dp[j].E = max(prev_dp[j].H + GAP_O, prev_dp[j].E) + GAP_E;
			int d = prev_dp[j-1].H + (que_seq[i-1] == pat_seq[j-1] ?MAT_S : MIS_P);
			curr_dp[j].H = -INF;
			if (curr_dp[j].F > curr_dp[j].H) {
				curr_dp[j].H = curr_dp[j].F;
				bt[i][j] = HORIZONTAL;
			}
			if (curr_dp[j].E > curr_dp[j].H) {
				curr_dp[j].H = curr_dp[j].E;
				bt[i][j] = VERTICAL;
			}
			if (d > curr_dp[j].H) {
				curr_dp[j].H = d;
				bt[i][j] = DIAGONAL;
			}
		}
		// Second pass
		tmp = curr_dp[pat_len].H + GAP_O + GAP_E;
		if (tmp > curr_dp[1].F) {
			curr_dp[1].F = tmp;
		}
		if (curr_dp[1].F > curr_dp[1].H) {
			curr_dp[1].H = curr_dp[1].F;
			bt[i][1] = WRAP_2;
		}
		for (int j = 2; j <= pat_len; j++) {
			tmp = max(curr_dp[j-1].H + GAP_O, curr_dp[j-1].F) + GAP_E;
			curr_dp[j].F = max(tmp, curr_dp[j].F);
			if (curr_dp[j].F > curr_dp[j].H) {
				curr_dp[j].H = curr_dp[j].F;
				bt[i][j] = HORIZONTAL;
			}
		}
		swap(prev_dp, curr_dp);
	}
	vector<RepInterval> segments;
	RepInterval tmp;
	tmp.end = que_len;
	int ti = que_len, tj = pat_len;
	while (ti > 0 or tj > 0) {
		switch (bt[ti][tj]) {
			case HORIZONTAL:
				tmp.gap++;
				tj--;
				break;
			case VERTICAL:
				tmp.gap++;
				ti--;
				break;
			case DIAGONAL:
				tmp.mis += (que_seq[ti-1] != pat_seq[tj-1]);
				tmp.mat += (que_seq[ti-1] == pat_seq[tj-1]);
				ti--;
				tj--;
				break;
			case WRAP_1:
				assert(tj == 1);
				tmp.mis += (que_seq[ti-1] != pat_seq[tj-1]);
				tmp.mat += (que_seq[ti-1] == pat_seq[tj-1]);
				tmp.beg = ti - 1;
				ti--;
				tj = pat_len;
				segments.push_back(tmp);
				tmp = RepInterval();
				tmp.end = ti;
				break;
			case WRAP_2:
				assert(tj == 1);
				tmp.gap++;
				tmp.beg = ti - 1;
				tj = pat_len;
				segments.push_back(tmp);
				tmp = RepInterval();
				tmp.end = ti - 1;
				break;
			default:
				break;
		}
	}
	tmp.beg = ti;
	segments.push_back(tmp);
	reverse(segments.begin(), segments.end());
	return segments;
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

struct SgResult {
	int beg, end;
	int score, mis, gap;
	SgResult() {
		beg = end = -1;
		score = -INF;
		mis = gap = 0;
	}
};

// Banding is not applicable
SgResult semi_global(const int n, const char *a, const int m, const char *b)
{
	// Classical SI score matrix
	const int MAT_S = 1;
	const int MIS_P = -4;
	const int GAP_O = -6;
	const int GAP_E = -1;
	const int VERTICAL = 1;
	const int HORIZONTAL = 2;
	const int DIAGONAL = 3;

	vector<int> prev_H(m + 1, -INF), curr_H(m + 1, -INF);
	vector<int> prev_E(m + 1, -INF), curr_E(m + 1, -INF);
	for (int i = 0; i <= m; i++) {
		prev_H[i] = prev_E[i] = 0;
	}
	vector<vector<uint8_t>> bt(n + 1);
	for (int i = 0; i <= n; i++) {
		bt[i].resize(m + 1, 0);
	}

	for (int i = 1; i <= n; i++) {
		curr_H[0] = curr_E[0] = GAP_O + i * GAP_E;
		int F = -INF;
		for (int j = 1; j <= m; j++) {
			F = max(F, curr_H[j-1] + GAP_O) + GAP_E;
			curr_E[j] = max(prev_E[j], prev_H[j] + GAP_O) + GAP_E;
			int M = prev_H[j-1] + (a[i-1] == b[j-1] ?MAT_S :MIS_P);
			if (F > M) {
				curr_H[j] = F;
				bt[i][j] = HORIZONTAL;
			} else {
				curr_H[j] = M;
				bt[i][j] = DIAGONAL;
			}
			if (curr_E[j] > curr_H[j]) {
				curr_H[j] = curr_E[j];
				bt[i][j] = VERTICAL;
			}
		}
		swap(prev_H, curr_H);
		swap(prev_E, curr_E);
	}

	SgResult ret;
	for (int j = 1; j <= m; j++) {
		if (prev_H[j] > ret.score) {
			ret.score = prev_H[j];
			ret.end = j; // 0-based open
		}
	}

	int pi = n, pj = ret.end;
	while (bt[pi][pj] != 0) {
		switch (bt[pi][pj]) {
		case DIAGONAL:
			pi--;
			pj--;
			ret.mis += (a[pi] != b[pj]);
			break;
		case VERTICAL:
			pi--;
			ret.gap++;
			break;
		case HORIZONTAL:
			pj--;
			ret.gap++;
			break;
		default:
			break;
		}
	}
	ret.gap += pi;
	ret.beg = pj; // 0-based closed
	return ret;
}

vector<RepInterval> extend_pattern(const int p_len, const char *p_seq, const int t_len, const char *t_seq)
{
	// Find the start position
	int len = min(t_len, (int)(p_len * 1.5));
	SgResult sg = semi_global(p_len, p_seq, len, t_seq);

	int pos = sg.beg;
	vector<RepInterval> ret = wraparound_dp(p_len,p_seq,t_len - pos, t_seq + pos);
	ret.resize(ret.size() - 1); // Kick off the last one
	for (RepInterval &r: ret) { // Reset coordinates
		// fprintf(stderr, "%d %d %d\n", r.beg, r.end, r.end - r.beg);
		r.beg += pos;
		r.end += pos;
	}

	return ret;
}

RepInterval pickout_pattern(const vector<RepInterval> &reps) {
	const int MIN_LENGTH = 120;
	int sum = 0;
	for (const RepInterval &r : reps) {
		sum += r.end - r.beg;
	}
	double ave_len = 1.0 * sum / reps.size();

	int min_len = INF, max_len = -INF;
	double min_div = INF;
	int min_id = -1;
	for (int i = 0; i < reps.size(); i++) {
		const RepInterval &r = reps[i];
		int len = r.end - r.beg;
		// if (len < MIN_LENGTH) continue;
		double div = abs(len - ave_len);
		if (div < min_div) {
			min_div = div;
			min_id = i;
		}
		min_len = min(min_len, len);
		max_len = max(max_len, len);
	}
	RepInterval ret = reps[min_id];

	fprintf(stderr, "Picked pattern: [%d,%d), len=%d from %ld repeats (%d,%d)\n",
		ret.beg, ret.end, ret.end - ret.beg, reps.size(), min_len, max_len);
	return ret;
}

struct LongRepeats {
	int pid; // Pattern ID
	int base_idx; // Base index
	string pattern;
	vector<RepInterval> repeats;
};

vector<LongRepeats> collect_long_repeats(const ZigOptions &opt, int n, const char *seq)
{
	assert(opt.n_threads > 0);
	int global_os = 0;
	vector<RepInterval> reps_bin[opt.n_threads];
	vector<vector<Bt1Cell>> bt_mat;
	string pattern;
	vector<RepInterval> all_reps;
	vector<LongRepeats> ret;

	// TODO: use thread pipeline to achieve better parallelization
	int batch_id = 0;
	while (global_os < n) {
		double r_start = realtime(), c_start = cputime();
		int m = n - global_os;
		const char *s0 = seq + global_os;
		if (pattern.empty()) {
			// Discover repeat pattern
			int k = min(PART_LEN, m);
			vector<RepInterval> pool = self_alignment2(opt, k, s0, bt_mat); // Multiple pattern?
			vector<RepInterval> reps;
			// Unwind nested repeats recursively to discover as short units as possible
			while (not pool.empty()) {
				vector<RepInterval> next_pool;
				for (const RepInterval &r: pool) {
					if (r.end - r.beg < 1000) {
						reps.push_back(r);
						continue;
					}
					vector<RepInterval> tmp = self_alignment2(opt, r.end - r.beg, s0 + r.beg, bt_mat);
					if (tmp.empty()) reps.push_back(r);
					else {
						for (RepInterval &e: tmp) {
							e.beg += r.beg;
							e.end += r.beg;
							next_pool.push_back(e);
						}
					}
				}
				pool = next_pool;
			}

			if (reps.empty()) {
				// Non-repetitive region (caution for serial execution)
				global_os += k;
				continue;
			}
			RepInterval p = pickout_pattern(reps);
			if (p.beg == -1) {
				global_os += k;
				continue;
			}
			int len = p.end - p.beg;
			pattern.resize(len);
			memcpy((char*)pattern.data(), s0 + p.beg, len);
		}

		// Multi-threaded WDP
		int n_threads = min((m + PART_LEN - 1) / PART_LEN, opt.n_threads);
		#pragma omp parallel for num_threads(n_threads)
		for (int i = 0; i < n_threads; i++) {
			reps_bin[i].clear();
			int os = i * PART_LEN;
			int t_len = min(m - os, PART_LEN);
			// fprintf(stderr, "thread %d, len = %d\n", i + 1, t_len);
			const char *s = s0 + os;
			reps_bin[i] = extend_pattern(pattern.length(), pattern.data(), t_len, s);
			for (RepInterval &r: reps_bin[i]) { // Add local offset
				r.beg += os;
				r.end += os;
			}
		}
		for (int i = 0; i < n_threads; i++) {
			// fprintf(stderr, "[%d, %d)\n", reps_bin[i].front().beg, reps_bin[i].back().end);
			for (RepInterval &r: reps_bin[i]) { // Add global offset
				r.beg += global_os;
				r.end += global_os;
			}
		}

		// Single-thread Stitching
		vector<RepInterval> sti_reps;
		sti_reps.insert(sti_reps.end(), reps_bin[0].begin(), reps_bin[0].end());
		for (int i = 1; i < n_threads; i++) {
			if (reps_bin[i-1].size() > 0 and reps_bin[i].size() > 0) {
				// Fill the gap in between
				int beg = reps_bin[i-1].back().end;
				int end = reps_bin[i].front().beg;
				vector<RepInterval> mid = wraparound_dp(
					pattern.length(), pattern.data(), end - beg, s0 + beg);
				for (RepInterval&r : mid) {
					r.beg += beg;
					r.end += beg;
					sti_reps.push_back(r);
				}
			}
			sti_reps.insert(sti_reps.end(), reps_bin[i].begin(), reps_bin[i].end());
		}

		// Sanity check: partitioned repeats must be sorted
		for (int i = 0; i < sti_reps.size(); i++) {
			assert(sti_reps[i].end > sti_reps[i].beg);
			if (i > 0) {
				assert(sti_reps[i].beg >= sti_reps[i-1].end);
			}
		}
		if (sti_reps.empty()) break; // FIXME
		assert(sti_reps.size() > 0); // At least the pattern itself will be identified

		// Poorly aligned regions
		const int MAX_DIV_CNT = 10; // Too lenient?
		vector<pair<int,int>> poor_range;
		int last_p = -1;
		for (int i = 0; i < sti_reps.size(); i++) {
			const RepInterval &r = sti_reps[i];
			double match_ratio = 1.0 -  (double)(r.mis + r.gap) / pattern.length();
			if (match_ratio < MIN_MATCH_RATIO) {
				// fprintf(stderr, "%d [%d, %d) mat=%d, mis=%d, gap=%d\n", i, r.beg, r.end, r.mat, r.mis, r.gap);
				if (last_p == -1) {
					last_p = i;
				}
			} else {
				if (last_p != -1 and i - last_p > MAX_DIV_CNT) {
					poor_range.push_back(make_pair(last_p, i));
				}
				last_p = -1;
			}
		}
		if (last_p != -1 and sti_reps.size() - last_p > MAX_DIV_CNT) {
			poor_range.push_back(make_pair(last_p, sti_reps.size()));
		}
		for (const pair<int,int> &pair: poor_range) {
			// fprintf(stderr, "[%d, %d)\n", pair.first, pair.second);
		}

		int old_n = all_reps.size();
		last_p = 0;
		for (const pair<int,int> &pair: poor_range) {
			for (int i = last_p; i < pair.first; i++) {
				all_reps.push_back(sti_reps[i]);
			}
			last_p = pair.second;
		}
		for (int i = last_p; i < sti_reps.size(); i++) {
			all_reps.push_back(sti_reps[i]);
		}

		int added_reps = all_reps.size() - old_n;
		// assert(added_reps > 0);
		fprintf(stderr, "Batch %d, offset: %d, added_reps: %d, real_time: %.2f, CPU_time: %.2f\n",
			++batch_id, global_os, added_reps, realtime() - r_start, cputime() - c_start);

		// End position of the current repeat pattern
		global_os = all_reps.back().end;
		if (not poor_range.empty()) {
			LongRepeats lr;
			lr.pattern = pattern;
			lr.repeats = all_reps;
			ret.push_back(lr);
			pattern.clear();
			all_reps.clear();
		}

		// TODO: post-process the poorly aligned regions
	}
	if (not pattern.empty()) {
		LongRepeats lr;
		lr.pattern = pattern;
		lr.repeats = all_reps;
		ret.push_back(lr);
	}

	fprintf(stderr, "Found %ld patterns\n", ret.size());
	for (int i = 0; i < ret.size(); i++) {
		const LongRepeats &r = ret[i];
		const string &s = r.pattern;
		int beg = r.repeats.front().beg;
		int end = r.repeats.back().end;
		int len = end - beg;
		fprintf(stderr, "Pattern %d: len=%ld, range=[%d,%d), len=%d\n", i+1, s.length(), beg, end, len);
	}
	return ret;
}

LongRepeats normalize_repeats(const ZigOptions &opt, int pat_len, const char *pat_seq, int n, const char *seq)
{
	assert(opt.n_threads > 0);
	int n_threads = opt.n_threads;
	int global_os = 0;
	vector<RepInterval> reps_bin[n_threads];
	string rep_pat;
	rep_pat.resize(pat_len);
	memcpy((char*)rep_pat.data(), pat_seq, pat_len);
	vector<RepInterval> all_reps;

	int batch_id = 0;
	while (global_os < n) {
		double r_start = realtime(), c_start = cputime();
		int m = n - global_os;
		const char *s0 = seq + global_os;

		// Multi-threading DP
		#pragma omp parallel for
		for (int i = 0; i < n_threads; i++) {
			reps_bin[i].clear();
			int os = i * PART_LEN;
			int t_len = min(m - os, PART_LEN);
			if (t_len <= 0) continue;
			const char *s = s0 + os;
			reps_bin[i] = extend_pattern(pat_len, pat_seq, t_len, s);
			for (RepInterval &r: reps_bin[i]) {
				r.beg += os;
				r.end += os;
			}
		}

		// Single-thread Stitching
		int old_n = all_reps.size();
		for (int i = 0; i < n_threads; i++) {
			for (RepInterval &r: reps_bin[i]) {
				r.beg += global_os;
				r.end += global_os;
				all_reps.push_back(r);
			}
			if (i > 0) {
				if (reps_bin[i-1].empty() or reps_bin[i].empty()) continue;
				// Fill the gap in between
				int beg = reps_bin[i-1].back().end;
				int end = reps_bin[i].front().beg;
				vector<RepInterval> mid = wraparound_dp(pat_len, pat_seq, end - beg, s0 + beg);
				for (RepInterval&r : mid) {
					r.beg += beg; // Untidy style
					r.end += beg;
					all_reps.push_back(r);
				}
			}
		}
		sort(all_reps.begin(), all_reps.end());

		int added_reps = all_reps.size() - old_n;
		if (added_reps == 0) break;
		// fprintf(stderr, "Batch %d, offset: %d, added_reps: %d, real_time: %.2f, CPU_time: %.2f\n",
		// 	++batch_id, global_os, added_reps, realtime() - r_start, cputime() - c_start);

		// Find the end position of the repeat pattern
		global_os = all_reps.back().end;
	}
	LongRepeats ret;
	ret.pattern = rep_pat;
	ret.repeats = all_reps;
	return ret;
}

// ------------------------ //

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

	const int NORMAL = 0;
	const int START_REP = 1;
	const int NEW_COPY = 2;
	const int WITHIN_REP = 3;
	const int END_REP = 4;

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

	if (opt.log_prefix) {
		string pair_vis_fn = string(opt.log_prefix) + "_p.txt";
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

void extended_paf_format(const string &t_name, int t_len, const string &ext_t, const string &q_name, int q_len, const string ext_q)
{
	assert(ext_t.length() == ext_q.length());
	char strand = '+';
	int matches_n = 0, mismatches_n = 0, len = ext_t.length();
	int mapq = 60;

	string cigar;
	char op_type = '@';
	int op_cnt = 0;
	for (int i = 0; i < len; i++) {
		if (ext_t[i] == '[' or ext_t[i] == ']') {
			if (op_type != '@') {
				cigar += to_string(op_cnt);
				cigar += op_type;
			}
			cigar += ext_t[i];
			op_type = '@';
			op_cnt = 0;
		} else {
			int curr_ot;
			if (ext_t[i] == '-') {
				curr_ot = 'I';
			} else if (ext_q[i] == '-') {
				curr_ot = 'D';
			} else {
				curr_ot = 'M'; // Match and mismatches
				if (ext_t[i] == ext_q[i]) matches_n++;
				else mismatches_n++;
			}

			if (curr_ot == op_type) op_cnt++;
			else {
				if (op_type != '@') {
					cigar += to_string(op_cnt);
					cigar += op_type;
				}
				op_type = curr_ot;
				op_cnt = 1;
			}
		}
	}
	if (op_type != '@' and op_cnt > 0) {
		cigar += to_string(op_cnt);
		cigar += op_type;
	}

	fprintf(stdout, "%s\t%d\t%d\t%d\t%c\t", q_name.c_str(), q_len, 0, q_len, strand);
	fprintf(stdout, "%s\t%d\t%d\t%d\t", t_name.c_str(), t_len, 0, t_len);
	fprintf(stdout, "%d\t%d\t%d\t", matches_n + mismatches_n, len, mapq);
	fprintf(stdout, "cg:Z:%s\n", cigar.c_str());
}

void align_with_dups(const ZigOptions &opt, const char *fn1, const char *fn2) {
	pair<string, string> pair1 = input_fasta_seq(fn1);
	pair<string, string> pair2 = input_fasta_seq(fn2);
	string name1 = pair1.first, seq1 = pair1.second;
	string name2 = pair2.first, seq2 = pair2.second;
	int t_len = seq1.length(), q_len = seq2.length();
	const char *t = seq1.data(), *q = seq2.data();
	if (t_len > SA_MAX_LEN or q_len > SA_MAX_LEN) {
		vector<int> t_bp = (t_len > SA_MAX_LEN) ?process_long(opt, t_len, t) :self_alignment(opt, t_len, t, "");
		vector<int> q_bp = (q_len > SA_MAX_LEN) ?process_long(opt, q_len, q) :self_alignment(opt, q_len, q, "");
		// It is slow because the partial matrix is still large
		vector<int> cv2 = banded_pairwise(opt, t_len, t, t_bp, q_len, q, q_bp);
		paf_format(name2, seq2, name1, seq1, cv2);
	} else {
		// banded_pairwise(opt, t_len, t, q_len, q);
		string t_vis_fn = "", q_vis_fn = "";
		if (opt.log_prefix) {
			t_vis_fn = string(opt.log_prefix) + "_s1.txt";
			q_vis_fn = string(opt.log_prefix) + "_s2.txt";
		}
		vector<int> t_bp = self_alignment(opt, t_len, t, t_vis_fn);
		vector<int> q_bp = self_alignment(opt, q_len, q, q_vis_fn);
		vector<int> cv2 = banded_pairwise(opt, t_len, t, t_bp, q_len, q, q_bp);
		paf_format(name2, seq2, name1, seq1, cv2);
	}
}

// ------------------------ //

struct LocalResult {
	int max_score;
	int beg_a, end_a;
	int beg_b, end_b;
};

// Use SIMD or other methods to accelerate it (banded alignment generates poor results)
LocalResult local_alignment(const int n, const char *a, const int m, const char *b)
{
	// Classical SI score matrix
	const int MAT_S = 1;
	const int MIS_P = -4;
	const int GAP_O = -6;
	const int GAP_E = -1;

	const int STOP = 0;
	const int VERTICAL = 1;
	const int HORIZONTAL = 2;
	const int DIAGONAL = 3;

	vector<int> prev_H(m + 1, 0), curr_H(m + 1, 0);
	vector<int> prev_E(m + 1, 0), curr_E(m + 1, 0);
	vector<vector<uint8_t>> bt(n + 1);
	for (int i = 0; i <= n; i++) bt[i].resize(m + 1, STOP);
	int max_score = -1;
	int end_i = 0, end_j = 0;

	for (int i = 1; i <= n; i++) {
		int F = 0;
		memset(curr_H.data(), 0, (m + 1) * sizeof(int));
		for (int j = 1; j <= m; j++) {
			F = max(F, curr_H[j-1] + GAP_O) + GAP_E;
			curr_E[j] = max(prev_E[j], prev_H[j] + GAP_O) + GAP_E;
			int M = prev_H[j-1] + (a[i-1] == b[j-1] ?MAT_S :MIS_P);
			if (F > curr_H[j]) {
				curr_H[j] = F;
				bt[i][j] = HORIZONTAL;
			}
			if (curr_E[j] > curr_H[j]) {
				curr_H[j] = curr_E[j];
				bt[i][j] = VERTICAL;
			}
			if (M > curr_H[j]) {
				curr_H[j] = M;
				bt[i][j] = DIAGONAL;
			}
			if (curr_H[j] > max_score) {
				max_score = curr_H[j];
				end_i = i;
				end_j = j;
			}
		}
		swap(prev_H, curr_H);
		swap(prev_E, curr_E);
	}

	int pi = end_i, pj = end_j;
	string ext_a, ext_b, align;
	while (bt[pi][pj] != STOP) {
		switch (bt[pi][pj]) {
		case DIAGONAL:
			ext_a += a[--pi];
			ext_b += b[--pj];
			align += (a[pi] == b[pj] ?' ' :'X');
			break;
		case VERTICAL:
			ext_a += a[--pi];
			ext_b += '-';
			align += ' ';
			break;
		case HORIZONTAL:
			ext_a += ' ';
			ext_b += b[--pj];
			align += ' ';
			break;
		default:
			break;
		}
	}
	reverse(ext_a.begin(), ext_a.end());
	reverse(ext_b.begin(), ext_b.end());
	reverse(align.begin(), align.end());

	// cout << ext_a << endl;
	// cout << align << endl;
	// cout << ext_b << endl;
	// printf("[%d, %d) aligns with [%d, %d), max_score=%d\n", pi, end_i, pj, end_j, max_score);

	LocalResult ret;
	ret.max_score = max_score;
	ret.beg_a = pi;
	ret.end_a = end_i;
	ret.beg_b = pj;
	ret.end_b = end_j;
	return ret;
}

int global_alignment(const ZigOptions &opt, const int n, const char *a, const int m, const char *b)
{
	const int UNIT_MAT_S = opt.mat_score;
	const int UNIT_MIS_P = opt.mis_pen;
	const int UNIT_GAP_O = opt.gap_o;
	const int UNIT_GAP_E = opt.gap_e;
	const double GAP_RATIO = 0.10;
	const int w = max(n, m) * GAP_RATIO;

	vector<int> prev_H(m + 1, -INF), curr_H(m + 1, -INF);
	vector<int> prev_E(m + 1, -INF), curr_E(m + 1, -INF);
	prev_H[0] = 0;
	for (int i = 1; i <= m; i++) {
		prev_H[i] = prev_E[i] = UNIT_GAP_O + i * UNIT_GAP_E;
	}

	int ret = -INF;
	for (int i = 1; i <= n; i++) {
		int beg = max(i - w, 1), end = min(i + w, m);
		if (beg > end) break;
		if (beg == 1) curr_H[beg-1] = UNIT_GAP_O + i * UNIT_GAP_E;
		else curr_H[beg-1] = -INF;
		int F = -INF;
		int max_score = -INF;
		for (int j = beg; j <= end; j++) {
			F = max(F, curr_H[j-1] + UNIT_GAP_O) + UNIT_GAP_E;
			curr_E[j] = max(prev_E[j], prev_H[j] + UNIT_GAP_O) + UNIT_GAP_E;
			int M = prev_H[j-1] + (a[i-1] == b[j-1] ?UNIT_MAT_S :UNIT_MIS_P);
			curr_H[j] = max(F, M);
			curr_H[j] = max(curr_H[j], curr_E[j]);
			max_score = max(max_score, curr_H[j]);
		}
		if (end < m) curr_H[end+1] = curr_E[end+1] = -INF;
		if (max_score == -INF) break;
		if (i == n and end == m) ret = curr_H[m];

		swap(prev_H, curr_H);
		swap(prev_E, curr_E);
	}
	// if (ret < 0) ret = -INF;
	return ret;
}

struct AlnSta {
	int match, mismatch, ins, del;
	string ext_a, ext_b;
	AlnSta() {
		match = mismatch = ins = del = 0;
	}
	void operator += (const AlnSta &b) {
		this->match += b.match;
		this->mismatch += b.mismatch;
		this->ins += b.ins;
		this->del += b.del;
	}
};

// Banded alignment is not used
AlnSta global_cigar(const ZigOptions &opt, const int n, const char *a, const int m, const char *b)
{
	const int UNIT_MAT_S = opt.mat_score;
	const int UNIT_MIS_P = opt.mis_pen;
	const int UNIT_GAP_O = opt.gap_o;
	const int UNIT_GAP_E = opt.gap_e;
	const int VERTICAL = 1;
	const int HORIZONTAL = 2;
	const int DIAGONAL = 3;

	vector<int> prev_H(m + 1, -INF), curr_H(m + 1, -INF);
	vector<int> prev_E(m + 1, -INF), curr_E(m + 1, -INF);
	prev_H[0] = 0;
	for (int i = 1; i <= m; i++) {
		prev_H[i] = prev_E[i] = UNIT_GAP_O + i * UNIT_GAP_E;
	}
	vector<vector<uint8_t>> bt(n + 1);
	for (int i = 0; i <= n; i++) {
		bt[i].resize(m + 1, 0);
	}

	for (int i = 1; i <= n; i++) {
		int beg = 1, end = m;
		if (beg == 1) curr_H[beg-1] = UNIT_GAP_O + i * UNIT_GAP_E;
		else curr_H[beg-1] = -INF;
		int F = -INF;
		int max_score = -INF;
		for (int j = beg; j <= end; j++) {
			F = max(F, curr_H[j-1] + UNIT_GAP_O) + UNIT_GAP_E;
			curr_E[j] = max(prev_E[j], prev_H[j] + UNIT_GAP_O) + UNIT_GAP_E;
			int M = prev_H[j-1] + (a[i-1] == b[j-1] ?UNIT_MAT_S :UNIT_MIS_P);
			if (F > M) {
				curr_H[j] = F;
				bt[i][j] = HORIZONTAL;
			} else {
				curr_H[j] = M;
				bt[i][j] = DIAGONAL;
			}
			if (curr_E[j] > curr_H[j]) {
				curr_H[j] = curr_E[j];
				bt[i][j] = VERTICAL;
			}
			max_score = max(max_score, curr_H[j]);
		}
		if (end < m) curr_H[end+1] = curr_E[end+1] = -INF;
		if (max_score == -INF) break;

		swap(prev_H, curr_H);
		swap(prev_E, curr_E);
	}

	assert(prev_H[m] != -INF);

	AlnSta ret;
	int pi = n, pj = m;
	while (bt[pi][pj] != 0) {
		switch (bt[pi][pj]) {
		case DIAGONAL:
			pi--;
			pj--;
			if (a[pi] == b[pj]) ret.match++;
			else ret.mismatch++;
			ret.ext_a += a[pi];
			ret.ext_b += b[pj];
			break;
		case VERTICAL:
			pi--;
			ret.del++;
			ret.ext_a += a[pi];
			ret.ext_b += '-';
			break;
		case HORIZONTAL:
			pj--;
			ret.ins++;
			ret.ext_a += '-';
			ret.ext_b += b[pj];
			break;
		default:
			break;
		}
	}
	// Be careful with the preceding gaps, which are introduced by wrong splitting
	ret.ins += pi;
	while (pi > 0) {
		ret.ext_a += a[--pi];
		ret.ext_b += '-';
	}
	ret.del += pj;
	while (pj > 0) {
		ret.ext_a += '-';
		ret.ext_b += b[--pj];
	}
	reverse(ret.ext_a.begin(), ret.ext_a.end());
	reverse(ret.ext_b.begin(), ret.ext_b.end());

	// Sanity check
	{
		assert(ret.ext_a.length() == ret.ext_b.length());
		int i = 0;
		for (char c: ret.ext_a) {
			if (c != '-') {
				assert(c == a[i++]);
			}
		}
		assert(i == n);

		i = 0;
		for (char c: ret.ext_b) {
			if (c != '-') {
				assert(c == b[i++]);
			}
		}
		assert(i == m);
	}
	return ret;
}

void align_long_seq(const ZigOptions &opt, const char *fn1, const char *fn2)
{
	pair<string, string> pair1 = input_fasta_seq(fn1);
	pair<string, string> pair2 = input_fasta_seq(fn2);
	string name1 = pair1.first, seq1 = pair1.second;
	string name2 = pair2.first, seq2 = pair2.second;
	int t_len = seq1.length(), q_len = seq2.length();
	const char *t_seq = seq1.data(), *q_seq = seq2.data();

	vector<LongRepeats> lr_t = collect_long_repeats(opt, t_len, t_seq);
	vector<LongRepeats> lr_q = collect_long_repeats(opt, q_len, q_seq);

	if (opt.log_prefix) {
		FILE* fo = fopen((string(opt.log_prefix) + "_t_pattern.tsv").c_str(), "w");
		assert(fo);
		fprintf(fo, "%s\t%s\t%s\n", "ID", "Length", "Pattern");
		for (int i = 0; i < lr_t.size(); i++) {
			fprintf(fo, "%d\t%ld\t%s\n", i + 1, lr_t[i].pattern.length(), lr_t[i].pattern.data());
		}
		fclose(fo);

		fo = fopen((string(opt.log_prefix) + "_t_reps.tsv").c_str(), "w");
		assert(fo);
		fprintf(fo, "%s\t%s\t%s\t%s\t%s\t%s\n", "ID", "beg", "end", "len", "mis", "gap");
		for (int i = 0; i < lr_t.size(); i++) {
			for (const RepInterval &r: lr_t[i].repeats) {
				fprintf(fo, "%d\t%d\t%d\t%d\t%d\t%d\n", i + 1, r.beg, r.end, r.end - r.beg, r.mis, r.gap);
			}
		}
		fclose(fo);

		fo = fopen((string(opt.log_prefix) + "_q_pattern.tsv").c_str(), "w");
		assert(fo);
		fprintf(fo, "%s\t%s\t%s\n", "ID", "Length", "Pattern");
		for (int i = 0; i < lr_q.size(); i++) {
			fprintf(fo, "%d\t%ld\t%s\n", i + 1, lr_q[i].pattern.length(), lr_q[i].pattern.data());
		}
		fclose(fo);

		fo = fopen((string(opt.log_prefix) + "_q_reps.tsv").c_str(), "w");
		assert(fo);
		fprintf(fo, "%s\t%s\t%s\t%s\t%s\t%s\n", "ID", "beg", "end", "len", "mis", "gap");
		for (int i = 0; i < lr_q.size(); i++) {
			for (const RepInterval &r: lr_q[i].repeats) {
				fprintf(fo, "%d\t%d\t%d\t%d\t%d\t%d\n", i + 1, r.beg, r.end, r.end - r.beg, r.mis, r.gap);
			}
		}
		fclose(fo);
	}
	int t_sum = 0, q_sum = 0;
	for (const LongRepeats &lr: lr_t) {
		for (const RepInterval &r: lr.repeats) {
			t_sum += r.end - r.beg;
		}
	}
	for (const LongRepeats &lr: lr_q) {
		for (const RepInterval &r: lr.repeats) {
			q_sum += r.end - r.beg;
		}
	}
	fprintf(stderr, "Target repeat fraction[%%]: %.2f\n", 100.0 * t_sum / t_len);
	fprintf(stderr, "Query repeat fraction[%%]: %.2f\n", 100.0 * q_sum / q_len);

	// Compare patterns pairwise
	int sum_pattern = lr_q.size() + lr_t.size();
	vector<int> parent_set(sum_pattern);
	for (int i = 0; i < sum_pattern; i++) {
		parent_set[i] = i;
	}
	const int MAX_PATTERN_DIS = 500000;
	for (int i = 0; i < lr_t.size(); i++) {
		const string &a = lr_t[i].pattern;
		int max_value = -INF, max_id = -1;
		int t_beg = lr_t[i].repeats.front().beg;
		int t_end = lr_t[i].repeats.back().end;
		for (int j = 0; j < lr_q.size(); j++) {
			const string &b = lr_q[j].pattern;
			int q_beg = lr_q[j].repeats.front().beg;
			int q_end = lr_q[j].repeats.back().end;
			int dis = max(t_beg, q_beg) - min(t_end, q_end);
			if (dis > MAX_PATTERN_DIS) {
				continue;
			}

			// fprintf(stderr, "%d(%ld) -> %d(%ld)\t", i, a.length(), j, b.length());
			string c = b + b;
			SgResult sg = semi_global(a.length(), a.data(), c.length(), c.data());
			double mat_ratio = 1.0 - 1.0 * (sg.mis + sg.gap) / a.length();
			// fprintf(stderr, "[%d,%d) sim=%.2f\n", sg.beg, sg.end, mat_ratio);
			if (mat_ratio < MIN_MATCH_RATIO) continue;
			if (sg.score > max_value) {
				max_value = sg.score;
				max_id = j;
			}
		}
		if (max_id != -1) {
			// fprintf(stderr, "Union %d %d\n", i, max_id);
			us_union(parent_set, i, lr_t.size() + max_id);
		}
	}
	for (int i = 0; i < lr_q.size(); i++) {
		const string &a = lr_q[i].pattern;
		int max_value = -INF, max_id = -1;
		int q_beg = lr_q[i].repeats.front().beg;
		int q_end = lr_q[i].repeats.back().end;
		for (int j = 0; j < lr_t.size(); j++) {
			const string &b = lr_t[j].pattern;
			int t_beg = lr_t[j].repeats.front().beg;
			int t_end = lr_t[j].repeats.back().end;
			int dis = max(t_beg, q_beg) - min(t_end, q_end);
			if (dis > MAX_PATTERN_DIS) {
				continue;
			}

			// fprintf(stderr, "%d(%ld) -> %d(%ld)\t", i, a.length(), j, b.length());
			string c = b + b;
			SgResult sg = semi_global(a.length(), a.data(), c.length(), c.data());
			double mat_ratio = 1.0 - 1.0 * (sg.mis + sg.gap) / a.length();
			// fprintf(stderr, "[%d,%d) sim=%.2f\n", sg.beg, sg.end, mat_ratio);
			if (mat_ratio < MIN_MATCH_RATIO) continue;
			if (sg.score > max_value) {
				max_value = sg.score;
				max_id = j;
			}
		}
		if (max_id != -1) {
			// fprintf(stderr, "Union %d %d\n", max_id, i);
			us_union(parent_set, i + lr_t.size(), max_id);
		}
	}
	vector<int> same_pat[sum_pattern];
	for (int i = 0; i < sum_pattern; i++) {
		int k = us_find(parent_set, i);
		// fprintf(stderr, "parent[%d] = %d\n", i, k);
		same_pat[k].push_back(i);
		if (i < lr_t.size()) lr_t[i].pid = k;
		else lr_q[i - lr_t.size()].pid = k;
	}

	// Normalize repeats
	int count = 0;
	vector<bool> kept(sum_pattern, false);
	for (const vector<int> &x: same_pat) {
		if (x.empty()) continue;
		// fprintf(stderr, "cluster %d:\n", ++count);
		int max_range = 0, max_id = -1;
		for (int i: x) {
			int len = 0;
			if (i < lr_t.size()) {
				const string &p = lr_t[i].pattern;
				const vector<RepInterval> &r = lr_t[i].repeats;
				len = r.back().end - r.front().beg;
				// fprintf(stderr, "i=%d, len=%ld, pattern:%s\n", i, p.length(), p.data());
			} else {
				int j = i - lr_t.size();
				const string &p = lr_q[j].pattern;
				const vector<RepInterval> &r = lr_q[j].repeats;
				len = r.back().end - r.front().beg;
				// fprintf(stderr, "j=%d, len=%ld, pattern:%s\n", j, p.length(), p.data());
			}
			if (len > max_range) {
				max_range = len;
				max_id = i;
			}
		}
		if (max_id == -1) continue;
		kept[max_id] = true;

		// Use the pattern with the longest repeat length to normalize other repeats
		string p;
		if (max_id < lr_t.size()) p = lr_t[max_id].pattern;
		else p = lr_q[max_id - lr_t.size()].pattern;
		int pid = us_find(parent_set, max_id);
		for (int i: x) {
			if (i == max_id) continue;
			LongRepeats lr;
			if (i < lr_t.size()) {
				vector<RepInterval> &r = lr_t[i].repeats;
				int os = r.front().beg;
				int len = r.back().end - r.front().beg;
				const char *s = t_seq + os;
				lr = normalize_repeats(opt, p.length(), p.data(), len, s);
				lr_t[i].pid = pid;
				lr_t[i].pattern = p;
				lr_t[i].repeats = lr.repeats;
				for (RepInterval &t: lr_t[i].repeats) {
					t.beg += os;
					t.end += os;
				}
			} else {
				int j = i - lr_t.size();
				const vector<RepInterval> &r = lr_q[j].repeats;
				int os = r.front().beg;
				int len = r.back().end - r.front().beg;
				const char *s = q_seq + os;
				lr = normalize_repeats(opt, p.length(), p.data(), len, s);
				lr_q[j].pid = pid;
				lr_q[j].pattern = p;
				lr_q[j].repeats = lr.repeats;
				for (RepInterval &t: lr_q[j].repeats) {
					t.beg += os;
					t.end += os;
				}
			}
		}
	}

	if (opt.log_prefix) {
		fprintf(stderr, "Outputting updated patterns and repeats\n");
		FILE* fo = fopen((string(opt.log_prefix) + "_all_pattern.tsv").c_str(), "w");
		assert(fo);
		fprintf(fo, "%s\t%s\t%s\n", "ID", "Length", "Pattern");
		for (int i = 0; i < sum_pattern; i++) {
			if (not kept[i]) continue;
			int pid = us_find(parent_set, i);
			if (i < lr_t.size()) {
				fprintf(fo, "%d\t%ld\t%s\n", pid, lr_t[i].pattern.length(), lr_t[i].pattern.data());
			} else {
				int j = i - lr_t.size();
				fprintf(fo, "%d\t%ld\t%s\n", pid, lr_t[j].pattern.length(), lr_t[j].pattern.data());
			}
		}
		fclose(fo);

		fo = fopen((string(opt.log_prefix) + "_t_reps.tsv").c_str(), "w");
		assert(fo);
		fprintf(fo, "%s\t%s\t%s\t%s\t%s\t%s\n", "ID", "beg", "end", "len", "mis", "gap");
		for (int i = 0; i < lr_t.size(); i++) {
			int pid = lr_t[i].pid;
			for (const RepInterval &r: lr_t[i].repeats) {
				fprintf(fo, "%d\t%d\t%d\t%d\t%d\t%d\n", pid, r.beg, r.end, r.end - r.beg, r.mis, r.gap);
			}
		}
		fclose(fo);

		fo = fopen((string(opt.log_prefix) + "_q_reps.tsv").c_str(), "w");
		assert(fo);
		fprintf(fo, "%s\t%s\t%s\t%s\t%s\t%s\n", "ID", "beg", "end", "len", "mis", "gap");
		for (int i = 0; i < lr_q.size(); i++) {
			int pid = lr_q[i].pid;
			for (const RepInterval &r: lr_q[i].repeats) {
				fprintf(fo, "%d\t%d\t%d\t%d\t%d\t%d\n", pid, r.beg, r.end, r.end - r.beg, r.mis, r.gap);
			}
		}
		fclose(fo);
	}

	// TODO: how to process singleton patterns?

	int n = 0, m = 0;
	t_sum = 0; q_sum = 0;
	for (LongRepeats &lr: lr_t) {
		lr.base_idx = n;
		n += lr.repeats.size();
		for (const RepInterval &r: lr.repeats) {
			t_sum += r.end - r.beg;
		}
	}
	for (LongRepeats &lr: lr_q) {
		lr.base_idx = m;
		m += lr.repeats.size();
		for (const RepInterval &r: lr.repeats) {
			q_sum += r.end - r.beg;
		}
	}
	fprintf(stderr, "T: %d repeats, fraction=%.2f %%\n", n, 100.0 * t_sum / t_len);
	fprintf(stderr, "Q: %d repeats, fraction=%.2f %%\n", m, 100.0 * q_sum / q_len);

	// Build score matrix of repeat units
	double t_cpu = cputime(), t_real = realtime();
	vector<vector<int>> matrix(n); // FIXME: memory inefficient for this sparse matrix
	for (int i = 0; i < n; i++) matrix[i].resize(m, -INF);
	for (const LongRepeats &lt : lr_t) {
		int idx_t = lt.base_idx;
		int pid_t = lt.pid;
		for (const LongRepeats &lq : lr_q) {
			int idx_q = lq.base_idx;
			int pid_q = lq.pid;
			if (pid_q != pid_t) continue;

			int os_t = lt.repeats.front().beg;
			int os_q = lq.repeats.front().beg;
			#pragma omp parallel for
			for (int i = 0; i < lt.repeats.size(); i++) {
				const RepInterval &t = lt.repeats[i];
				for (int j = 0; j < lq.repeats.size(); j++) {
					const RepInterval &q = lq.repeats[j];
					int dis = abs((t.beg - os_t) - (q.beg - os_q)); // Ignore the global offset
					if (dis < MAX_UNIT_DIS) {
						matrix[idx_t + i][idx_q + j] = global_alignment(
							opt, t.end - t.beg, t_seq + t.beg, q.end - q.beg, q_seq + q.beg);
					}
				}
			}
		}
	}
	fprintf(stderr, "Build scoring matrix: %.2f real time, %.2f CPU time\n", realtime() - t_real, cputime() - t_cpu);


	// for (int i = 0; i < n; i++) {
	// 	for (int j = 0; j < m; j++) {
	// 		if (matrix[i][j] != -INF) {
	// 			fprintf(stdout, "%d -> %d: %d\n", i, j, matrix[i][j]);
	// 		}
	// 	}
	// }

	// TODO: set the penalty of deleting units
	// sqrt(2 * L * max_div * mis_penalty)
	t_cpu = cputime(); t_real = realtime();
	const int DEL_UNIT = opt.del_unit;
	const int VERTICAL = 1;
	const int HORIZONTAL = 2;
	const int DIAGONAL = 3;
	vector<vector<int>> dp(n + 1);
	vector<vector<uint8_t>> bt(n + 1);
	for (int i = 0; i <= n; i++) {
		dp[i].resize(m + 1, -INF);
		bt[i].resize(m + 1, 0);
	}
	dp[0][0] = 0;
	for (int j = 1; j <= m; j++) {
		dp[0][j] = DEL_UNIT * j;
	}
	for (int i = 1; i <= n; i++) {
		dp[i][0] = DEL_UNIT * i;
		for (int j = 1; j <= m; j++) {
			int v = dp[i-1][j] + DEL_UNIT;
			int h = dp[i][j-1] + DEL_UNIT;
			int d = dp[i-1][j-1] + matrix[i-1][j-1];
			if (v > dp[i][j]) {
				dp[i][j] = v;
				bt[i][j] = VERTICAL;
			}
			if (h > dp[i][j]) {
				dp[i][j] = h;
				bt[i][j] = HORIZONTAL;
			}
			if (d > dp[i][j]) {
				dp[i][j] = d;
				bt[i][j] = DIAGONAL;
			}
		}
	}
	fprintf(stderr, "DP: %.2f real time, %.2f CPU time\n", realtime() - t_real, cputime() - t_cpu);

	// Construct alignment topology
	t_cpu = cputime(); t_real = realtime();
	vector<pair<int,int>> aln;
	vector<int> t_del;
	vector<int> q_del;
	int pi = n, pj = m;
	while (bt[pi][pj] != 0) {
		switch (bt[pi][pj]) {
		case DIAGONAL:
			pi--;
			pj--;
			aln.emplace_back(make_pair(pi, pj));
			break;
		case VERTICAL:
			pi--;
			t_del.push_back(pi);
			break;
		case HORIZONTAL:
			pj--;
			q_del.push_back(pj);
			break;
		default:
			break;
		}
	}
	reverse(aln.begin(), aln.end());
	reverse(t_del.begin(), t_del.end());
	reverse(q_del.begin(), q_del.end());

	int aln_len = 0, t_del_len = 0, q_del_len = 0;
	FILE *fo = opt.log_prefix ? fopen((string(opt.log_prefix) + "_aln.tsv").c_str(), "w") :nullptr;
	if (fo) {
		fprintf(fo, "%s\t%s\t%s\t%s\t", "t_id", "t_beg", "t_end", "t_len");
		fprintf(fo, "%s\t%s\t%s\t%s\t", "q_id", "q_beg", "q_end", "q_len");
		fprintf(fo, "%s\t%s\t%s\n", "match", "mismatch", "gap");
	}
	vector<RepInterval> agg_t;
	for (const LongRepeats &lr: lr_t) {
		agg_t.insert(agg_t.end(), lr.repeats.begin(), lr.repeats.end());
	}
	vector<RepInterval> agg_q;
	for (const LongRepeats &lr: lr_q) {
		agg_q.insert(agg_q.end(), lr.repeats.begin(), lr.repeats.end());
	}

	AlnSta mut; // Small mutations
	for (const pair<int,int> &p: aln) {
		int i = p.first, j = p.second;
		const RepInterval &t = agg_t[i];
		const RepInterval &q = agg_q[j];
		aln_len += min(t.end - t.beg, q.end - q.beg);
		AlnSta res = global_cigar(opt, t.end - t.beg, t_seq + t.beg, q.end - q.beg, q_seq + q.beg);
		mut += res;
		if (fo) {
			fprintf(fo, "%d\t%d\t%d\t%d\t", i, t.beg, t.end, t.end - t.beg);
			fprintf(fo, "%d\t%d\t%d\t%d\t", j, q.beg, q.end, q.end - q.beg);
			fprintf(fo, "%d\t%d\t%d\n", res.match, res.mismatch, res.ins + res.del);
		}
	}
	if (fo) fclose(fo);

	double mut_rate = 100.0 * (mut.mismatch + mut.ins + mut.del) / aln_len;
	fprintf(stderr, "Alignment between repeats:\n");
	fprintf(stderr, "    aligned length: %d, map ratio: %.2f %%\n", aln_len, 100.0 * aln_len / min(t_len, q_len));
	fprintf(stderr, "    matches: %d, mismatches: %d, insertions: %d, deletions: %d, mutation rate: %.2f %%\n",
		mut.match, mut.mismatch, mut.ins, mut.del, mut_rate);

	// Target deletion
	fo = opt.log_prefix ? fopen((string(opt.log_prefix) + "_t_del.tsv").c_str(), "w") :nullptr;
	if (fo) fprintf(fo, "%s\t%s\t%s\t%s\n", "ID", "beg", "end", "len");
	for (int i: t_del) {
		const RepInterval &t = agg_t[i];
		if (fo) fprintf(fo, "%d\t%d\t%d\t%d\n", i, t.beg, t.end, t.end - t.beg);
		t_del_len += t.end - t.beg;
	}
	if (fo) fclose(fo);
	fprintf(stderr, "    target deletion length: %d, percentage: %.2f %%\n", t_del_len, 100.0 * t_del_len / t_len);

	// Query deletion
	fo = opt.log_prefix ? fopen((string(opt.log_prefix) + "_q_del.tsv").c_str(), "w") :nullptr;
	if (fo) fprintf(fo, "%s\t%s\t%s\t%s\n", "ID", "beg", "end", "len");
	for (int j: q_del) {
		const RepInterval &q = agg_q[j];
		if (fo) fprintf(fo, "%d\t%d\t%d\t%d\n", j, q.beg, q.end, q.end - q.beg);
		q_del_len += q.end - q.beg;
	}
	if (fo) fclose(fo);
	fprintf(stderr, "    query deletion length: %d, percentage: %.2f %%\n", q_del_len, 100.0 * q_del_len / q_len);

	fprintf(stderr, "Construct topology: %.2f real time, %.2f CPU time\n", realtime() - t_real, cputime() - t_cpu);

	// Finalizing the alignment
	vector<RepInterval> aln_t, aln_q;
	vector<RepInterval> del_t, del_q;
	for (auto &pair: aln) {
		aln_t.push_back(agg_t[pair.first]);
		aln_q.push_back(agg_q[pair.second]);
	}
	for (int i: t_del) {
		del_t.push_back(agg_t[i]);
	}
	for (int i: q_del) {
		del_q.push_back(agg_q[i]);
	}

	{
		for (int i = 1; i < aln_t.size(); i++) {
			assert(aln_t[i].beg >= aln_t[i-1].end);
			assert(aln_q[i].beg >= aln_q[i-1].end);
		}
		for (int i = 1; i < del_t.size(); i++) {
			assert(del_t[i].end >= del_t[i-1].beg);
		}
		for (int i = 1; i < del_q.size(); i++) {
			assert(del_q[i].end >= del_q[i-1].beg);
		}
	}

	string final_ct, final_cq;
	int jt = 0, jq = 0;
	int long_gap_cnt = 0;
	for (int i = 0; i <= aln_q.size(); i++) {
		int t_beg = i == 0 ?0 :aln_t[i-1].end, t_end = i == aln_q.size() ?t_len :aln_t[i].beg;
		int q_beg = i == 0 ?0 :aln_q[i-1].end, q_end = i == aln_q.size() ?q_len :aln_q[i].beg;
		int gap_t_len = t_end - t_beg, gap_q_len = q_end - q_beg;
		// fprintf(stderr, "[%d,%d) t_len=%d -> [%d,%d) q_len=%d\n", t_beg, t_end, len_t, q_beg, q_end, len_q);

		// Find deleted segments within the gap
		vector<RepInterval> sub_t;
		for (; jt < del_t.size(); jt++) {
			int b = del_t[jt].beg, e = del_t[jt].end;
			if (b >= t_beg and e <= t_end) {
				RepInterval x = del_t[jt];
				x.beg -= t_beg;
				x.end -= t_beg;
				sub_t.push_back(x);
			} else if (b >= t_end) {
				break;
			}
		}
		vector<RepInterval> sub_q;
		for (; jq < del_q.size(); jq++) {
			int b = del_q[jq].beg, e = del_q[jq].end;
			if (b >= q_beg and e <= q_end) {
				RepInterval x = del_q[jq];
				x.beg -= q_beg;
				x.end -= q_beg;
				sub_q.push_back(x);
			} else if (b >= q_end) {
				break;
			}
		}
		// Caution: deletion-insertion runs might be caused by inaccurate partitioning

		// Collect the left segments
		string left_t;
		vector<int> map_t;
		const char *st = t_seq + t_beg;
		for (int j = 0; j <= sub_t.size(); j++) {
			int b = (j == 0) ?0 :sub_t[j-1].end;
			int e = (j == sub_t.size()) ?gap_t_len :sub_t[j].beg;
			for (int k = b; k < e; k++) {
				left_t.push_back(st[k]);
				map_t.push_back(k);
			}
		}
		left_t.push_back('$'); // Guard
		map_t.push_back(gap_t_len);

		string left_q;
		vector<int> map_q;
		const char *sq = q_seq + q_beg;
		for (int j = 0; j <= sub_q.size(); j++) {
			int b = (j == 0) ?0 :sub_q[j-1].end;
			int e = (j == sub_q.size()) ?gap_q_len :sub_q[j].beg;
			for (int k = b; k < e; k++) {
				left_q.push_back(sq[k]);
				map_q.push_back(k);
			}
		}
		left_q.push_back('$');
		map_q.push_back(gap_q_len);

		// Be careful with long gaps
		if (left_t.length() > 5000 or left_q.length() > 5000) {
			long_gap_cnt++;
		}

		// FIXME: use SI scoring matrix
		AlnSta as = global_cigar(opt, left_t.length(), left_t.data(), left_q.length(), left_q.data());
		const string &ext_t = as.ext_a;
		const string &ext_q = as.ext_b;

		int last_t = 0, last_q = 0;
		int pnt_t = 0, pnt_q = 0;
		assert(final_ct.size() == final_cq.size());
		int kt = 0, kq = 0;
		for (int j = 0; j < ext_t.length(); j++) {
			if (ext_t[j] != '-') {
				assert(ext_t[j] == left_t[pnt_t]);
				int pos_t = map_t[pnt_t];
				if (ext_t[j] != '$') assert(ext_t[j] == st[pos_t]);
				// Duplication deletions between bases
				for (; kt < sub_t.size(); kt++) {
					int b = sub_t[kt].beg, e = sub_t[kt].end;
					if (b >= last_t and e <= pos_t) {
						final_ct.push_back('[');
						final_cq.push_back('[');
						for (int c = b; c < e; c++) {
							final_ct.push_back(st[c]);
							final_cq.push_back('-');
						}
						final_ct.push_back(']');
						final_cq.push_back(']');
					} else if (b >= pos_t) break;
				}
				last_t = pos_t;
				pnt_t++;
			}
			if (ext_q[j] != '-') {
				assert(ext_q[j] == left_q[pnt_q]);
				int pos_q = map_q[pnt_q];
				if (ext_q[j] != '$') assert(ext_q[j] == sq[pos_q]);
				for (; kq < sub_q.size(); kq++) {
					int b = sub_q[kq].beg, e = sub_q[kq].end;
					if (b >= last_q and e <= pos_q) {
						final_ct.push_back('[');
						final_cq.push_back('[');
						for (int c = b; c < e; c++) {
							final_ct.push_back('-');
							final_cq.push_back(sq[c]);
						}
						final_ct.push_back(']');
						final_cq.push_back(']');
					} else if (b >= pos_q) break;
				}
				last_q = pos_q;
				pnt_q++;
			}
			if (ext_t[j] != '$') final_ct.push_back(ext_t[j]);
			if (ext_q[j] != '$') final_cq.push_back(ext_q[j]);
		}
		assert(final_cq.size() == final_ct.size());

		if (i < aln_q.size()) {
			as = global_cigar(opt, aln_t[i].end - aln_t[i].beg, t_seq + aln_t[i].beg,
				aln_q[i].end - aln_q[i].beg, q_seq + aln_q[i].beg);
			final_ct.push_back('[');
			final_ct += as.ext_a;
			final_ct.push_back(']');
			final_cq.push_back('[');
			final_cq += as.ext_b;
			final_cq.push_back(']');
		}
	}
	fprintf(stderr, "Found %d long gaps\n", long_gap_cnt);

	// Sanity check
	int i = 0;
	for (char c: final_ct) {
		if (c == '[' or c == ']' or c == '-') {
			continue;
		}
		assert(c == t_seq[i++]);
	}
	assert(i == t_len);

	i = 0;
	for (char c: final_cq) {
		if (c == '[' or c == ']' or c == '-') {
			continue;
		}
		assert(c == q_seq[i++]);
	}
	assert(i == q_len);

	// fprintf(stdout, "%s\n", final_ct.data());
	// fprintf(stdout, "%s\n", final_cq.data());

	int cnt_mat = 0, cnt_mis = 0, cnt_del = 0, cnt_ins = 0;
	for (i = 0; i < final_ct.length(); i++) {
		if (final_ct[i] == '-') cnt_ins++;
		else if (final_cq[i] == '-') cnt_del++;
		else if (final_cq[i] == final_ct[i]) cnt_mat++;
		else cnt_mis++;
	}

	fprintf(stderr, "Finalized alignment: ");
	int cnt_aln = cnt_mis + cnt_mat;
	fprintf(stderr, "    aligned length: %d, map ratio: %.2f %%\n", cnt_aln, 100.0 * cnt_aln / min(t_len, q_len));
	fprintf(stderr, "    mismatch: %d, percentage: %.2f %%\n", cnt_mis, 100.0 * cnt_mis / cnt_aln);
	fprintf(stderr, "    deletion: %d, percentage: %.2f %%\n", cnt_del, 100.0 * cnt_del / t_len);
	fprintf(stderr, "    insertion: %d, percentage: %.2f %%\n", cnt_ins, 100.0 * cnt_ins / q_len);
	fprintf(stderr, "\n");

	// Output CIGAR
	extended_paf_format(name1, t_len, final_ct, name2, q_len, final_cq);
}

int usage(const ZigOptions &o) {
	fprintf(stderr, "Usage: zigalign [options] seq1.fa seq2.fa > aln.paf\n");
	fprintf(stderr, "  Common options:\n");
	fprintf(stderr, "    -t [INT]  number of threads\n");
	fprintf(stderr, "    -v [STR]  intermediate results prefix\n");
	fprintf(stderr, "  Scoring parameters for pairwise alignment:\n");
	fprintf(stderr, "    -A [INT]  match score [%d]\n", o.mat_score);
	fprintf(stderr, "    -B [INT]  mismatch penalty [%d]\n", o.mis_pen);
	fprintf(stderr, "    -O [INT]  open gap(indel) penalty [%d]\n", o.gap_o);
	fprintf(stderr, "    -E [INT]  extend gap penalty [%d]\n", o.gap_e);
	fprintf(stderr, "    -D [INT]  repeat unit deletion penalty [%d]\n", o.del_unit);
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
	while ((c = getopt(argc, argv, "A:B:O:E:D:u:d:p:a:b:o:e:k:v:t:")) >= 0) {
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
				opt.del_unit = -abs(str2int(optarg));
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
			default:
				fprintf(stderr, "Unrecognized option `%c`\n", c);
				return 1;
		}
	}
	omp_set_num_threads(opt.n_threads);

	if (argc - optind == 2) {
		// align_with_dups(opt, argv[optind], argv[optind+1]);
		align_long_seq(opt, argv[optind], argv[optind+1]);
	} else {
		fprintf(stderr, "Two FASTA files are required\n");
		return 1;
	}
	fprintf(stderr, "Program finishes in %.3f CPU seconds, %.3f real seconds\n", cputime()-ctime, realtime()-rtime);
	return 0;
}