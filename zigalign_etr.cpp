#include <cmath>
#include <memory>
#include <sstream>

#include "zigalign.h"
#include "wraparound.h"

vector<RepInterval> self_alignment(const ZigOptions &opt, int n, const char *seq, vector<vector<Bt1Cell>> &bt)
{
	if (n > SA_MAX_LEN) {
		fprintf(stderr, "Sequence length %d exceeds the limit of self alignment %d\n", n, SA_MAX_LEN);
		abort();
	}

	const int MIN_UNIT = opt.min_unit_size;
	const int MAX_UNIT = opt.max_unit_size;
	const int MAT_SCORE = 0;
	const int OPEN_TR = opt.open_tr_pen;
	const int CLOSE_TR = opt.close_tr_pen;
	const int SA_MAT_SCORE = opt.self_aln.mat;
	const int SA_MIS_PEN = opt.self_aln.mis;
	const int SA_GAP_O = opt.self_aln.gap_o;
	const int SA_GAP_E = opt.self_aln.gap_e;

	const int BT_NON_REP = 0;
	const int BT_START_REP = 1;
	const int BT_NEW_COPY = 2;
	const int BT_WITHIN_REP = 3;
	const int BT_END_REP = 4;

	if (MAT_SCORE >= SA_MAT_SCORE) {
		fprintf(stderr, "Warning: setting higher score for diagonal match can't detect tandem repeats\n");
	}

	if (n < MIN_UNIT) {
		vector<RepInterval> ret;
		return ret;
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
		bt[i][i].event = BT_NON_REP;
		bt[i][i].pi = 1;
		bt[i][i].pj = i-1;
	}
	prev_dp[MIN_UNIT].H = MIN_UNIT * MAT_SCORE;
	prev_dp[MIN_UNIT].beg = 1;

	// Main loop
	for (int i = MIN_UNIT + 1; i <= n; i++) {
		for (int j = 0; j <= i; j++) curr_dp[j] = Dp1Cell();
		curr_dp[i].H = prev_dp[i-1].H + MAT_SCORE;
		bt[i][i].event = BT_NON_REP;
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

struct TandemRepeats {
	int pid; // Pattern ID
	int base_idx; // Base index
	string pattern;
	double ave_len;
	double ave_mis;
	double ave_gap;
	vector<RepInterval> repeats;
};

vector<TandemRepeats> group_repeats(const vector<RepInterval> &reps)
{
	vector<TandemRepeats> ret;
	if (reps.empty()) return ret;

	TandemRepeats tr;
	tr.ave_len = reps[0].end - reps[0].beg;
	tr.ave_mis = tr.ave_gap = 0;
	tr.repeats.push_back(reps[0]);
	for (int i = 1; i < reps.size(); i++) {
		const RepInterval &r = reps[i];
		if (r.mat == -1) { // New group
			tr.ave_len /= tr.repeats.size();
			tr.ave_mis /= (tr.repeats.size() - 1);
			tr.ave_gap /= (tr.repeats.size() - 1);
			ret.push_back(tr);
			tr.ave_len = r.end - r.beg;
			tr.ave_mis = tr.ave_gap = 0;
			tr.repeats.push_back(r);
		} else {
			tr.ave_len += r.end - r.beg;
			tr.ave_gap += r.gap;
			tr.ave_mis += r.mis;
			tr.repeats.push_back(r);
		}
	}
	tr.ave_len /= tr.repeats.size();
	tr.ave_mis /= (tr.repeats.size() - 1);
	tr.ave_gap /= (tr.repeats.size() - 1);
	ret.push_back(tr);
	return ret;
}

void write_repeats(FILE *fo, const vector<TandemRepeats> &tr)
{
	fprintf(fo, "%s\t%s\t%s\t%s\t%s\t%s\n", "ID", "beg", "end", "len", "mis", "gap");
	for (int i = 0; i < tr.size(); i++) {
		for (const RepInterval &r: tr[i].repeats) {
			fprintf(fo, "%d\t%d\t%d\t%d\t%d\t%d\n", i + 1, r.beg, r.end, r.end - r.beg, r.mis, r.gap);
		}
	}
}

vector<TandemRepeats> read_repeats(const string &fn)
{
	vector<TandemRepeats> ret;
	ifstream in(fn);
	assert(in.is_open());
	string line;
	TandemRepeats tr;
	int last_id = -1;
	getline(in, line); // Skip metaline
	while (getline(in, line)) {
		int id, len;
		stringstream ss(line);
		RepInterval r;
		ss >> id >> r.beg >> r.end >> len >> r.mis >> r.gap;
		if (last_id != -1 and id != last_id) {
			ret.push_back(tr);
			tr.repeats.clear();
		}
		tr.repeats.push_back(r);
		last_id = id;
	}
	ret.push_back(tr);
	in.close();

	for (int i = 0; i < ret.size();i ++) {
		TandemRepeats &tr = ret[i];
		tr.pid = -1;
		tr.base_idx = 0;
		tr.ave_len = tr.ave_gap = tr.ave_mis = 0;
		int n_reps = 0;
		for (const RepInterval &r: tr.repeats) {
			tr.ave_len += r.end - r.beg;
			if (r.mis != -1) {
				n_reps++;
				tr.ave_gap += r.gap;
				tr.ave_mis += r.mis;
			}
		}
		tr.ave_len /= tr.repeats.size();
		if (n_reps) tr.ave_mis /= n_reps;
		if (n_reps) tr.ave_gap /= n_reps;
		fprintf(stderr, "TR Group %d, reps=%ld, len=%.2f, mis=%.2f, gap=%.2f, range=[%d,%d)\n",
			i + 1, tr.repeats.size(), tr.ave_len, tr.ave_mis, tr.ave_gap, tr.repeats.front().beg, tr.repeats.back().end);
	}
	fprintf(stderr, "\n");
	return ret;
}

// Test if group B can be merged into group A
bool merge_tr_test(const TandemRepeats &tr_a, const char *seq_a, const TandemRepeats &tr_b, const char *seq_b)
{
	const int MAX_COMPARE = 300;
	// The average length is more reliable
	double len_diff = abs(tr_a.ave_len - tr_b.ave_len) / tr_a.ave_len;
	if (len_diff < MAX_TR_DIV) {
		// Template of B: U + U
		const RepInterval &b = tr_b.repeats.front();
		int len_b = b.end - b.beg;
		string sb;
		sb.resize(len_b);
		memcpy((char*) sb.data(), seq_b + b.beg, len_b);
		sb += sb;
		// Use average divergence to reduce disturbance of outliers
		int n_comp = 0;
		double ave_div = 0;
		for (const RepInterval &a: tr_a.repeats) {
			int len_a = a.end - a.beg;
			double len_div = 1.0 * abs(len_b - len_a) / min(len_b, len_a);
			if (len_div > MAX_TR_DIV) continue;
			// Use semi-global alignment to handle repeat rotations
			SgResult sg = semi_global(len_a, seq_a + a.beg, sb.length(), sb.data());
			n_comp++;
			ave_div += 1.0 * (sg.mis + sg.gap) / len_a;
			if (n_comp > MAX_COMPARE) break;
		}
		ave_div /= n_comp;
		if (ave_div < MAX_TR_DIV) {
			return true;
		}
		// fprintf(stderr, "ave_div = %.2f, merge = %d\n", ave_div, to_merge);
	}
	return false;
}

vector<TandemRepeats> collect_long_repeats2(const ZigOptions &opt, int n, const char *seq)
{
	int n_threads = opt.n_threads;
	vector<vector<Bt1Cell>> bt_bin[n_threads];
	vector<RepInterval> reps_bin[n_threads];
	vector<TandemRepeats> tr_bin[n_threads];
	vector<TandemRepeats> ret;
	int global_os = 0;
	int batch_id = 0;
	while (global_os < n) {
		double r_start = realtime(), c_start = cputime();
		int m = n - global_os;
		const char *s0 = seq + global_os;

		// Multi-threading DP
		#pragma omp parallel for
		for (int i = 0; i < n_threads; i++) {
			int beg = i * SA_PART_LEN;
			int end = min(m, beg + SA_PART_LEN);
			if (beg >= end) {
				reps_bin[i].clear();
				tr_bin[i].clear();
			} else {
				const char *s = s0 + beg;
				vector<RepInterval> curr_reps, next_reps;
				curr_reps = self_alignment(opt, end - beg, s, bt_bin[i]);
				reps_bin[i].clear();
				while (not curr_reps.empty()) {
					next_reps.clear();
					for (const RepInterval &r: curr_reps) {
						// fprintf(stderr, "Break down [%d, %d) len=%d\n", r.beg, r.end, r.end - r.beg);
						vector<RepInterval> tmp = self_alignment(opt, r.end - r.beg, s + r.beg, bt_bin[i]);
						if (tmp.empty()) {
							// Non-nested
							reps_bin[i].push_back(r);
						} else {
							// Recursively break down the nested repeat
							for (RepInterval &t: tmp) {
								t.beg += r.beg;
								t.end += r.beg;
								next_reps.push_back(t);
							}
						}
					}
					swap(curr_reps, next_reps);
				}
				sort(reps_bin[i].begin(), reps_bin[i].end());
				// Remove the last unit because it might be incomplete
				if (not reps_bin[i].empty()) {
					reps_bin[i].resize(reps_bin[i].size() - 1);
				}
				for (RepInterval &r: reps_bin[i]) {
					r.beg += global_os + beg;
					r.end += global_os + beg;
				}
				tr_bin[i] = group_repeats(reps_bin[i]);
			}
		}
		double c_time = cputime() - c_start, r_time = realtime() - r_start;
		int n_reps = 0;
		for (int i = 0; i < n_threads; i++) {
			n_reps += reps_bin[i].size();
		}
		fprintf(stderr, "Batch %d, offset: %d, reps: %d; CPU: %.2f, Wall: %.2f, Core: %.2f\n",
			++batch_id, global_os, n_reps, c_time, r_time, c_time / r_time);

		if (0) {
			for (int i = 0; i < n_threads; i++) {
				fprintf(stdout, "thread %d\n", i);
				for (const TandemRepeats &g: tr_bin[i]) {
					fprintf(stdout, "ave len=%.2f, mis=%.2f, gap=%.2f, range=[%d,%d)\n",
						g.ave_len, g.ave_mis, g.ave_gap, g.repeats.front().beg, g.repeats.back().end);
					// for (const RepInterval &r: g.repeats) {
					// 	fprintf(stdout, "[%d, %d) len=%d, mat=%d, mis=%d, gap=%d\n",
					// 		r.beg, r.end, r.end - r.beg, r.mat, r.mis, r.gap);
					// }
					// fprintf(stdout, "\n");
				}
			}
		}

		// Bridge adjacent batches
		for (int i = 0; i < n_threads; i++) {
			for (const TandemRepeats &g: tr_bin[i]) {
				if (ret.empty()) {
					ret.push_back(g);
				} else {
					TandemRepeats &t = ret.back();
					bool to_merge = merge_tr_test(t, seq, g, seq);
					if (to_merge) {
						t.repeats.insert(t.repeats.end(), g.repeats.begin(), g.repeats.end());
					} else {
						ret.push_back(g);
					}
				}
			}
		}
		global_os += n_threads * SA_PART_LEN;
	}

	for (int i = 0; i < ret.size();i ++) {
		// Recalculate meta information for the merged tandem repeats
		TandemRepeats &tr = ret[i];
		tr.pid = -1;
		tr.base_idx = 0;
		tr.ave_len = tr.ave_gap = tr.ave_mis = 0;
		int n_reps = 0;
		for (const RepInterval &r: tr.repeats) {
			tr.ave_len += r.end - r.beg;
			if (r.mis != -1) {
				n_reps++;
				tr.ave_gap += r.gap;
				tr.ave_mis += r.mis;
			}
		}
		tr.ave_len /= tr.repeats.size();
		if (n_reps) tr.ave_mis /= n_reps;
		if (n_reps) tr.ave_gap /= n_reps;
		fprintf(stderr, "TR Group %d, reps=%ld, len=%.2f, mis=%.2f, gap=%.2f, range=[%d,%d)\n",
			i + 1, tr.repeats.size(), tr.ave_len, tr.ave_mis, tr.ave_gap, tr.repeats.front().beg, tr.repeats.back().end);
	}

	return ret;
}

vector<TandemRepeats> collect_long_repeats(const ZigOptions &opt, int n, const char *seq)
{
	assert(opt.n_threads > 0);
	int global_os = 0;
	vector<RepInterval> reps_bin[opt.n_threads];
	vector<vector<Bt1Cell>> bt_mat;
	string pattern;
	vector<RepInterval> all_reps;
	vector<TandemRepeats> ret;

	// TODO: use thread pipeline to achieve better parallelization
	int batch_id = 0;
	while (global_os < n) {
		double r_start = realtime(), c_start = cputime();
		int m = n - global_os;
		const char *s0 = seq + global_os;
		if (pattern.empty()) {
			// Discover repeat pattern
			int k = min(SA_PART_LEN, m);
			vector<RepInterval> pool = self_alignment(opt, k, s0, bt_mat); // Multiple pattern?
			vector<RepInterval> reps;
			// Unwind nested repeats recursively to discover as short units as possible
			while (not pool.empty()) {
				vector<RepInterval> next_pool;
				for (const RepInterval &r: pool) {
					if (r.end - r.beg < 1000) {
						reps.push_back(r);
						continue;
					}
					vector<RepInterval> tmp = self_alignment(opt, r.end - r.beg, s0 + r.beg, bt_mat);
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
		int n_threads = min((m + SA_PART_LEN - 1) / SA_PART_LEN, opt.n_threads);
		#pragma omp parallel for num_threads(n_threads)
		for (int i = 0; i < n_threads; i++) {
			reps_bin[i].clear();
			int os = i * SA_PART_LEN;
			int t_len = min(m - os, SA_PART_LEN);
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
			TandemRepeats lr;
			lr.pattern = pattern;
			lr.repeats = all_reps;
			ret.push_back(lr);
			pattern.clear();
			all_reps.clear();
		}

		// TODO: post-process the poorly aligned regions
	}
	if (not pattern.empty()) {
		TandemRepeats lr;
		lr.pattern = pattern;
		lr.repeats = all_reps;
		ret.push_back(lr);
	}

	fprintf(stderr, "Found %ld patterns\n", ret.size());
	for (int i = 0; i < ret.size(); i++) {
		const TandemRepeats &r = ret[i];
		const string &s = r.pattern;
		int beg = r.repeats.front().beg;
		int end = r.repeats.back().end;
		int len = end - beg;
		fprintf(stderr, "Pattern %d: len=%ld, range=[%d,%d), len=%d\n", i+1, s.length(), beg, end, len);
	}
	return ret;
}

TandemRepeats normalize_repeats(const ZigOptions &opt, int pat_len, const char *pat_seq, int n, const char *seq)
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
			int os = i * SA_PART_LEN;
			int t_len = min(m - os, SA_PART_LEN);
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
	TandemRepeats ret;
	ret.pattern = rep_pat;
	ret.repeats = all_reps;
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

void extended_paf_format(const string &t_name, int t_len, const string &ext_t,
	const string &q_name, int q_len, const string ext_q, const vector<int> &rid_map)
{
	assert(ext_t.length() == ext_q.length());
	char strand = '+';
	int matches_n = 0, mismatches_n = 0, len = ext_t.length();
	int mapq = 60;

	string cigar;
	char op_type = '@';
	int op_cnt = 0;
	int t_pointer = 0;
	for (int i = 0; i < len; i++) {
		if (ext_t[i] == 'A' or ext_t[i] == 'C' or ext_t[i] == 'G' or ext_t[i] == 'T') {
			t_pointer++;
		}
		if (ext_t[i] == '[' or ext_t[i] == ']') {
			if (op_type != '@') {
				cigar += to_string(op_cnt);
				cigar += op_type;
			}
			if (ext_t[i] == '[') {
				int id = -1;
				for (int j = 0; j < 20; j++) {
					if (t_pointer - j >= 0 and rid_map[t_pointer - j] != -1) {
						id = rid_map[t_pointer - j];
						break;
					}
					if (t_pointer + j < t_len and rid_map[t_pointer + j] != -1) {
						id = rid_map[t_pointer + j];
						break;
					}
				}
				if (id != -1) cigar += '(' + to_string(id) + ')';
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

void align_etr2(const ZigOptions &opt, const char *fn_t, const char *fn_q)
{
	FastaSeq target = input_fasta_seq(fn_t);
	FastaSeq query = input_fasta_seq(fn_q);
	const string &name_t = target.name, &name_q = query.name;
	const string &seq_t = target.seq, &seq_q = query.seq;
	int t_len = seq_t.length(), q_len = seq_q.length();
	const char *t_seq = seq_t.data(), *q_seq = seq_q.data();

	vector<TandemRepeats> lr_t = collect_long_repeats2(opt, t_len, t_seq);
	vector<TandemRepeats> lr_q = collect_long_repeats2(opt, q_len, q_seq);

	if (opt.log_prefix) {
		FILE *fo = fopen((string(opt.log_prefix) + "_t_reps.tsv").c_str(), "w");
		assert(fo);
		write_repeats(fo, lr_t);
		fclose(fo);

		fo = fopen((string(opt.log_prefix) + "_q_reps.tsv").c_str(), "w");
		assert(fo);
		write_repeats(fo, lr_q);
		fclose(fo);
	}

	int t_sum = 0, q_sum = 0;
	for (const TandemRepeats &lr: lr_t) {
		for (const RepInterval &r: lr.repeats) {
			t_sum += r.end - r.beg;
		}
	}
	for (const TandemRepeats &lr: lr_q) {
		for (const RepInterval &r: lr.repeats) {
			q_sum += r.end - r.beg;
		}
	}
	fprintf(stderr, "Target repeat fraction[%%]: %.2f\n", 100.0 * t_sum / t_len);
	fprintf(stderr, "Query repeat fraction[%%]: %.2f\n", 100.0 * q_sum / q_len);
	exit(1);

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
			TandemRepeats lr;
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
	for (TandemRepeats &lr: lr_t) {
		lr.base_idx = n;
		n += lr.repeats.size();
		for (const RepInterval &r: lr.repeats) {
			t_sum += r.end - r.beg;
		}
	}
	for (TandemRepeats &lr: lr_q) {
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
	for (const TandemRepeats &lt : lr_t) {
		int idx_t = lt.base_idx;
		int pid_t = lt.pid;
		for (const TandemRepeats &lq : lr_q) {
			int idx_q = lq.base_idx;
			int pid_q = lq.pid;
			if (pid_q != pid_t) continue;

			#pragma omp parallel for
			for (int i = 0; i < lt.repeats.size(); i++) {
				const RepInterval &t = lt.repeats[i];
				for (int j = 0; j < lq.repeats.size(); j++) {
					const RepInterval &q = lq.repeats[j];
					matrix[idx_t + i][idx_q + j] = global_alignment(
						opt.self_aln, t.end - t.beg, t_seq + t.beg, q.end - q.beg, q_seq + q.beg);
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
	const int DEL_UNIT = opt.del_unit_pen;
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
	for (const TandemRepeats &lr: lr_t) {
		agg_t.insert(agg_t.end(), lr.repeats.begin(), lr.repeats.end());
	}
	vector<RepInterval> agg_q;
	for (const TandemRepeats &lr: lr_q) {
		agg_q.insert(agg_q.end(), lr.repeats.begin(), lr.repeats.end());
	}

	vector<int> rep_id_map(t_len + 1, -1);
	for (const auto &lr :lr_t) {
		for (const auto &r: lr.repeats) {
			rep_id_map[r.beg] = lr.pid; // Mark the position one base before the repeat
		}
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
	extended_paf_format(name_t, t_len, final_ct, name_q, q_len, final_cq, rep_id_map);
}

void categorize_repeats(vector<TandemRepeats> &lr_t, const char *seq_t, vector<TandemRepeats> &lr_q, const char *seq_q)
{
	int pid = 0;
	for (TandemRepeats &tr: lr_t) {
		tr.pid = pid++;
	}
	for (TandemRepeats &tr: lr_q) {
		tr.pid = pid++;
	}
	vector<int> set_id(pid);
	for (int i = 0; i < pid; i++) {
		set_id[i] = i;
	}

	for (int i = 0; i < lr_t.size(); i++) {
		TandemRepeats &a = lr_t[i];
		for (int j = i - 1; j >= 0; j--) {
			TandemRepeats &b = lr_t[j];
			bool to_merge = merge_tr_test(b, seq_t, a, seq_t);
			if (to_merge) {
				us_union(set_id, a.pid, b.pid);
				a.pid = us_find(set_id, a.pid);
				b.pid = us_find(set_id, b.pid);
			}
		}
	}
	for (int i = 0; i < lr_q.size(); i++) {
		TandemRepeats &a = lr_q[i];
		for (int j = i - 1; j >= 0; j--) {
			TandemRepeats &b = lr_q[j];
			bool to_merge = merge_tr_test(b, seq_q, a, seq_q);
			if (to_merge) {
				us_union(set_id, a.pid, b.pid);
				a.pid = us_find(set_id, a.pid);
				b.pid = us_find(set_id, b.pid);
			}
		}
	}

	for (int i = 0; i < lr_t.size(); i++) {
		TandemRepeats &a = lr_t[i];
		for (int j = 0; j < lr_q.size(); j++) {
			TandemRepeats &b = lr_q[j];
			bool to_merge = merge_tr_test(b, seq_q, a, seq_t);
			if (to_merge) {
				us_union(set_id, a.pid, b.pid);
				a.pid = us_find(set_id, a.pid);
				b.pid = us_find(set_id, b.pid);
			}
		}
	}
	for (int i = 0; i < lr_q.size(); i++) {
		TandemRepeats &a = lr_q[i];
		for (int j = 0; j < lr_t.size(); j++) {
			TandemRepeats &b = lr_t[j];
			bool to_merge = merge_tr_test(b, seq_t, a, seq_q);
			if (to_merge) {
				us_union(set_id, a.pid, b.pid);
				a.pid = us_find(set_id, a.pid);
				b.pid = us_find(set_id, b.pid);
			}
		}
	}

	for (int i = 0; i < lr_t.size(); i++) {
		const TandemRepeats &tr = lr_t[i];
		fprintf(stderr, "TR Group %d, reps=%ld, len=%.2f, mis=%.2f, gap=%.2f, range=[%d,%d)\n",
			tr.pid, tr.repeats.size(), tr.ave_len, tr.ave_mis, tr.ave_gap, tr.repeats.front().beg, tr.repeats.back().end);
	}
	fprintf(stderr, "\n");
	for (int i = 0; i < lr_q.size(); i++) {
		const TandemRepeats &tr = lr_q[i];
		fprintf(stderr, "TR Group %d, reps=%ld, len=%.2f, mis=%.2f, gap=%.2f, range=[%d,%d)\n",
			tr.pid, tr.repeats.size(), tr.ave_len, tr.ave_mis, tr.ave_gap, tr.repeats.front().beg, tr.repeats.back().end);
	}
}

void align_etr(const ZigOptions &opt, const char *fn_t, const char *fn_q)
{
	FastaSeq target = input_fasta_seq(fn_t);
	FastaSeq query = input_fasta_seq(fn_q);
	const string &name_t = target.name, &name_q = query.name;
	const string &seq_t = target.seq, &seq_q = query.seq;
	int t_len = seq_t.length(), q_len = seq_q.length();
	const char *t_seq = seq_t.data(), *q_seq = seq_q.data();

	vector<TandemRepeats> lr_t = read_repeats(string(opt.log_prefix) + "_t_reps.tsv");
	vector<TandemRepeats> lr_q = read_repeats(string(opt.log_prefix) + "_q_reps.tsv");

	categorize_repeats(lr_t, t_seq, lr_q, q_seq);
}
